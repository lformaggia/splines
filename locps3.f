      subroutine locps3(n,m,keypa,aps,xp,ug,vg,u10,u20,v10,v20,
     1                  r,z,dist,imes,iter,tole,distl)
      parameter(itmax=100,zerom=1.e-6,tol=zerom,eps=zerom)
c
      real xp(3),r(3),z(3),xi(2,2),g(2)
      real ru(3),rv(3),xinit(2,2,2)
      real xit(2),p(2),pt(2),ptt(2)
      integer keypa(*)
      real aps(3,*)
c
      data xinit/1.,0.,0.,1.,
     1           0.707106,-0.707106,-0.707106,+0.707106/
c
      ipa(a) = max(1,min(int(a),n-1))
      jpa(a) = max(1,min(int(a),m-1))
      ugl(q) = max(u1,min(q,u2))
      vgl(q) = max(v1,min(q,v2))
c
      imes  = 0
      disto = distl
      toler = tole
c
c     ************ limit bounds  ***********
c
      u1 = u10
      v1 = v10
      u2 = u20
      v2 = v20
      ug = min(u2,max(u1,ug))
      vg = min(v2,max(v1,vg))
      i  = ipa(ug)
      j  = jpa(vg)
      u  = ug - float(i)
      v  = vg - float(j)
c
c set iteration counters
c
      iter   = 0
      itry   = 1
c
c check distance
c
      npa    = (i-1)*(m-1)+j
      kpa    = keypa(npa)
      ndu    = aps(1,kpa)
      ndv    = aps(2,kpa)
c      xl     = aps(3,kpa)
      xl     = 1.0
      p(1)   = ug
      p(2)   = vg
      xit(1) = 0.
      xit(2) = 0.
      fret = distps(xp,ug,vg,xit,0.,keypa,aps,n,m,z,
     1               r,der,g,ugn,vgn,0)
      if(fret.lt.disto*xl)then
         goto 1000
      endif
c
987   continue
c
c     ************* initialise ***************
c
c
      p(1)    = ug
      p(2)    = vg
c
c set initial directions
c
      xi(1,1) = xinit(1,1,itry)
      xi(2,1) = xinit(2,1,itry)
      xi(1,2) = xinit(1,2,itry)
      xi(2,2) = xinit(2,2,itry)
      pt(1)   = p(1)
      pt(2)   = p(2)
c
c       ****************** main loop ***************
c
1      iter = iter +1
c
c sweep over direction set
c
       fp = fret
       ibig =0
       del =0.
       do 13 i=1,2
         xit(1) = xi(1,i)
         xit(2) = xi(2,i)
         ugs   = p(1)
         vgs   = p(2)
         fret0 = xminli(tol,xmin,xp,r,z,ugs,vgs,xit,ru,rv,
     1                 der,g,u1,v1,u2,v2,ugn,vgn,n,m,
     1                 keypa,aps,xl, ierr)
         p(1) = ugn
         p(2) = vgn
c
c check distance
c
      if(fret0.lt.disto*xl)then
         fret = fret0
         goto 1000
      endif
c
c get the largest descent direction
c
         if(abs(fret-fret0).gt.del)then
           del  = abs(fret-fret0)
           ibig = i
         endif
         fret = fret0
13    continue
c
c          ********** convergence tests *************
c
      if(fret.lt.disto*xl)then
c
C test   distance
c
        goto 1000
      endif
      if(max(abs(p(1)-pt(1)),abs(p(2)-pt(2))).lt.eps*max(ug,vg))then
c
c test extrema
c
        imes = -2
        go to 1000
      endif
c
c test convergence
c
      if(2*abs(fp-fret).le.toler*(abs(fp)+abs(fret)))then
         imes = -1
         go to 1000
      endif
c
c test n. of iterations
c
      if(iter.eq.itmax)then
         imes = 2
         go to 1000
      endif
c
c      ************* end convergence tests  **************
c
c construct the extrapolated point and the average direction
c saving  the old starting point into pt
c
      ptt(1) =   2*p(1)-pt(1)
      ptt(2) =   2*p(2)-pt(2)
      xit(1) =   p(1) - pt(1)
      xit(2) =   p(2) - pt(2)
      pt(1)  =      p(1)
      pt(2)  =      p(2)
c
c   ****** check wether largest descent direction should be discarded
c
      ugs = ptt(1)
      vgs = ptt(2)
      fptt   = distps(xp,ugs,vgs,xit,0.,keypa,aps,n,m,z,
     1               r,der,g,ugk,vgk,0)
      if(fptt.ge.fp)goto 1
c                         ** use the old set
c
c
      t = 2.*(fp-2.*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
      if(t.ge.0)goto 1
c                         ** use the old set
c
c              ********** discard ***********
c
c normalize direction
c
      aa = sqrt(xit(1)**2+xit(2)**2)
      aa = 1./max(zerom,aa)
      xit(1) = aa*xit(1)
      xit(2) = aa*xit(2)
c
c move to the minimum along new direction xit
c
      ugs  = p(1)
      vgs  = p(2)
      fret = xminli(tol,xmin,xp,r,z,ugs,vgs,xit,ru,rv,
     1              der,g,u1,v1,u2,v2,ugn,vgn,n,m,
     1              keypa,aps,xl,ierr)
      p(1) = ugn
      p(2) = vgn
      if(fret.lt.disto*xl)then
         goto 1000
      endif
c
c
c discard largest decent direction
c
      xi(1,ibig) = xit(1)
      xi(2,ibig) = xit(2)
      go to 1
c
c normal end conditions
c
1000  continue
c
c try other directions
c
      if(imes.ne.0.and.itry.eq.1.and.iter.lt.8)then
        itry=itry+1
        goto 987
      endif
c
      dist =  fret
      ug   = p(1)
      vg   = p(2)
      return
      end
