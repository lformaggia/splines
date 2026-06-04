C=NAME locps4
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/locps4.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:38 $
C=TYPE subroutine
C=PURPOSE general gensur
C=KEYWORDS Parametric_polynomila_surfaces 2D_Minimization
C=BLOCK ABSTRACT
c
c This routine finds the minimum of the function
c
c          DISTPS = || r - xp ||
c
c where xp is a point on the 3d space and r a point on a parametric surface
c whose coefficents are contained into aps. It makes use of a modified version
c of the Powell discardinding largest descent direction method published on
c Numerical Recipies. The routine carries out a series of line search, on each
c line a minimization algorithm is employed, which makes uses of the same methodology
c adopted for the 1D routine xmin1d. The search is limited to a prescribed interval.
c the starting search point should be near to the desired minimum, so a rough work
c should be done beforehand in order to give a 'good' starting point. The method
c is a 'zero order method', which makes no use of derivative information (slow
c but relatively robust).
C
C=END ABSTRACT
C=BLOCK USAGE
C
c   call locps4(n,m,keypa,aps,xp,ug,vg,u10,u20,v10,v20,r,z,dist,imes,iter,
c               xinit,utl,distl,ftl,zerom)
C
C
c   INPUT     DIMENSION       DESCRIPTION
c
c n                         n. of knots on the surface along u-coordinate
c m                         n. of knots on the surface along v-coordinate
c aps         3,*           surface coefficent in standard format.
c keypa       (n-1)(m-1)    array of pointers to aps.
c xp          3             coordinates of given point.
c ug,vg                     starting value of the global parametric coordinates
c u10,u20                   limit of interval where the minimum is searched (u-coord)
c v10,v20                   limit of interval where the minimum is searched (v-coord)
c xinit       2,2           unit vector on the parameter plane containing the initial
c                           search direction. The two vector must be ORTHOGONAL.
c utl                       tolerance on the variation of the parametric coordinates
c                           between two iteration of the procedure (see NOTES).
c distl                     tolerance on the minimum value of distance (see NOTES).
c ftl                       tolerance on the variation of the value of the distance
c                           between two iteration of the procedure (see NOTES).
c zerom                     machine zero (suggested: 1.e-5 for 32 bit machines)
c
c
c   OUTPUT     DIMENSION       DESCRIPTION
c
C ug,vg                     parametric coordinates u at found minimum point.
C r           3             found minimum r=r(ug,vg).
c dist                      distance ||r(ug,vg) - xp||
c z           3             r(ug,vg) - xp
c imes                      error message (see MESSAGES and DIAGNOSTICS).
c iter                      N. of iterations used for convergence (if imes=/-1)
c
c OTHERS
c
c itmax       parameter     max n. of iterations
c
c
c NOTES:
c
c    during the process the minimum is searched along set of 2 orthogonal directions
c    (in the parameter plane), contained in the array xi. The line search is carried
c    out by using a routine (xminli) similar to xmin1d, but which use a zero order method.
c
c   .........................
c   .Termination conditions .
c   .........................
c
c (1)  max(|u1-u0|,|v1-v0|)<max(|u1|,|v1|)      Variation of parametric coordinates below
c                                            given tolerance
c (2)  |fv - fw | < 0.5*ftl*(|fv|+|fw|)         Variation of the distance between two
c                                               successive iterations below tolerance
c (3)  ||r - xp|| < distl                     Distance below given tolerance
c (4)  iter > itmax                          N. of iterations exceeded
c
c
c
c
c   MESSAGES and DIAGNOSTICS
c
c   imes =0  termination condition (3): normal convergence reached,
c            POINT IS ON the curve: ||r -xp||<distl
c   imes =1  termination condition (2), but ||r-xp||>=distl: normal convergence
c            but point is NOT on the curve
c   imes =2  termination condition (1) but condition (3) is not satisfied:
c            stationary point reached, point is NOT on the curve
c   imes =-1 termination condition (4): ERROR CONDITION: N. OF ITERATION EXCEEDED
c                                       ---------------
c
c
c   NOTES ON THE METHOD
c
c   It uses a modification of a canned routine  which adopts a strategy which
c   makes no use of the derivative information. It is slow but relatively robust.
c   It keeps adjornig two search directions, devined by 2 versors stored in the array xit.
c   The routine xminli carries out the line search. Details on the method may be found
c   on Numerical Recipes.
C
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine locps4(n,m,keypa,aps,xp,ug,vg,u10,u20,v10,v20,
     1                  r,z,dist,imes,iter,xinit,utl,distl,ftl,zerom)
      parameter(itmax=100)
c
c nuova versione: le direzioni iniziali vengono date in input in xinit
c
c
      real xp(3),r(3),z(3),xi(2,2),g(2)
      real ru(3),rv(3),xinit(2,2)
      real xit(2),p(2),pt(2),ptt(2)
      integer keypa(*)
      real aps(3,*)
c
      ipa(a) = max(1,min(int(a),n-1))
      jpa(a) = max(1,min(int(a),m-1))
      ugl(q) = max(u1,min(q,u2))
      vgl(q) = max(v1,min(q,v2))
c
      imes  = 0
      disto = distl
      toler = ftl
      tol   = zerom
      eps   = utl
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
c
c check distance
c
      npa    = (i-1)*(m-1)+j
      kpa    = keypa(npa)
      ndu    = aps(1,kpa)
      ndv    = aps(2,kpa)
      xl     =1.
      p(1)   = ug
      p(2)   = vg
      xit(1) = 0.
      xit(2) = 0.
      fret = distps(xp,ug,vg,xit,0.,keypa,aps,n,m,z,
     1               r,der,g,ru,rv,ugn,vgn,0)
      if(fret.lt.disto)then
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
      xi(1,1) = xinit(1,1)
      xi(2,1) = xinit(2,1)
      xi(1,2) = xinit(1,2)
      xi(2,2) = xinit(2,2)
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
         fret0 = xminli(utl,xmin,xp,r,z,ugs,vgs,xit,ru,rv,
     1                 der,g,u1,v1,u2,v2,ugn,vgn,n,m,
     1                 keypa,aps, ierr,zerom,distl,0)
         p(1) = ugn
         p(2) = vgn
c
c check distance
c
      if(fret0.lt.disto)then
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
      if(fret.lt.disto)then
c
c test   distance
c
        goto 1000
      endif
c
c test convergence
c
      if(2*abs(fp-fret).le.toler*(abs(fp)+abs(fret)))then
         imes = 1
         go to 1000
      endif
c
      if(max(abs(p(1)-pt(1)),abs(p(2)-pt(2))).lt.eps*max(ug,vg))then
c
c test extrema
c
        imes = -2
        go to 1000
      endif
c
c test n. of iterations
c
      if(iter.eq.itmax)then
         imes = -2
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
     1               r,der,g,ru,rv,ugk,vgk,0)
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
      fret = xminli(utl,xmin,xp,r,z,ugs,vgs,xit,ru,rv,
     1              der,g,u1,v1,u2,v2,ugn,vgn,n,m,
     1              keypa,aps,ierr,zerom,distl,0)
      p(1) = ugn
      p(2) = vgn
      if(fret.lt.disto)then
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
      dist =  fret
      ug   = p(1)
      vg   = p(2)
      return
      end
C
C=END SOURCE
