C=NAME xminli
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/xminli.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:41 $
C=TYPE real function
C=PURPOSE Local_to locps4
C=KEYWORDS Parametric_polynimial_Surfaces Line_search_Minimization
C=BLOCK ABSTRACT
c
c This routine finds the minimum of the function
c
c          DISTPS = || r - xp ||
c
c where xp is a point on the 3d space and r a point on a parametric surface
c whose coefficents are contained into aps along a line starting from a given point
c in the parametric plane. That is r is here a function of a single parameter:
c
c     r = R(x) = r(u(x),v(x)) where
c
c     u(x) = ugs + x*xi(1)   v(x)  = vgs + x*xi(2)
c
c ugs,vgs are the coordinate of the startinmg point, while xi is a unitary vector
c containing the search direction. That is the distance function we here seek to mimimize
c is
c
c         F(x) = ||R(x) -xp|| = || r( ugs+x*xi(1),vgs+x*xi(2) ) - xp||
C
C=END ABSTRACT
C=BLOCK USAGE
C
c dist= xminli(utl,xmin,xp,r,z,ugs,vgs,xi,ru,rv,der,g,u10,u20,v10,v20,ug,vg,n,m,
c             keypa,aps,imes,zerom,distl,ind)
C
C
c   INPUT     DIMENSION       DESCRIPTION
c
c n                         n. of knots on the surface along u-coordinate
c m                         n. of knots on the surface along v-coordinate
c aps         3,*           surface coefficent in standard format.
c keypa       (n-1)(m-1)    array of pointers to aps.
c xp          3             coordinates of given point.
c ugs,vgs                   starting value of the global parametric coordinates
c xi          2             search direction
c x                         coordinate along the search direction
c u10,u20                   limit of interval where the minimum is searched (u-coord)
c v10,v20                   limit of interval where the minimum is searched (v-coord)
c utl                       tolerance on the variation of the parametric coordinates
c                           between two iteration of the procedure (see NOTES).
c distl                     tolerance on the minimum value of distance (see NOTES).
c zerom                     machine zero (suggested: 1.e-5 for 32 bit machines)
c ind                         0 -> derivative information not required
c                           /=0 -> derivative information given
c
c
c   OUTPUT     DIMENSION       DESCRIPTION
c
C ug,vg                     parametric coordinates u at found minimum point.
C r           3             found minimum r=r(ug,vg).
c xminli                    distance ||r(ug,vg) - xp||
c z           3             r(ug,vg) - xp
c ru          3             partial derivativo of position vector: r,u(ug,vg)
c rv          3             partial derivative of position vector: r,v*ug,vg)
c imes                      error message (see MESSAGES and DIAGNOSTICS).
c g           2             gradient of distance function ||r - xp|| at point (ug,vg)
c                           (if ind /=0)
c der                       derivative of distance function along search direction
c                           (if ind /=0)
c iter                      N. of iterations used for convergence (if imes=/-1)
c
c OTHERS
c
c itmax       parameter     max n. of iterations
c zeps        parameter     tiny number (It prevents zero division and it is used
c                           in the termination condition 1 - see below)
c
c
c NOTES:
c
c    during the process the minimum is searched along direction xi. It uses a zero order
c    method. The minimum is bracketed by the routine brauv, which is a modification of
c    routine brasf of Numerical Recipies.
c
c   .........................
c   .Termination conditions .
c   .........................
c
c (1)  |x-0.5(a+b)|<utl*|x| +zeps -0.5*(b-a)    Variation of parametric coordinates below
c                                            given tolerance
c (2)  ||r - xp|| < distl                     Distance below given tolerance
c (3)  iter > itmax                          N. of iterations exceeded
c
c
c  During the process, b and a (b>=a) delimit the interval where the minimun is
c  seeked, i.e. a<=xminli<=b
c
c   MESSAGES and DIAGNOSTICS
c
c   imes =0  termination condition (2): normal convergence reached,
c            POINT IS ON the curve: ||r -xp||<distl
c   imes =2  termination condition (1) but condition (2) is not satisfied:
c            stationary point reached, point is NOT on the curve
c   imes =-1 termination condition (3): ERROR CONDITION: N. OF ITERATION EXCEEDED
c                                       ---------------
c
C
C=END USAGE
C
C=BLOCK SOURCE
C
      real*8 function xminli(utl,xmin,xp,r,z,ugs,vgs,xi,ru,rv,
     1                      der,g,u10,v10,u20,v20,ug,vg,n,m,keypa,
     1                      aps, imes,zerom,distl,ind)
      implicit none
c
c this routine finds the minimum of the distance function
c
c          distps  = || r - xp || along the line starting from point
c                                 ugs,vgs and of direction xi
c
c r = r(xmin) = r(u(xmin),v(xmin))  u(xmin) = ugs + xmin*xi(1)
c                                   v(xmin) = vgs + xmin*xi(2)
c
c xp is fixed a point in the 3d space and r a point on the spline surface
c it is a modificationof the routine for 1d search
c in numerical recipes (brent)
c
      integer*8 imes,ind,itmax,iter,keypa(*),m,n
      real*8 a,aps(3,*),ax,b,bl1,bl2,blx,bx,cgold,cx,d,der,distl
      real*8 distps,dx,e,etemp,fa,fb,fc,fu,fv,fw,fx,g(2),p,q,r(3)
      real*8 rr,tol1,tol2,u,u10,u20,ug,ugs,utl,v,v10,v20,vg,vgs,w,x
      real*8 xi(2),xmin,xm,xp(3),z(3),ru(3),rv(3),zerom,zeps
      parameter (itmax=100,zeps=1.0d-10,cgold=.3819660d0)
c
      imes=0
c
c first find limits for beta
c
      call limiuv(ugs,vgs,xi,u10,u20,v10,v20,bl1,bl2)
c
c now bracket the minimum
c
      blx = 0.05*(bl2-bl1)
      ax = max(-blx,bl1)
      bx = min(blx,bl2)
      call brauv(ax,bx,cx,fa,fb,fc,ugs,vgs,keypa,aps,n,m,
     1                 xp,xi,bl1,bl2)
c
c now dist(bx ) <= dist(ax) and dist(bx) <= dist(cx)
c
      ax = max(bl1,min(bl2,ax))
      bx = max(bl1,min(bl2,bx))
      cx = max(bl1,min(bl2,cx))
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.d0
      fx=distps(xp,ugs,vgs,xi,x,keypa,aps,n,m,z,r,der,g,
     1                     ru,rv,ug,vg,0)
c
      if(fx.lt.distl)then
         goto 3
      endif
c
      if(2.d0*(b-a).le.utl*(abs(a)+abs(b)))then
c
c we have reached an extrema
c
         imes =2
         goto 3
      endif
c
      fv=fx
      fw=fx
c
      do 11 iter=1,itmax
        xm=0.5d0*(a+b)
        tol1=utl*abs(x)+zeps
        tol2=2.d0*tol1
c
c convergence up to tol
c
        if(abs(x-xm).le.(tol2-.5d0*(b-a)))then
           imes =2
           goto 3
        endif
        if(abs(e).gt.tol1) then
c
c parabolic fit
c
          rr=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*rr
          q=2.d0*(q-rr)
          if(q.gt.0.d0) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(.5d0*q*etemp).or.p.le.q*(a-x).or.
     *        p.ge.q*(b-x)) goto 1
c
c if parabolic step is ok: take it going to 1
c
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
c
c parabolic step no ok: take golden step
c
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=cgold*e
2       if(abs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
        u=max(bl1,min(bl2,u))
        fu=distps(xp,ugs,vgs,xi,u,keypa,aps,n,m,z,r,dx,g,
     1                     ru,rv,ug,vg,0)

c
c housekeeping
c
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
        if(fx.lt.distl)goto 3
11    continue
c
c number of iteration exceeded
c
      imes =1
c
c normal termination
c
3     xmin=max(bl1,min(x,bl2))
      ug=max(u10,min(ugs+xmin*xi(1),u20))
      vg=max(v10,min(vgs+xmin*xi(2),v20))
      xminli=fx
      return
      end
C
C=END SOURCE
