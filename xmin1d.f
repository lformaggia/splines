C=NAME xmin1d
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/xmin1d.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:41 $
C=TYPE real function
C=PURPOSE general gensur
C=KEYWORDS Parametric_Polynomial_Curves 1D_Minimization
C=BLOCK ABSTRACT
c
c This routine finds the minimum of the function
c
c          FUNDIS  = || r - xp ||
c
c where xp is a point on the 3d space and r a point on a parametric curve
c whose coefficents are contained into xln1d. It makes use of the routine
c for 1D search DBRENT published on Numerical Recipies. The routine requires to give a
c triplet of points ax,bx,cx, which form the so called bracketing absissae.
c THE MINIMUM IS THEN SEARCHED IN THE INTERVAL [ax,cx].
c -----------------------------------------------------
C
C=END ABSTRACT
C=BLOCK USAGE
C
c     real  function xmin1d(ax,bx,cx,tol,xmin,n,xln1d,isn,xp,r,
c    1                      du,imes,utl,distl,ftl,zerom)
C
C
c   INPUT     DIMENSION       DESCRIPTION
c
c ax,bx,cx                  minimum bracketing absissae, such that ax<=bx<cx and
c                           FUNDIS(ax)>FUNDIS(bx) and FUNDIS(cx)>FUNDIS(bx).
c n                         n. of knots on the curve (n. of arcs +1).
c xln1d       3,*           curve coefficent in standard format.
c isn         n-1           array of pointers to xln1d.
c xp          3             coordinates of given point.
c utl                       tolerance on the variation of the parametric coordinate
c                           between two iteration of the procedure(see NOTES).
c distl                     tolerance on the minimum value of FUNDIS(see NOTES).
c ftl                       tolerance on the variation of the value of FUNDIS
c                           between two iteration of the procedure(see NOTES).
c zerom                     machine accuracy (normally put to 1.e-5 for 32 bit)
c
c
c   OUTPUT     DIMENSION       DESCRIPTION
c
c xmin1d                    min(|| r(u) - xp||) found.
C xmin                      parametric coordinate u at minimum point.
C r           3             found minimum r=r(xmin).
c du                        derivative of FUNDIS=|| r(u) - xp|| at minimum point.
c imes                      error message (see MESSAGES and DIAGNOSTICS).
c
c OTHERS
c
c
c zeps        parameter     tiny value ( to avoid division by zero)
c itmax       parameter     max n. of iterations
c
c
c NOTES:
c
c    during the process the variables a,b and x form a bracketing
c    triplet, such that a<x<b and FUNDIS(x) <= FUNDIS(a) and
c    FUNDIS(x) <= FUNDIS(b)
c
c   .........................
c   .Termination conditions .
c   .........................
c
c (1)  | x - 0.5*(a+b)  | <= 2*utl - 0.5*(b-a)  Variation of parametric coordinate below
c                                             given tolerance
c (2)  |fv - fw | < 0.5*ftl*(|fv|+|fw|)         Variation of the distance function FUNDIS
c                                               between two successive iterations
c (3)  ||r - xp|| < distl                     Distance below given tolerance
c (4)  iter > itmax                          N. of iterations exceeded
c
c
c    where a and b and x are the latest values of the bracketing abscissae so that
c    the minimum is in the interval spanned by those 3 values.
c
c   * beware that while utl and ftl are relative values, distl rapresent a absolute
c     tolerance
c
c
c   MESSAGES and DIAGNOSTICS
c
c   imes =0  termination condition (3): normal convergence reached,
c            point on the curve: ||r -xp||<distl
c   imes =1  termination condition (2), but ||r-xp||>=distl: normal convergence
c            but point is NOT on the curve
c   imes =2  termination condition (1) but condition (3) is not satisfied
c            stationary point reached, point is NOT on the curve
c   imes =-1 termination condition (4): ERROR CONDITION: N. OF ITERATION EXCEEDED
c                                       ---------------
c
c
c   NOTES ON THE METHOD
c
c   It uses a modification of a canned routine (DBRENT) which adopts a strategy which
c   makes a limited use of the derivative information. It is slow but relatively robust.
c   It keeps adjornig two variables a,b  which bracket the minimum i.e.
c   a local minimum for ||r - xp|| is contained into [a,b].
c   The routine requires to give in input  a triplet of points ax,bx,cx, so that
c   FUNDIS(ax)>FUNDIS(bx) and FUNDIS(cx)>FUNDIS(bx). This point may be found by using the
c   routine brau.
c
c   THE MINIMUM IS THEN SEARCHED IN THE INTERVAL SPANNED BY ax,bx,cx
C
C
C=END USAGE
C
C=BLOCK SOURCE
C
      real*8 function xmin1d(ax,bx,cx,xmin,n,xln1d,isn,xp,r,
     1                      du,imes,utl,distl,ftl,zerom)
      implicit none
c
c this routine finds the minimum of the function
c
c          fundis  = || r - xp ||
c
c where xp is a point on the 3d space and r a point on the spline
c whose patches are contained into xln1d it makes use of the routine
c for 1d search on numerical recipes (dbrent)
c
      integer*8 imes,isn(*),iter,itmax,n
      real*8 a,ax,b,bx,cx,d,d1,d2,distl,du,dv,dw,dx,e,ftl,fu,fundis
      real*8 fv,fw,fx,olde,tol1,tol2,u,u1,u2,utl,uu0,v,w,x,xln1d(3,*)
      real*8 xmin,xm,xp(3),r(3),zerom,zeps
      parameter (itmax=100,zeps=1.0d-10)
      logical ok1,ok2
      a=min(ax,cx)
      b=max(ax,cx)
      imes =0
      v=bx
      w=v
      x=v
      e=0.d0
      fx=fundis(x,xln1d,isn,n,dx,xp,r)
      fv=fx
      fw=fx
      dv=dx
      dw=dx
c
      uu0=-1.d30
      do 11 iter=1,itmax
        xm=0.5d0*(a+b)
        tol1=utl*abs(x)+zeps
        tol2=2.d0*tol1
        if(fx.lt.distl)then
           goto 3
        endif
        if(abs(x-xm).le.(tol2-.5d0*(b-a)))then
            imes=2
            goto 3
        endif
        if(abs(e).gt.tol1) then
          d1=2.d0*(b-a)
          d2=d1
          if(dw.ne.dx) d1=(w-x)*dx/(dx-dw)
          if(dv.ne.dx) d2=(v-x)*dx/(dx-dv)
          u1=x+d1
          u2=x+d2
          ok1=((a-u1)*(u1-b).gt.0.d0).and.(dx*d1.le.0.d0)
          ok2=((a-u2)*(u2-b).gt.0.d0).and.(dx*d2.le.0.d0)
          olde=e
          e=d
          if(.not.(ok1.or.ok2))then
            go to 1
          else if (ok1.and.ok2)then
            if(abs(d1).lt.abs(d2))then
              d=d1
            else
              d=d2
            endif
          else if (ok1)then
            d=d1
          else
            d=d2
          endif
          if(abs(d).gt.abs(0.5d0*olde))go to 1
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(dx.ge.0.d0) then
          e=a-x
        else
          e=b-x
        endif
        d=0.5d0*e
2       if(abs(d).ge.tol1) then
          u=x+d
          u=max(0.d0,min(u,dble(n)))
          fu=fundis(u,xln1d,isn,n,du,xp,r)
c          if(2*abs(u-uu0).lt.zerom*(abs(u)+abs(uu0)))then
c             imes=-2
c             go to 99
c          endif
        else
          u=x+sign(tol1,d)
          u=max(0.d0,min(u,dble(n)))
          fu=fundis(u,xln1d,isn,n,du,xp,r)
          if(fu.gt.fx)go to 3
c          if(2*abs(u-uu0).lt.zerom*(abs(u)+abs(uu0)))then
c             imes=-2
c             go to 99
c          endif
        endif
        uu0 =u
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          dv=dw
          w=x
          fw=fx
          dw=dx
          x=u
          fx=fu
          dx=du
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            dv=dw
            w=u
            fw=fu
            dw=du
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
            dv=du
          endif
        endif
11    continue
99    continue
c
      if(imes.eq.0)then
         if(fx.lt.distl)goto3
         if(2.d0*abs(fv-fw).lt.ftl*(abs(fv)+abs(fw)))then
           imes =1
         else
           imes=-1
         endif
      endif
c
      xmin1d=fx
      xmin =x
      return
3     xmin=x
      xmin1d=fx
      return
      end
C
C=END SOURCE
