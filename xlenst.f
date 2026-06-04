C=NAME xlenst
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/xlenst.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:40 $
C=TYPE real function
C=PURPOSE general
C=KEYWORDS parametric_Polynomial_curvea
C=BLOCK ABSTRACT
C
c It returns the length of a segment of parametric curve between the points
c r(u1) and r(u2)
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
c real function xlenst(xln1a,u1,u2,eps,prec,in,ierr)
c
C  INPUT       DIMENSION      DESCRIPTION
c
c xln1a       *             array of coefficient for the curve, in standard
c                           format, FOR THE GIVEN ARC
c u1                        local parametric coordinate of first point
c u2                        local parametric coordinate of second point
c eps                       tolerance (see convergence criterion)
c in                        convergence type
c                           (see convergence criterion)
c ierr                      error indicator:
c                           0 -> no error
c                           n -> max n. iteration exceeded (n=nit+1)
c
C OUTPUT      DIMENSION     DESCRIPTION
c
c xlenst                    lenght of curve between u1 and u2
c prec                      difference of the value of the approximation in the last
c                           2 consecutive iterations
c
c PARAMETERS
c
c nit                       max number of iterations
c
c
c ** WARNING **
c
c the routine computes the integral only within a given arc. This because only
c within an arc the expression for the arc length is regular enough to assure a good
c accuracy of the result. If we need to compute the lenght of a curve within two given
c GLOCBAL COORDINATES ws1 and ws2 we may use the following loop:
c
c
c          ns1 = int(w1)
c          ns2 = int(w2)
c          tol= relative_tolerance
c          tlen=0.
c          do 500 is=ns1,ns2
c            ista = isn(is)
c            u1 = max(0.,w1-is)
c            u2 = min(1.,w2-is)
c            tlen = tlen + xlenst(xln1d(1,ista),u1,u2,tol,prec,0,ierr)
c500       continue
c
c isn
c
c  ** NOTE **
c
c It uses Romberg adaptive integration.
c
c BEWARE: if u2 < u1 then xlenst<0.
c
c
c ** CONVERGENCE CRITERION **
c
c |prec| = |s-os| <= tolerance
c
c where os and s are the values of the lenght computed in the last two iterations
c
c if in  =0 then tolerance =eps*|os|    (relative tolerance)
c if in /=0 then tolerance =eps         (absolute tolerance)
c
c ERROR CONDITIONS
c
c If number of iterations is exceeded ierr returns nit+1, where nit is a parameter
c setting the maximum number of iterations. Normally a value of nit=10 suffices.
c If number of iteration is exceeded then we have a very irregular curve. In this case it
c is better to split the integral into two integrals over (u1,um) amd (um,u2),
c (um=(u1+u2)/2).
c
C=END USAGE
C
C=BLOCK SOURCE
C
      real*8 function xlenst(xln1d,u1,u2,eps,prec,in,ierr)
      implicit none
c
c *** this sub. computes the length of a segment of a spline.
c     using a romberg adaptive integration
c
      integer*8 ierr,in,it,jt,kt,nit
      parameter(nit=20)
      real*8 aa,eps,eps1,f1,f2,os,ost,prec,s,st,sum,tnm,u1,u2,u21
      real*8 x,xd(3),xddd(3),xdd(3),xln1d(*),xp(3),del
c
      eps1 = eps
      ierr = 0
      os   = -1.d30
c
      call evps1d(xln1d,u1,xp,xd,xdd,xddd,f1,-1)
      call evps1d(xln1d,u2,xp,xd,xdd,xddd,f2,-1)
c
      u21 = u2-u1
      st = 0.5d0*u21*(f1+f2)
      ost = st
      kt = 1
      do 200 it=1,nit
      tnm = dble(kt)
      del = u21/tnm
      x = u1+0.5d0*del
      sum = 0.d0
      do 100 jt=1,kt
      call evps1d(xln1d,x,xp,xd,xdd,xddd,aa,-1)
      sum = sum+aa
      x = x+del
  100 continue
      st = 0.5d0*(st+del*sum)
      kt = kt*2
      s = (4.d0*st-ost)/3.d0
      if(in.eq.0) eps1 = eps*abs(os)
      prec = abs(s-os)
      if(prec.le.eps1) goto 300
      os = s
      ost = st
  200 continue
      ierr = nit+1
  300 xlenst=s
      return
      end
C
C=END SOURCE
