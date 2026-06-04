C=NAME getpoi
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/getpoi.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:37 $
C=TYPE subroutine
C=PURPOSE general
C=KEYWORDS parametric_cubic_splines
C=BLOCK ABSTRACT
C
C It computes the position vector, 1st derivative
c  and an estimate of the curvature at a point in a cubic spline curve.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C call getpoi(ndimn,n,q,t,cs,is,u,x,xu,xk)
C
C  INPUT       DIMENSION      DESCRIPTION
c
c ndimn                     N. of dimensions
c n                         n. of knots
c q           ndimn,n       Knots coordinates
c cs          n             Arc chord length
C t           ndimn,n       Tangent at knots
c is                        Number of arc ( 1<= is <=n-1)
c u                         Parametric coordinate, local to arc is (0.<=u<=1.)
c
C OUTPUT      DIMENSION     DESCRIPTION
c
c xu          ndimn         x,u 1st derivative with respect to parameter u at the
c                           point
c xk                        estimate of curvature at the point
c
c NOTE :
c This routine is somehow obsolete. It should be subsituted by a call to getp2 and getk2.
c the value of the curvature is only an estimate. The major difference with getk2 is that
c it can be computed for an arbitrary number of dimension ndimn (>=2) and that it
c is faster.
c The curvature is estimated by xk = ||x,uu||/(c**2), where x,uu is the 2nd derivative and
c c the chord length. It can be a good estimate for planar curves.
c
c
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine getpoi(ndimn,n,q,t,cs,is,u,x,tp,xk)
c
c  gives the pint corresponding to parameter u of
c  segment is of the cubic spline
c
c   x-> point      tp -> tangent at point
c
      real      q(ndimn,*),t(ndimn,*),cs(*),x(*),tp(*)
c
      xk = 0.
      clen = cs(is)
      ir = is+1
      do 10 id=1,ndimn
         r12 = q(id,ir-1) - q(id,ir )
         t1  = t(id,ir-1)
         t2  = t(id,ir  )
         a1  =  2*r12 + clen*(t1+t2)
         a2  = -3*r12 - clen*(2*t1+t2)
         a3  =  clen*t1
         x(id)=
     1         ( ( a1                        *u +
     1             a2                       )*u +
     1             a3                       )*u +
     1              q(id,ir-1)
         tp(id)= ( 3*a1                      *u +
     1             2*a2                     )*u +
     1               a3
         xk    =   xk + (6*a1*u + 2*a2)**2
 10   continue
      xk = sqrt(xk)/(clen*clen)
      return
      end
C
C=END SOURCE
