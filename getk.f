C=NAME getk
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/getk.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:37 $
C=TYPE subroutine
C=PURPOSE general
C=KEYWORDS parametric_cubic_splines
C=BLOCK ABSTRACT
C
C It computes the curvature a point on a cubic spline
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
C call getk(ndimn,n,q,t,cs,is,u,xk)
C
C
C  INPUT       DIMENSION      DESCRIPTION
c
c ndimn                     N. of dimensions (IT MUST BE equal to 2 or 3 !!)
c n                         n. of knots
c q           ndimn,n       Knots coordinates
c cs          n             Arc chord length
C t           ndimn,n       Tangent at knots
c is                        Number of arc ( 1<= is <=n-1)
c u                         Parametric coordinate, local to arc is (0.<=u<=1.)
c
C OUTPUT      DIMENSION     DESCRIPTION
c
c xk                        curve curvature
c
c NOTE:
c
c  The values of the curvature is computed by using the Frenet-Serret formula.
c
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine getk(ndimn,n,q,t,cs,is,u,xk)
c
c  gives the curvature at the point corresponding to u of
c  segment is of a  cubic spline
c
c
      real      q(ndimn,*),t(ndimn,*),cs(*)
c
      clen = cs(is)
      xk = 0.
      ir = is+1
      do 10 id=1,ndimn
         r12 = q(id,ir-1) - q(id,ir )
         t1  = t(id,ir-1)
         t2  = t(id,ir  )
         a1  = 12*r12 + 6*clen*(t1+t2  )
         a2  = -6*r12 - 2*clen*(2*t1+t2)
         xk  = xk  + (a1*u +  a2)**2
 10   continue
      xk =sqrt(xk)/(clen*clen)
      return
      end
C
C=END SOURCE
