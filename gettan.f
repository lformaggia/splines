C=NAME gettan
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/gettan.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:37 $
C=TYPE subroutine
C=PURPOSE general
C=KEYWORDS parametric_cubic_splines
C=BLOCK ABSTRACT
C
C It computes the  1st derivative at a point in a cubic spline curve.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C call gettan(ndimn,n,q,t,cs,is,u,xu)
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
c xu          ndimn         x,u = 1st derivative with respect to parameter u at the
c                           point
c
c
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine gettan(ndimn,n,q,t,cs,is,u,tp)
c
c  tp -> tangent at point
c
      real      q(ndimn,*),t(ndimn,*),cs(*),tp(*)
c
      clen = cs(is)
      ir = is+1
      do 10 id=1,ndimn
         r12 = q(id,ir-1) - q(id,ir )
         t1  = t(id,ir-1)
         t2  = t(id,ir  )
         a1  =  6*r12 + 3*clen*(t1+t2  )
         a2  = -6*r12 - 2*clen*(2*t1+t2)
         tp(id)= (   a1                      *u +
     1               a2                     )*u +
     1               clen*t1
 10   continue
      return
      end
C
C=END SOURCE
