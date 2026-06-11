C=NAME gett2
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/gett2.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:37 $
C=TYPE subroutine
C=PURPOSE general
C=KEYWORDS parametric_cubic_splines
C=BLOCK ABSTRACT
C
C It computes the tangent vector at a knot of a spline: T = t/s, where s = ||r,u||
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
C call gett2(ndimn,t,xd,s)
C
C
C  INPUT       DIMENSION      DESCRIPTION
c
c ndimn                     N. of dimensions (IT MUST BE equal to 2 or 3 !!)
C xd          ndimn,n       Derivative at knots
c s                         ||xd||
c
c
C OUTPUT      DIMENSION     DESCRIPTION
c
c t                         tangent at knots
c
c
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine gett2(ndimn,t,xd,s)
c
c get tangent t= dx/ds
c
      implicit none
      integer*8 i,ndimn
      real*8 s,s1,t(ndimn),xd(ndimn)
      s1 = 1.d0/s
      do 10 i=1,ndimn
10    t(i) = xd(i)*s1
      return
      end
C
C=END SOURCE
