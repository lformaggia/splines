C=NAME patc1d
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/patc1d.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:39 $
C=TYPE subroutine
C=PURPOSE local
C=KEYWORDS parametric_cubic_splines
C=BLOCK ABSTRACT
C
c It gives the array of coefficients of a cubic spline arc
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
c call patc1d(ndimn,apatch,q,t,cs,is)
c
C  INPUT       DIMENSION      DESCRIPTION
c
c ndimn                     N. of dimensions
c q           ndimn,n       Knots coordinates
c cs          n             Arc chord length
C t           ndimn,n       Tangent at knots
c is                        Number of arc ( 1<= is <=n-1)
C OUTPUT      DIMENSION     DESCRIPTION
c
c apatch      ndimn,4       array of coefficents A so that
c                           r(u) = UA
c
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine patc1d(ndimn,apatch,q,t,cs,is)
c
c  gives the array of coefficents corresponding to
c  a cubic spline segments
c
      real      q(ndimn,*),t(ndimn,*),cs(*)
      real      apatch(ndimn,4)
c
      clen = cs(is)
      ir = is+1
      do 10 id=1,ndimn
         r12 = q(id,ir-1) - q(id,ir )
         t1  = t(id,ir-1)
         t2  = t(id,ir  )
         apatch(id,4)  =  2*r12 + clen*(t1+t2)
         apatch(id,3)  = -3*r12 - clen*(2*t1+t2)
         apatch(id,2)  =  clen*t1
         apatch(id,1)  =  q(id,ir-1)
 10   continue
      return
      end
C
C=END SOURCE
