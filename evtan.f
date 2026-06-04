C=NAME evtan
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/evtan.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:36 $
C=TYPE subroutine
C=PURPOSE special to psplin
C=KEYWORDS parametric_cubic_splines
C=BLOCK ABSTRACT
C
C It computes tangents at knots points for a parametric cubic spline by
C solving a tri-diagonal system of equations
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
C call evtan(ndimn,n,q,cs,a,b,c,t,ispty,timp)
C
C
C  INPUT       DIMENSION      DESCRIPTION
c
c ndimn                     n. of dimensions
c n                         n. of knots
c q           ndimn,n       knots coordinates
c cs          n             arc chord length
c ispty       2             end condition marker for spline:
c                           =0 -> natural spline !=0 ->tangent imposed
c timp        ndimn,2       value of tangent to be imposed to the corresponding
c                           end of the spline (if the relative value of ispty !=0)
c
c HELP ARRAYS
c
c a           n
c b           n             (a,b and c are use by the tridiagonal solver)
c c           n
c
c
C OUTPUT      DIMENSION     DESCRIPTION
c
C t           ndimn,n       tangent at knots
c
c
c ROUTINES CALLED
C
C  trcoe  - set coefficents for the tridiagonal system
C  trif   - L-U decomposition of tri-diagonal systems (Thomas Algorithm)
C  tris2  - Backsubstitution for tri-diagonal systems
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine evtan(ndimn,n,q,cs,a,b,c,t,ispty,timp)
      implicit none
      integer*8 i,id,ispty(*),n,ndimn
      real*8 a(*),b(*),c(*),cmax,cs(*),q(ndimn,*),t(ndimn,*)
      real*8 timp(ndimn,2),zero
      parameter(zero=1.d-10)
c
      cmax = cs(1)
      do 1 i=2,n-1
       cmax = max(cmax,cs(i))
 1    continue
      if(abs(cmax).lt.zero)then
        do 4 i=1,n
        do 4 id=1,ndimn
                         t(id,i)=0.d0
4       continue
        return
      endif
c
c calculate tridiagonal system coeficent
c
      call trcoe(ndimn,n,q,cs,a,b,c,t,ispty,timp)
c
c evaluate tangents
c
      call trif(n,a,b,c)
      call tris2(n,ndimn,a,b,c,t)
      return
      end
C
C=END SOURCE
