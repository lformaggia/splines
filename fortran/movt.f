C=NAME movt
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/movt.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:39 $
C=TYPE subroutine
C=PURPOSE local_to getta3
C=KEYWORDS parametric_bicubic_splines
C=BLOCK ABSTRACT
C
C It inserts a value into the imposed tangent vector timp
C
C=END ABSTRACT
C=BLOCK USAGE
C
C call movet(ndimn,timp,l,tn,k,i,j,n,m)
C
C  INPUT       DIMENSION      DESCRIPTION
c
c ndimn                     N. of dimensions
c n,m                       N. of knots along u and v directions respectively
c tn          ndimn,3,n,m   array that will contain the derivative value we want to impose
c l                         index
c i                         index
c j                         index
c k                         derivative considered:
c                           1 -> x,u   2-> x,v   3 -> x,uv
c
C OUTPUT      DIMENSION     DESCRIPTION
c
c timp        ndimn,2       imposed tangent value
c
c
c This routine inserts into timp a value for the tangent. It is used by GETTA3 for
c the imposition of the boundary value of the tangent to the iso-curve which will
c form the basis of the tensor-product surface.
c
c      timp(*,l) = tn (*,k,i,j)
c
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine movt(ndimn,timp,l,tn,k,i,j,n,m)
      implicit none
      integer*8 id,i,j,k,l,m,n,ndimn
      real*8 timp(ndimn,2),tn(ndimn,3,n,m)
      do 10 id=1,ndimn
         timp(id,l)=tn(id,k,i,j)
10    continue
      return
      end
C
C=END SOURCE
