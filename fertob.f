C=NAME fertob
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/fertob.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:36 $
C=TYPE subroutine
C=PURPOSE general
C=KEYWORDS parametric_bicubic_splines
C=BLOCK ABSTRACT
C
C It computes the Bezier control points for a patch from the value of the
C coefficents of the Ferguson tensor-product rapresentation
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
C call fertob(ndimn,qferg,bbez)
C
C
C  INPUT       DIMENSION      DESCRIPTION
c
c ndimn                     n. of dimensions
c qferg       ndimn,4,4     coefficents of the ferguson rapresentation:
c
c                                   ! r  (0,0) r  (0,1) r,v (0,0) r,v (0,1) !
c                           QFERG=  ! r  (1,0) r  (1,1) r,v (1,0) r,v (1,1) !
c                                   ! r,u(0,0) r,u(0,1) r,uv(0,0) r,uv(0,1) !
c                                   ! r,u(1,0) r,u(1,1) r,uv(1,0) r,uv(1,1) !
c
C OUTPUT      DIMENSION     DESCRIPTION
c
c bbez        ndimn,4,4     array of control points for the Beziet rapresentation
c                           of the patch
c
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine fertob(ndimn,qferg,bbez)
      parameter(c13=1./3.,c19=1./9.)
      real qferg(ndimn,4,4),bbez(ndimn,4,4)
      integer i
      q(k,j) = qferg(i,k,j)
c
      do 10 i=1,ndimn
        bbez(i,1,1) = q(1,1)
        bbez(i,1,2) = q(1,1)+c13*q(1,3)
        bbez(i,1,3) = q(1,2)-c13*q(1,4)
        bbez(i,1,4) = q(1,2)
c
        bbez(i,2,1) = q(1,1)+c13*q(2,1)
        bbez(i,2,2) = q(1,1)+c13*(q(1,3)+q(3,1))+c19*q(3,3)
        bbez(i,2,3) = q(1,2)+c13*(q(3,2)-q(1,4))-c19*q(3,4)
        bbez(i,2,4) = q(1,2)+c13*q(3,2)
c
        bbez(i,3,1) = q(2,1)-c13*q(4,1)
        bbez(i,3,2) = q(2,1)+c13*(q(2,3)-q(4,1))-c19*q(4,3)
        bbez(i,3,3) = q(2,2)-c13*(q(2,4)+q(4,2))+c19*q(4,4)
        bbez(i,3,4) = q(2,2)-c13*q(4,2)
c
        bbez(i,4,1) = q(2,1)
        bbez(i,4,2) = q(2,1)+c13*q(2,3)
        bbez(i,4,3) = q(2,2)-c13*q(2,4)
        bbez(i,4,4) = q(2,2)
10    continue
      return
      end
C
C=END SOURCE
