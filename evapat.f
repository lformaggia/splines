C=NAME evapat
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/evapat.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:36 $
C=TYPE subroutine
C=PURPOSE interface
C=KEYWORDS parametric_bicubic_splines Parametric_surfaces
C=BLOCK ABSTRACT
C
C It computes the coefficents of the polinomial expansion for a patch starting
C from the  Ferguson rapresentation of a  bicubic
C
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
C
C call EVAPAT(ndimn,qferg,apatch)
c
C
C  INPUT       DIMENSION      DESCRIPTION
C  ndimn                    number of dimensions
c  qferg      ndimn,4,4     Ferguson rapresentation of the bicubic:
c
c
C OUTPUT      DIMENSION     DESCRIPTION
c
c apatch      ndimn,4,4     coefficents of the patch:
c                           r(id) = SUM apatch(id,i,j)u^(i-1)v^(j-1)
C                                     i,j
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine evapat(ndimn,qferg,apatch)
      implicit none
c
c evaluates the matrix a of the bicubic patch from
c the matrix q of the ferguson formulation.
c
c                 a = cqc*
c
c (we have preferred a straightforward matrix moltiplication
c  type of implementation because easier to program, though more
c  expensive)
c
      integer*8 i,id,l,ndimn
      real*8 q(4,4),qferg(ndimn,4,4),apatch(ndimn,4,4)
      real*8 c(4,4)
      data c/1.d0,0.d0,-3.d0, 2.d0,
     2       0.d0,0.d0, 3.d0,-2.d0,
     3       0.d0,1.d0,-2.d0, 1.d0,
     4       0.d0,0.d0,-1.d0, 1.d0/
c
      do 10 id=1,ndimn
        do 9 i=1,4
          do 8 l=1,4
            q(i,l)=qferg(id,i,l)
 8        continue
 9      continue
        do 11 i=1,4
          do 12 l=1,4
             apatch(id,i,l)=
     1              c(i,1)*(q(1,1)*c(l,1)+
     1            q(1,2)*c(l,2)+q(1,3)*c(l,3)+
     1            q(1,4)*c(l,4))+
     1              c(i,2)*(q(2,1)*c(l,1)+
     1            q(2,2)*c(l,2)+q(2,3)*c(l,3)+
     1            q(2,4)*c(l,4))+
     1              c(i,3)*(q(3,1)*c(l,1)+
     1            q(3,2)*c(l,2)+q(3,3)*c(l,3)+
     1            q(3,4)*c(l,4))+
     1              c(i,4)*(q(4,1)*c(l,1)+
     1            q(4,2)*c(l,2)+q(4,3)*c(l,3)+
     1            q(4,4)*c(l,4))
12        continue
11      continue
10    continue
      return
      end
C
C=END SOURCE
