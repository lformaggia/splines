C=NAME evapa2
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/evapa2.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:35 $
C=TYPE subroutine
C=PURPOSE interface
C=KEYWORDS parametric_bicubic_splines Parametric_surfaces
C=BLOCK ABSTRACT
C
C It computes the coefficents of the polinomial expansion for a patch starting
C from the  data rapresentation of a  bicubic
C
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
C
C call EVAPA2(n,m,ndimn,tanret,choret,coor,apatch,i,j)
c
C
C  INPUT       DIMENSION      DESCRIPTION
C  n                        N. of knots along the u direction
C  m                        N. of knots along the v direction
C  ndimn                    number of dimensions
c  tanret     ndimn,3,n,m    tangent vector at the knots of the bicubic
c  choret     2,n,m         chord length for each patch
c  coor       ndimn,*       coordinates of knots
c  i,j                      position of the patch in the lattice on the
c                           parametric plane
c
C OUTPUT      DIMENSION     DESCRIPTION
c
c apatch      ndimn,4,4     coefficents of the patch:
c                           r(id) = SUM apatch(id,i,j)u^(i-1)v^(j-1)
C                                     i,j
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine evapa2(n,m,ndimn,tanret,choret,coor,apatch,
     1                  i,j)
      implicit none
c
c evaluates the matrix a of the bicubic patch from
c the matrix q of the ferguson formulation.
c or reading it from the common area surpa
c
c
c                 a = cqc*
c
c (we have preferred a straightforward matrix moltiplication
c  type of implementation because easier to program, though more
c  expensive)
c
      integer*8 i,i00,i01,i1,i10,i11,id,j,k1,k2,l,m,n,ndimn,iretp
      real*8 a,ab,ad,b,bc,c,cd,d
      real*8 q(4,4),tanret(ndimn,3,n,m),choret(2,n,m)
      real*8 apatch(ndimn,4,4),coor(ndimn,*)
      real*8 cc(4,4)
c
      iretp(k1,k2)=(k2-1)*n+k1
      data cc/1.d0,0.d0,-3.d0, 2.d0,
     2       0.d0,0.d0, 3.d0,-2.d0,
     3       0.d0,1.d0,-2.d0, 1.d0,
     4       0.d0,0.d0,-1.d0, 1.d0/
c
c     if(ispa.eq.0)then
      i00= iretp(i,j)
      i10= iretp(i+1,j)
      i01= iretp(i,j+1)
      i11= iretp(i+1,j+1)
      a  = choret(1,i,j)
      b  = choret(2,i+1,j)
      c  = choret(1,i,j+1)
      d  = choret(2,i,j)
      ab = a*b
      bc = b*c
      cd = c*d
      ad = a*d
c
      do 10 id=1,ndimn
        q(1,1) =    coor(id,i00)
        q(2,1) =    coor(id,i10)
        q(3,1) =  a*tanret(id,1,i  ,j)
        q(4,1) =  a*tanret(id,1,i+1,j)
        q(1,2) =    coor(id,i01)
        q(2,2) =    coor(id,i11)
        q(3,2) =  c*tanret(id,1,i  ,j+1)
        q(4,2) =  c*tanret(id,1,i+1,j+1)
        q(1,3) =  d*tanret(id,2,i  ,j)
        q(2,3) =  b*tanret(id,2,i+1,j)
        q(3,3) = ad*tanret(id,3,i  ,j)
        q(4,3) = ab*tanret(id,3,i+1,j)
        q(1,4) =  d*tanret(id,2,i  ,j+1)
        q(2,4) =  b*tanret(id,2,i+1,j+1)
        q(3,4) = cd*tanret(id,3,i  ,j+1)
        q(4,4) = bc*tanret(id,3,i+1,j+1)
        do 11 i1=1,4
          do 12 l=1,4
             apatch(id,i1,l)=
     1              cc(i1,1)*(q(1,1)*cc(l,1)+
     1            q(1,2)*cc(l,2)+q(1,3)*cc(l,3)+
     1            q(1,4)*cc(l,4))+
     1              cc(i1,2)*(q(2,1)*cc(l,1)+
     1            q(2,2)*cc(l,2)+q(2,3)*cc(l,3)+
     1            q(2,4)*cc(l,4))+
     1              cc(i1,3)*(q(3,1)*cc(l,1)+
     1            q(3,2)*cc(l,2)+q(3,3)*cc(l,3)+
     1            q(3,4)*cc(l,4))+
     1              cc(i1,4)*(q(4,1)*cc(l,1)+
     1            q(4,2)*cc(l,2)+q(4,3)*cc(l,3)+
     1            q(4,4)*cc(l,4))
12        continue
11      continue
10    continue
c     else
c       ist=(j-1)*16*ndimn*n+1
c       do 112 jj=1,4
c       do 112 ii=1,4
c       do 112 id=1,ndimn
c          apatch(id,ii,jj)=aps(ist)
c          ist = ist+1
c12   continue
c     endif
      return
      end
C
C=END SOURCE
