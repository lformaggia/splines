C=NAME evaqfg
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/evaqfg.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:36 $
C=TYPE subroutine
C=PURPOSE interface
C=KEYWORDS parametric_bicubic_splines
C=BLOCK ABSTRACT
C
C From the internal representation of a bicubic surface patch it
c builds the matrix Q of the Ferguson patch
C
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
C call evaqfg(ndimn,tanret,choret,n,m,i,j,qferg,coor)
C
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
c qferg       ndimn,4,4     Matrix Q of the Ferguson representation
C
C
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine evaqfg(ndimn,tanret,choret,n,m,i,j,qferg,
     1                  coor)
c
c it evaluates the matrix q for the ferguson patch (i,j) and
c stores it into qferg
c
      real qferg(ndimn,4,4),tanret(ndimn,3,n,m),choret(2,n,m)
      real coor(ndimn,*)
      iretp(k1,k2)=(k2-1)*n +k1
c
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
      do 10 id=1,ndimn
c
c get the derivatives and the points
c
        qferg(id,1,1) =    coor(id,i00)
        qferg(id,2,1) =    coor(id,i10)
        qferg(id,3,1) =  a*tanret(id,1,i  ,j)
        qferg(id,4,1) =  a*tanret(id,1,i+1,j)
        qferg(id,1,2) =    coor(id,i01)
        qferg(id,2,2) =    coor(id,i11)
        qferg(id,3,2) =  c*tanret(id,1,i  ,j+1)
        qferg(id,4,2) =  c*tanret(id,1,i+1,j+1)
        qferg(id,1,3) =  d*tanret(id,2,i  ,j)
        qferg(id,2,3) =  b*tanret(id,2,i+1,j)
        qferg(id,3,3) = ad*tanret(id,3,i  ,j)
        qferg(id,4,3) = ab*tanret(id,3,i+1,j)
        qferg(id,1,4) =  d*tanret(id,2,i  ,j+1)
        qferg(id,2,4) =  b*tanret(id,2,i+1,j+1)
        qferg(id,3,4) = cd*tanret(id,3,i  ,j+1)
        qferg(id,4,4) = bc*tanret(id,3,i+1,j+1)
10    continue
      return
      end
C
C=END SOURCE
