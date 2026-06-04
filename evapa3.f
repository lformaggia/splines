C=NAME evapa3
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/evapa3.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:35 $
C
C=TYPE subroutine
C=PURPOSE general
C
C=KEYWORDS Parametric_bicubic_splines
C
C=BLOCK ABSTRACT
C It provides the matrix A containing the patch coefficent for
C a given patch. It is a modification of EVAPA2. The difference here is that the
C routine fills
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
C   call evapa3(n,m,ndimn,tanret,choret,coor,apatch,i,j)
c
c INPUT       DIMENSION     DESCRIPTION
c  n
c  m
c  ndimn
c  tanret     ndimn,3,n,m
c  choret     2,n,m
c  coor       ndimn,*
c  i,j                      patch identifiers: they correspond to the identifiers
c                           of the bottom-left knot in the patch under consideration
c
c   OUTPUT    DIMENSION     DESCRIPTION
c
c  apatch     ndimn,0:16    coefficent matrix (If ndimn=3 we have in the header:
c                           apatch(1,0)=n; apatch(2,0)=m; apatch(3,0)=xl
c                           where xl is a characteristic dimension of the patch)
c
c
c NOTES
c a point of parametric coordinate (u,v) will be computed as
c
c    r(u,v) = [1,u,u**2,u**3] * apatch *  | 1  |
c                                         | v  |
c                                         |v**2|
c                                         |v**3|
C
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine evapa3(n,m,ndimn,tanret,choret,coor,apatch,
     1                  i,j)
c
c evaluates the matrix a of the bicubic patch and stores it into the
c standard format
c
c
      real tanret(ndimn,3,n,m),choret(2,n,m)
      real apatch(ndimn,0:*),coor(ndimn,*)
c
      iretp(k1,k2)=(k2-1)*n+k1
c
      call evapa2(n,m,ndimn,tanret,choret,coor,apatch(1,1),
     1            i,j)
c
      if(ndimn.gt.2)then
c
c add header to apatch
c
        apatch(1,0) = 4
        apatch(2,0) = 4
        if(ndimn.eq.3)then
          i00= iretp(i,j)
          i10= iretp(i+1,j)
          i01= iretp(i,j+1)
          i11= iretp(i+1,j+1)
          xl1 = sqrt((coor(1,i11)-coor(1,i00))**2+
     1               (coor(1,i11)-coor(1,i00))**2+
     1               (coor(1,i11)-coor(1,i00))**2)
          xl2 = sqrt((coor(1,i01)-coor(1,i10))**2+
     1               (coor(1,i01)-coor(1,i10))**2+
     1               (coor(1,i01)-coor(1,i10))**2)
          apatch(3,0) = max(xl1,xl2)
        endif
      endif
      return
      end
C
C=END SOURCE
