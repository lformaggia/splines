C=NAME buiaps
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/buiaps.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:34 $
C=TYPE subroutine
C=AUTHOR L.formaggia
C=PURPOSE General
C=KEYWORDS Parametric_bicubic_spline Parametric_polynomial_surfaces Data_structure_interface
C
C=BLOCK ABSTRACT
C
c It converts the rapresentation of a bicubic spline surface into the general
c rapresentation for a polynomial surface
C
C=END ABSTRACT
C=BLOCK USAGE
C
C    call buiaps(n,m,tanret,choret,coor,keypa,lastpa,aps,mpatch,msul,ierr)
C
C  INPUT      DIMENSION     DESCRIPTION
C  n                        n. of knots along u-coordinate
c  m                        n. of knots along v-coordinate
c  tanret     3,3,n,m       derivative vector at knots (see routine getta3)
C  choret     2,n,m         chord length for each patch (see routine getta3)
C  coor       3,n*m          knots coordinates
C  mpatch                   max. number of patches allowed on surface
C  msul                     max. dimension of array of scoefficents aps
c
C
C OUTPUT      DIMENSION     DESCRIPTION
c
C  keypa      mpatch        array of pointers to aps
C  aps        3,msul        array of coefficents
C  lastpa                   last location in aps used
C  ierr                                   error condition
C
C   * DIAGNOSTICS *
c
c ierr = 0    no error
c ierr = 1    number of patches exceeded (mpatch<(n-1)*(m-1))
c ierr = 2     number of coefficents exceeded (msul<17*(n-1)*(m-1))
c ierr = 3    1+2
C
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine buiaps(n,m,tanret,choret,coor,keypa,lastpa,aps,mpatch,
     1                  msul,ierr)
      implicit none
      integer*8 i,id,ierr,ii,j,jj,keypa(*),lastpa,m,mpatch,msul,n,npa
      real*8 tanret(3,3,n,m),choret(2,n,m),coor(3,*)
      real*8 apatch(3,4,4),aps(3,*),xl
c
c check
c
      ierr =0
      if(mpatch.lt.(n-1)*(m-1))ierr=1
      if(msul.lt.17*(n-1)*(m-1))ierr=ierr+2
      if(ierr.ne.0)return
      lastpa = 1
      npa = 0
c
c start loop on patches
c
      do i=1,n-1
      do j=1,m-1
        npa = npa +1
        keypa(npa) = lastpa
c
        xl=max(sqrt(choret(1,i,j)**2+choret(2,i,j  )**2),
     1         sqrt(choret(1,i,j)**2+choret(2,i+1,j)**2))
        aps(1,lastpa)=4
        aps(2,lastpa)=4
        aps(3,lastpa)=xl
c
c evaluate coefficents
c
        call evapa2(n,m,3,tanret,choret,coor,apatch,i,j)
c
        lastpa = lastpa +1
c
c store coefficents
c
        do jj=1,4
        do ii=1,4
          lastpa = lastpa +1
          do id =1,3
            aps(id,lastpa)=apatch(id,ii,jj)
          end do
        end do
        end do
c
      end do
      end do
      return
      end
C
C=END SOURCE
