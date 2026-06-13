C=NAME loas2d
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/loas2d.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:38 $
C=TYPE subroutine
C=PURPOSE surface spline
C=KEYWORDS parametric_polynomial_surfaces Data_base_Interface
c
C It reads the geometrical information for a parametric  surface
c from the data base  and stores it into keypa and aps
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
C       call loas2d(keypa,aps,msul,mpat,isur,igsur,ib,n,m,
c                   kfio,ifio,ips,ito,ierr)
c
c  INPUT      DIMENSION     DESCRIPTION
c
c isur                      local surface numbering
c kfio                      i/o unit number associated to the file of keys
c ifio                      i/o unit number associated to the file of coefficents
c mpat                      max. number of patches allowed
c msel                      max dimension of aps
c ito                       i/o unit number associated to error i/o channel
c
c
c  OUTPUT     DIMENSION     DESCRIPTION
c n                         n. of knots along u
c m                         n. of knots along v
c ib                        surface indicator (normally used to indicte the bounday condition
c                           associated to the surface)
c igsur                     global numbering to be attributed to surface
c aps         3,msul        array of coefficents
c keypa       mpat          keys to array of  coefficents
c kfio                      i/o unit number associated to the file of keys
c ifio                      i/o unit number associated to the file of coefficents
c ierr                      error condition indicator (see DIAGNOSTICS)
c
c * DIAGNOSTICS*
c
c  ierr =0                  no error
c  ierr =1                  n. of patchess  exceeds mpat
c  ierr =2                  n. of   surface coefficents exceeds msel
c  ierr =3                  ierr=1 + ierr=2
c
c  * NOTE *
c
c  if ierr/=0 a message is sent to the unit associated to ito. In standard FORTRAN
c  if you want the message on the standard output you should put ito=6
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine loas2d(keypa,aps,msul,mpat,isur,igsur,ib,n,m,
     1                  kfio,ifio,ips,ito,ierr)
      implicit none
      integer*8 ax_dummy,ib,ifio, ierr,igsur,i,ii,if_dummy,ips,isur
      integer*8 ito,j,jj,key,keypa(mpat),kfio,m,mpat,msul,n,npa
      integer*8 ndu,ndv
      real*8 aps(3,msul),ax,ay,az,xl
      ierr =0
      ips=0
      read(kfio,rec=isur)key,n,m,igsur,ib
      if(n.eq.0)then
C
C DEGENERATE SURFACE
C
       read(ifio,rec=key)aps(1,1),aps(2,1),aps(3,1)
       ips =1
       return
      endif
c
      key = key-1
      npa = 0
      do i=1,n-1
      do j=1,m-1
         npa=npa+1
         ips = ips +1
         if(npa.le.mpat)keypa(npa)=ips
         read(ifio,rec=key+ips)ndu,ndv,xl
         aps(1,ips)=dble(ndu)
         aps(2,ips)=dble(ndv)
         aps(3,ips)=xl
         do jj=1,ndv
         do ii=1,ndu
            ips = ips+1
            if(ips.gt.msul)then
              exit
            endif
            read(ifio,rec=key+ips)ax,ay,az
            aps(1,ips)=ax
            aps(2,ips)=ay
            aps(3,ips)=az
         end do
         end do
      end do
      end do
      if(npa.gt.mpat)then
        write(ito,*)' ERROR LOAS2D : 001'
        write(ito,*)' too many patches on surface'
        write(ito,*)' global and local surface number:',isur,igsur
        write(ito,*)' reading aborted'
        write(ito,*)' increase max n. of patches at least to',npa
        ierr =1
      endif
      if(ips.gt.msul)then
        write(ito,*)' ERROR LOAS2D : 002'
        write(ito,*)' too many coefficents on surface'
        write(ito,*)' global and local surface number:',isur,igsur
        write(ito,*)' reading aborted'
        write(ito,*)' increase MSUL at least to',ips
        ierr =1 + ierr
      endif
      return
      end
C
C=END SOURCE
