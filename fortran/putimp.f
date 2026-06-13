      subroutine putimp(ndimn,timp,nn,t,nl)
      implicit none
      integer*8 id,ndimn,nn,nl
      real*8 t(ndimn,*),timp(ndimn,*)
      do id=1,ndimn
        timp(id,nn)=t(id,nl)
      end do
      return
      end
