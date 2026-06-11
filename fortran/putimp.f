      subroutine putimp(ndimn,timp,nn,t,nl)
      implicit none
      integer*8 id,ndimn,nn,nl
      real*8 t(ndimn,*),timp(ndimn,*)
      do 10 id=1,ndimn
10        timp(id,nn)=t(id,nl)
      return
      end
