      subroutine putimp(ndimn,timp,nn,t,nl)
      real t(ndimn,*),timp(ndimn,*)
      do 10 id=1,ndimn
10        timp(id,nn)=t(id,nl)
      return
      end
