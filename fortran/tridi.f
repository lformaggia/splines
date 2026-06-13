      subroutine trif (n,a,b,c)
c
c
c decomposition
c
      implicit none
      integer*8 i,n,n1
      real*8 a(*),b(*),c(*)
c
      i=1
      c(i)=c(i)/b(i)
c
      n1=n-1
c
      do i=2,n1
        b(i)=b(i)-a(i)*c(i-1)
        c(i)=c(i)/b(i)
      end do
c
      i=n
      b(i)=b(i)-a(i)*c(i-1)
c
      return
      end
      subroutine tris (n,a,b,c,d)
c
      implicit none
      integer*8 i,ib,n,n1
      real*8 a(*),b(*),c(*),d(*)
c
c back substitution .. single d
c
      i=1
      d(i)=d(i)/b(i)
c
      do i=2,n
        d(i)=(d(i)-a(i)*d(i-1))/b(i)
      end do
c
      n1=n-1
      do ib=1,n1
        i=n-ib
        d(i)=d(i)-c(i)*d(i+1)
      end do
c
      return
      end
      subroutine TRIS2(n,m,a,b,c,d)
c
      implicit none
      integer*8 i,ib,id,m,n,n1
      real*8 a(*),b(*),c(*),d(m,*)
c
c back substitution .. multiple d
c
      i=1
      do id=1,m
        d(id,i)=d(id,i)/b(i)
      end do
c
      do i=2,n
      do id=1,m
        d(id,i)=(d(id,i)-a(i)*d(id,i-1))/b(i)
      end do
      end do
c
      n1=n-1
      do ib=1,n1
        i=n-ib
        do id=1,m
          d(id,i)=d(id,i)-c(i)*d(id,i+1)
        end do
      end do
c
      return
      end
