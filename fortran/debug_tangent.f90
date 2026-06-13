!===============================================================================
! Test to Debug Tangent Vector Computation
!===============================================================================

program debug_tangent
  implicit none
  integer, parameter :: n = 5
  real(kind=8), parameter :: pi = 3.14159265358979323846d0
  real(kind=8) :: q(2, n), t(2, n), cs(n-1)
  real(kind=8) :: a(n), b(n), c(n), timp(4)
  real(kind=8) :: theta
  integer :: i, ispty(2)

  write(*,*) "====== DEBUG: Tangent Vector Computation ======"
  write(*,*)

  ! Create points on quarter circle
  do i = 1, n
    theta = (pi/2.0d0) * dble(i-1) / dble(n-1)
    q(1, i) = cos(theta)
    q(2, i) = sin(theta)
  end do

  ! Compute chord lengths
  call cholen(2_8, int(n, 8), q, cs)

  write(*,*) "Chord lengths computed:"
  do i = 1, n-1
    write(*,'(A,I2,A,F15.10)') "  cs(", i, ") = ", cs(i)
  end do
  write(*,*)

  ! Try natural spline (ispty = 0)
  write(*,*) "Testing NATURAL boundary conditions (ispty = 0)"
  ispty(1) = 0
  ispty(2) = 0

  write(*,*) "Before evtan, tangent array:"
  do i = 1, n
    write(*,'(A,I2,A,F12.8,A,F12.8)') "  t(:,", i, ") = (", t(1,i), ", ", t(2,i), ")"
  end do
  write(*,*)

  write(*,*) "Calling evtan..."
  call evtan(2_8, int(n, 8), q, cs, a, b, c, t, ispty, timp)

  write(*,*)
  write(*,*) "After evtan, tangent array:"
  do i = 1, n
    if (isnan(t(1,i)) .or. isnan(t(2,i))) then
      write(*,'(A,I2,A)') "  t(:,", i, ") = NaN detected!"
    else if (t(1,i) > 1.0d10 .or. t(2,i) > 1.0d10) then
      write(*,'(A,I2,A)') "  t(:,", i, ") = Very large value (likely Inf)"
    else
      write(*,'(A,I2,A,F12.8,A,F12.8)') "  t(:,", i, ") = (", t(1,i), ", ", t(2,i), ")"
    end if
  end do

end program debug_tangent

function isnan(x)
  real(kind=8), intent(in) :: x
  logical :: isnan
  isnan = x /= x
end function isnan
