!===============================================================================
! Detailed Debugging Test for Arc Length Computation
!===============================================================================

program debug_arc_length
  implicit none
  integer, parameter :: n = 3  ! Simple 2-segment line
  real(kind=8) :: q(2, n), t(2, n), cs(n-1), len(n-1)
  real(kind=8) :: a(n), b(n), c(n), timp(4)
  real(kind=8) :: expected_length, total_arc
  integer :: i, ispty(2)

  write(*,*) "====== DEBUG: Arc Length Computation ======"
  write(*,*)

  ! Test 1: Straight line from (0,0) to (2,0)
  write(*,*) "Test 1: Straight horizontal line (0,0) -> (2,0)"
  write(*,*) "Expected length: 2.0"
  write(*,*)

  ! Set up points on a straight line
  q(1, 1) = 0.0d0
  q(2, 1) = 0.0d0

  q(1, 2) = 1.0d0
  q(2, 2) = 0.0d0

  q(1, 3) = 2.0d0
  q(2, 3) = 0.0d0

  ! Compute spline with natural boundary conditions
  ispty(1) = 0
  ispty(2) = 0

  write(*,*) "Calling psplin..."
  call psplin(2_8, int(n, 8), q, t, cs, len, ispty, a, b, c)

  ! Print chord lengths
  write(*,*)
  write(*,*) "Chord lengths (cs):"
  do i = 1, n-1
    write(*,'(A,I2,A,F15.10)') "  cs(", i, ") = ", cs(i)
  end do

  ! Print tangent vectors
  write(*,*)
  write(*,*) "Tangent vectors (t):"
  do i = 1, n
    write(*,'(A,I2,A,F12.6,A,F12.6,A)') "  t(:,", i, ") = (", t(1,i), ", ", t(2,i), ")"
  end do

  ! Print computed arc lengths
  write(*,*)
  write(*,*) "Arc lengths (len):"
  total_arc = 0.0d0
  do i = 1, n-1
    write(*,'(A,I2,A,F15.10)') "  len(", i, ") = ", len(i)
    total_arc = total_arc + len(i)
  end do

  write(*,*)
  write(*,'(A,F15.10)') "Total arc length:    ", total_arc
  write(*,'(A,F15.10)') "Expected:            ", 2.0d0
  write(*,'(A,F15.10)') "Error:               ", abs(total_arc - 2.0d0)

  ! Test 2: Vertical line
  write(*,*)
  write(*,*) "====== Test 2: Vertical line (0,0) -> (0,3) ======"
  write(*,*) "Expected length: 3.0"
  write(*,*)

  q(1, 1) = 0.0d0
  q(2, 1) = 0.0d0

  q(1, 2) = 0.0d0
  q(2, 2) = 1.5d0

  q(1, 3) = 0.0d0
  q(2, 3) = 3.0d0

  call psplin(2_8, int(n, 8), q, t, cs, len, ispty, a, b, c)

  write(*,*) "Chord lengths (cs):"
  do i = 1, n-1
    write(*,'(A,I2,A,F15.10)') "  cs(", i, ") = ", cs(i)
  end do

  write(*,*) "Tangent vectors (t):"
  do i = 1, n
    write(*,'(A,I2,A,F12.6,A,F12.6,A)') "  t(:,", i, ") = (", t(1,i), ", ", t(2,i), ")"
  end do

  write(*,*) "Arc lengths (len):"
  total_arc = 0.0d0
  do i = 1, n-1
    write(*,'(A,I2,A,F15.10)') "  len(", i, ") = ", len(i)
    total_arc = total_arc + len(i)
  end do

  write(*,*)
  write(*,'(A,F15.10)') "Total arc length:    ", total_arc
  write(*,'(A,F15.10)') "Expected:            ", 3.0d0
  write(*,'(A,F15.10)') "Error:               ", abs(total_arc - 3.0d0)

end program debug_arc_length
