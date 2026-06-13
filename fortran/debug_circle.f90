!===============================================================================
! Detailed Circle Debugging Test
!===============================================================================

program debug_circle
  implicit none
  integer, parameter :: n = 5  ! Simpler: 5 points on quarter circle
  real(kind=8), parameter :: pi = 3.14159265358979323846d0
  real(kind=8) :: q(2, n), t(2, n), cs(n-1), len(n-1)
  real(kind=8) :: a(n), b(n), c(n), timp(4)
  real(kind=8) :: theta, total_arc, expected_arc
  real(kind=8) :: chord_sum
  integer :: i, ispty(2)

  write(*,*) "====== DEBUG: Circle Arc Length ======"
  write(*,*)

  ! Test: Quarter circle (0 to pi/2, radius = 1)
  write(*,*) "Quarter circle, radius = 1.0"
  write(*,'(A,F15.10)') "Expected arc length: ", pi/2.0d0
  write(*,*)

  ! Create points on unit circle for quarter arc
  do i = 1, n
    theta = (pi/2.0d0) * dble(i-1) / dble(n-1)
    q(1, i) = cos(theta)
    q(2, i) = sin(theta)
    write(*,'(A,I2,A,F10.6,A,F10.6,A)') "Point ", i, ": (", q(1,i), ", ", q(2,i), ")"
  end do

  write(*,*)
  write(*,*) "Computing spline with natural boundary conditions..."

  ! Compute spline
  ispty(1) = 0
  ispty(2) = 0

  call psplin(2_8, int(n, 8), q, t, cs, len, ispty, a, b, c)

  ! Print chord lengths
  write(*,*)
  write(*,*) "Chord lengths (cs):"
  chord_sum = 0.0d0
  do i = 1, n-1
    chord_sum = chord_sum + cs(i)
    write(*,'(A,I2,A,F15.10)') "  cs(", i, ") = ", cs(i)
  end do
  write(*,'(A,F15.10)') "  Sum of chords:      ", chord_sum

  ! Print tangent vectors
  write(*,*)
  write(*,*) "Tangent vectors (t):"
  do i = 1, n
    write(*,'(A,I2,A,F12.8,A,F12.8,A,F12.8)') "  t(:,", i, ") = (", t(1,i), ", ", t(2,i), "), mag=", &
        sqrt(t(1,i)**2+t(2,i)**2)
  end do

  ! Print computed arc lengths
  write(*,*)
  write(*,*) "Arc lengths (len):"
  total_arc = 0.0d0
  do i = 1, n-1
    write(*,'(A,I2,A,F15.10,A,F15.10)') "  len(", i, ") = ", len(i), ", ratio to chord = ", len(i)/cs(i)
    total_arc = total_arc + len(i)
  end do

  write(*,*)
  write(*,'(A,F15.10)') "Total arc length:    ", total_arc
  write(*,'(A,F15.10)') "Expected:            ", pi/2.0d0
  write(*,'(A,F15.10)') "Error:               ", abs(total_arc - pi/2.0d0)
  write(*,'(A,F15.10)') "Relative error:      ", abs(total_arc - pi/2.0d0) / (pi/2.0d0)
  write(*,'(A,F15.10)') "Ratio to chord sum:  ", total_arc / chord_sum

end program debug_circle
