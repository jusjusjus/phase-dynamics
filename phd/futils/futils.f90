
MODULE futils
IMPLICIT NONE
  ! PRIVATE

  double precision :: p  ! Period
  double precision :: p2 ! Period / 2

  public :: unmod
  public :: set_period
  public :: threshold_data
  public :: poincare_times_no_false_returns
  public :: poincare_times_interpolate


CONTAINS

  SUBROUTINE set_period( new_period )
  double precision, INTENT(in) :: new_period
    p = new_period
    p2 = 0.5D0 * p
  END SUBROUTINE set_period


  SUBROUTINE unmod (x, x_mod, num_samples)
  integer, INTENT(in) :: num_samples
  double precision, dimension(num_samples), INTENT(in) :: x_mod
  double precision, dimension(num_samples), INTENT(out) :: x
    integer :: i
    integer :: n
    double precision :: ph_diff
    n = 0
    x = x_mod
    DO i = 2, num_samples
      x(i) = x(i) + p * n
      ph_diff = x(i) - x(i-1)
      IF (ph_diff > p2) THEN
          n = n-1
          x(i) = x(i) - p
      ELSE IF (ph_diff < -p2) THEN
          n = n+1
          x(i) = x(i) + p
      END IF
    END DO
  END SUBROUTINE


  SUBROUTINE threshold_data (x, num_samples, threshold, n_min, segments, num_segments)
  double precision, dimension(num_samples), intent(in) :: x
  integer, intent(in) :: num_samples
  double precision, intent(in) :: threshold
  integer, intent(in) :: n_min
  integer, dimension(num_samples, 2), intent(out) :: segments
  integer, intent(out) :: num_segments
    integer :: n, n0
    num_segments = 0
    n0 = 1
    DO n = 1, num_samples
      IF (x(n) < threshold) THEN
        IF (n-n0 > n_min) THEN
          num_segments = num_segments + 1
          segments(num_segments, 1) = n0
          segments(num_segments, 2) = n-1
        END IF
        n0 = n
      END IF
    END DO
    IF (n-n0 > n_min) THEN
      num_segments = num_segments + 1
      segments(num_segments, 1) = n0
      segments(num_segments, 2) = n-1
    END IF
  END SUBROUTINE


  SUBROUTINE poincare_times_no_false_returns (cross_idx, ti, num_crossings, x_mod, x, num_samples)

  ! idx stores the indicies of the crossings
  ! ti stores the time indicies of return
  ! num_crossings stores the number of crossings
  ! x_mod is modded phase
  ! x is unmodded phase
  ! num_samples is length of x
  ! function stores positive times when x crosses 0 

  integer, dimension(num_samples), intent(in) :: x_mod, x
  integer, intent(in) :: num_samples
  integer, INTENT(out) :: num_crossings
  integer, dimension(num_samples), intent(out) :: cross_idx
  integer, dimension(num_samples), intent(out) :: ti
    integer :: i, last_crossing = 0
    double precision :: dx_mod, dx
    double precision :: threshold
    num_crossings = 1 ! number of crossings
    threshold = -p2


    DO i = 2, num_samples
    ! find first positive crossing of phase, x = 0
      dx_mod = x_mod(i)-x_mod(i-1)
      IF (dx_mod < threshold) THEN
      ! only computes positiv crossings
        cross_idx(num_crossings) = i-1 ! saves value before the crossing
        last_crossing = cross_idx(num_crossings)
        EXIT
      END IF
    END DO

    DO i = i+1, num_samples
    ! ... now for the rest
      dx_mod = x_mod(i)-x_mod(i-1)
      dx     = x(i-1)-x(last_crossing)
      IF (dx_mod < threshold .AND. dx > -threshold) THEN
      ! only computes positiv crossings, and only if one oscillation was performed
        ti(num_crossings) = i-1-last_crossing
        num_crossings = num_crossings + 1
        cross_idx(num_crossings) = i-1 ! saves value before the crossing
        last_crossing = cross_idx(num_crossings)
      END IF
    END DO
  END SUBROUTINE poincare_times_no_false_returns


  SUBROUTINE poincare_times_interpolate (idx, ti, num_crossings, x, num_samples, idx_interp, ti_interp)

  ! ti(i) is the discrete time between two crossings
  ! idx(i) is the index of the step before 0
  ! x is phase, x[idx(i)] > 0.0
  ! times are the output where we will store the correct return times (= approx index)
  ! n: number of approximate zeros that were found

  integer, dimension(num_crossings), intent(in) :: idx, ti
  integer, intent(in) :: num_crossings
  double precision, dimension(num_samples), intent(in) :: x
  integer, intent(in) :: num_samples
  double precision, dimension(num_crossings), intent(out) :: idx_interp, ti_interp

    integer :: i, j, j1, JJ, JJ1
    double precision :: x_j, x_j1, x_jj, x_jj1
    DO i = 1, num_crossings
      j   = idx(i)     ! index before crossing
      j1  = idx(i)+1   ! index after crossing
      JJ  = j+ti(i)    ! index before return
      JJ1 = j+ti(i)+1  ! index after return

      !! Corresponding phases
      ! Before crossing are close to 'p'
      x_j   = x(j)-p
      x_jj  = x(JJ)-p
      ! After crossing are just above zero
      x_j1  = x(j1)
      x_jj1 = x(JJ1)

      idx_interp(i) = idx(i)- x_j/(x_j1-x_j)
      ti_interp(i)  = ti(i) + x_j/(x_j1-x_j)-x_jj/(x_jj1-x_jj)
    END DO
  END SUBROUTINE poincare_times_interpolate

END MODULE futils
