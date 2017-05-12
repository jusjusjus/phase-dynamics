
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


  SUBROUTINE poincare_times_no_false_returns (idx, ti, n_idx, x_mod, x, num_samples)

  ! idx stores the indicies of the crossings
  ! ti stores the time indicies of return
  ! n_idx stores the number of crossings
  ! x_mod is modded phase
  ! x is phase (unmodded?)
  ! num_samples is num_samples of x
  ! function stores positive AND negative crossings of th0

  integer, dimension(num_samples), intent(in) :: x_mod, x
  integer, intent(in) :: num_samples
  integer, INTENT(out) :: n_idx
  integer, dimension(num_samples), intent(out) :: idx
  integer, dimension(num_samples), intent(out) :: ti
    integer :: i, previous_crossing = 0
    double precision :: dtht_mod, dtht
    n_idx = 1 ! number of crossings


    DO i = 2, num_samples
    ! find first positive crossing of phase, x = 0
      IF (x_mod(i-1)-x_mod(i) > 1.5*p2) THEN
      ! only computes positiv crossings
        idx(n_idx) = i-1 ! saves value before the crossing
        previous_crossing = idx(n_idx)
        EXIT
      END IF
    END DO

    DO i = i, num_samples
    ! ... now for the rest
      dtht_mod = x_mod(i-1)-x_mod(i)
      dtht = x(i-1)-x(previous_crossing)

      IF (dtht_mod > 1.5*p2 .AND. dtht > 1.5*p2) THEN
      ! only computes positiv crossings, and only if one oscillation was performed
        ti(n_idx) = i-1-previous_crossing
        n_idx = n_idx + 1
        idx(n_idx) = i-1 ! saves value before the crossing
        previous_crossing = idx(n_idx)
      END IF
    END DO
  END SUBROUTINE poincare_times_no_false_returns


  SUBROUTINE poincare_times_interpolate (idx, ti, num_zeros, x, num_samples, idx_interp, ti_interp)

  ! ti(i) is the discrete time between two crossings
  ! idx(i) is the index of the step before 0
  ! x is phase, x[idx(i)] > 0.0
  ! times are the output where we will store the correct return times (= approx index)
  ! n: number of approximate zeros that were found

  integer, dimension(num_zeros), intent(in) :: idx, ti
  integer, intent(in) :: num_zeros
  double precision, dimension(num_samples), intent(in) :: x
  integer, intent(in) :: num_samples
  double precision, dimension(num_zeros), intent(out) :: idx_interp, ti_interp

    integer :: i, j, j1, JJ, JJ1
    double precision :: x_j, x_j1, x_jj, x_jj1
    DO i = 1, num_zeros
      j1 = idx(i)+1
      j = idx(i)  ! starting index of one recurrence
      JJ1 = j+ti(i)+1
      JJ = j+ti(i) ! ending index of one recurrence

      x_j = x(j)
      ! corresponding phases, corrected for the jumps
      IF (x_j>p2) THEN
        x_j = x_j-p
      END IF
      x_j1 = x(j1)
      ! all phases are smaller then 2pi
      IF (x_j1>p2) THEN
        x_j1 = x_j1-p
      END IF
      x_jj = x(JJ)
      IF (x_jj>p2) THEN
        x_jj = x_jj-p
      END IF
      x_jj1 = x(JJ1)
      IF (x_jj1>p2) THEN
        x_jj1 = x_jj1 - p
      END IF
      idx_interp(i) = idx(i)- x_j/(x_j1-x_j)
      ti_interp(i)  = ti(i) + x_j/(x_j1-x_j)-x_jj/(x_jj1-x_jj)
    END DO
  END SUBROUTINE poincare_times_interpolate

END MODULE futils
