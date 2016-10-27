!> *****************************************************************************
!> Support tic/toc timing, as in MATLAB, and accurate wallclock() function.
!> The timing is generally much more accurate if this module is compiled with
!> OMP active.
!> *****************************************************************************
MODULE omp_timer_mod
  implicit none
  private
  real(8), external:: omp_get_wtime
  type omp_timer_t
  contains
    procedure, nopass:: wallclock
    procedure, nopass:: delay
  end type
  type(omp_timer_t), public:: omp_timer
  public wallclock
CONTAINS

!> =============================================================================
FUNCTION wallclock() result (delta)
  real(8):: time, delta
  real(8), save:: offset=-1.
  !.............................................................................
  time = omp_get_wtime()
  if (offset == -1.) then
    offset=time
  end if
  delta = time-offset
END FUNCTION wallclock

!> =============================================================================
!> Active spin delay
!> =============================================================================
SUBROUTINE delay (delta)
  real:: delta
  real(8):: wc
  !.............................................................................
  wc = wallclock()
  do while (wallclock()-wc < delta)
  end do
END SUBROUTINE delay

END MODULE omp_timer_mod
