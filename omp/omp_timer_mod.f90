!*******************************************************************************
!> Support tic/toc timing, as in MATLAB, and accurate wallclock() function.
!> The timing is generally much more accurate if this module is compiled with
!> OMP active.
!*******************************************************************************
MODULE omp_timer_mod
  implicit none
  private
  real(8), external:: omp_get_wtime
  type omp_timer_t
  contains
    procedure, nopass:: delay
    procedure, nopass:: get
    procedure, nopass:: set
  end type
  type(omp_timer_t), public:: omp_timer
  real(8), save:: offset=0d0
  public wallclock
CONTAINS

!===============================================================================
FUNCTION wallclock() result (time)
  real(8):: time
  !.............................................................................
  time = omp_get_wtime()
  if (offset == 0d0) then
    offset=time
  end if
  time = time-offset
END FUNCTION wallclock

!===============================================================================
!> Active spin delay
!===============================================================================
SUBROUTINE delay (delta)
  real:: delta
  real(8):: wc
  !.............................................................................
  wc = wallclock()
  do while (wallclock()-wc < delta)
  end do
END SUBROUTINE delay

!===============================================================================
!> Get offset
!===============================================================================
SUBROUTINE get (wc)
  real(8):: wc
  !.............................................................................
  wc = offset
END SUBROUTINE get

!===============================================================================
!> Set offset
!===============================================================================
SUBROUTINE set (wc)
  real(8):: wc
  !.............................................................................
  offset = wc
END SUBROUTINE set

END MODULE omp_timer_mod
