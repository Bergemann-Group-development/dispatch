!> *****************************************************************************
!> Support tic/toc timing, as in MATLAB, and accurate wallclock() function.
!> The timing is generally much more accurate if this module is compiled with
!> OMP active.
!> *****************************************************************************
MODULE omp_timer_mod
  implicit none
  real(8), save:: wt=0d0
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE tic
  implicit none
  !.............................................................................
  real(8):: omp_get_wtime
  wt = omp_get_wtime()
END SUBROUTINE tic

!> =============================================================================
!> =============================================================================
SUBROUTINE toc (label, n)
  implicit none
  character(len=*):: label
  integer:: n
  real(8):: omp_get_wtime
  !.............................................................................
  wt = (omp_get_wtime() - wt)
  if (wt/n < 1d-6) then
    print 1, label, wt,  'sec', wt/n*1e9, 'ns/cell'
  else if (wt/n < 1d-3) then
    print 1, label, wt,  'sec', wt/n*1e6, 'mus/cell'
  else if (wt/n < 1d0) then
    print 1, label, wt,  'sec', wt/n*1e3, 'ms/cell'
  else
    print 1, label, wt,  'sec', wt/n, 's/cell'
    1 format(1x,a,1h:,2(g20.3,1x,a))
  end if
END SUBROUTINE toc

!> =============================================================================
FUNCTION wallclock() result (delta)
  real(kind=8):: time, delta
  real(kind=8), save:: offset=-1.
  real(8):: omp_get_wtime
!...............................................................................
  time = omp_get_wtime()
  if (offset == -1.) then
    offset=time
  end if
  delta = time-offset
END FUNCTION wallclock

END MODULE omp_timer_mod
