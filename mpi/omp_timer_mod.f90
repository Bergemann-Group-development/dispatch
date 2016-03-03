!> *****************************************************************************
!> Support tic/toc timing, as in MATLAB, and accurate wallclock() function.
!> The timing is generally much more accurate if this module is compiled with
!> OMP active.
!> *****************************************************************************
MODULE omp_timer_mod
  USE omp_mod
  implicit none
  real(8), save:: wt=0d0
  interface toc
    module procedure toc_i4, toc_i8
  end interface
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
SUBROUTINE toc_i4 (label, n)
  implicit none
  character(len=*):: label
  integer:: n
  !.............................................................................
  call toc_r (label, real(n))
END SUBROUTINE toc_i4

!> =============================================================================
!> =============================================================================
SUBROUTINE toc_i8 (label, n)
  implicit none
  character(len=*):: label
  integer(8):: n
  !.............................................................................
  call toc_r (label, real(n))
END SUBROUTINE toc_i8

!> =============================================================================
!> =============================================================================
SUBROUTINE toc_r (label, n)
  implicit none
  character(len=*):: label
  real:: n
  real(8):: omp_get_wtime
  !.............................................................................
  wt = (omp_get_wtime() - wt)*omp_numthreads()
  if (wt/n < 1d-6) then
    print 1, label, wt,  'sec', wt/n*1e9, 'ns/upd,', n, 'updates'
  else if (wt/n < 1d-3) then
    print 1, label, wt,  'sec', wt/n*1e6, 'mus/upd,', n, 'updates'
  else if (wt/n < 1d0) then
    print 1, label, wt,  'sec', wt/n*1e3, 'ms/upd,', n, 'updates'
  else
    print 1, label, wt,  'sec', wt/n, 's/upd,', n, 'updates'
    1 format(1x,a16,1h:,2(g20.3,1x,a),1p,e10.2,1x,a)
  end if
  wt = omp_get_wtime()
END SUBROUTINE toc_r

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
