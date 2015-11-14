!> *****************************************************************************
!> Support tic/toc timing, as in MATLAB, and accurate wallclock() function.
!> The timing is generally much more accurate if this module is compiled with
!> OMP active.
!> *****************************************************************************
MODULE mpi_timer
  implicit none
#ifndef _OPENMP
  include 'mpif.h'
#endif
  real(8), save:: wt=0d0
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE tic
  implicit none
  !.............................................................................
#ifdef _OPENMP
  real(8):: omp_get_wtime
  wt = omp_get_wtime()
#else
  real(8):: mpi_wtime
  wt = MPI_wtime()
#endif
END SUBROUTINE tic

!> =============================================================================
!> =============================================================================
SUBROUTINE toc (label, n)
  implicit none
  character(len=*):: label
  integer:: n
  !.............................................................................
#ifdef _OPENMP
  real(8):: omp_get_wtime
  wt = (omp_get_wtime() - wt)
#else
  real(8):: mpi_wtime
  wt = (MPI_wtime() - wt)
#endif
  if (wt/n < 1d-6) then
    print 1, label, wt,  'sec', wt/n*1e9, 'ns/call'
  else if (wt < 1d-3) then
    print 1, label, wt,  'sec', wt/n*1e6, 'mus/call'
  else if (wt < 1d0) then
    print 1, label, wt,  'sec', wt/n*1e3, 'ms/call'
  else
    print 1, label, wt,  'sec', wt/n, 's/call'
    1 format(1x,a,1h:,2(g20.3,1x,a))
  end if
END SUBROUTINE toc

!> =============================================================================
FUNCTION wallclock() result (delta)
  real(kind=8):: time, delta
  real(kind=8), save:: offset=-1.
!...............................................................................
#ifdef _OPENMP
  real(kind=8):: OMP_Get_Wtime, OMP_Get_Wtick
  logical:: OMP_in_parallel
  if (OMP_in_parallel()) then
    time = OMP_Get_Wtime()
  else
    !$omp parallel
    time = OMP_Get_Wtime()
    !$omp end parallel
  end if
  if (offset == -1.) then
    !print *,'tick:', OMP_Get_Wtick()
#else
  time = MPI_Wtime()
  if (offset == -1.) then
    !print *,'tick:', MPI_Wtick()
#endif
    offset=time
  end if
  delta = time-offset
END FUNCTION wallclock

END MODULE mpi_timer
