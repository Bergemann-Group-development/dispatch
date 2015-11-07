!> *****************************************************************************
!> Support tic/toc timing, as in MATLAB
!> *****************************************************************************
MODULE timer_mod
  implicit none
#ifndef _OPEMP
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
  wt = (omp_get_wtime() - wt)/n
#else
  real(8):: mpi_wtime
  wt = (MPI_wtime() - wt)/n
#endif
  if (wt < 1d-3) then
    print 1, label, wt*1e6, 'mus/call'
  else if (wt < 1d0) then
    print 1, label, wt*1e3, 'ms/call'
  else
    print 1, label, wt, 's/call'
    1 format(1x,a,g11.3,1x,a)
  end if
END SUBROUTINE toc

END MODULE timer_mod
