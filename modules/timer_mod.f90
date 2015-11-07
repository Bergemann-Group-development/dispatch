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

END MODULE timer_mod
