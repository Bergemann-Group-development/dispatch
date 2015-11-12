! $Id$
!> =============================================================================
MODULE omp_mod
  implicit none
  integer, save:: omp_nthreads
  logical, save:: omp_trace=.false.
PUBLIC
CONTAINS

#if ! defined (_OPENMP)
#if ! defined (__xlc__)
!> =============================================================================
LOGICAL FUNCTION omp_in_parallel()
  implicit none
  omp_in_parallel = .false.
END FUNCTION
#endif
#endif

!> =============================================================================
INTEGER FUNCTION omp_mythread()
  implicit none
  integer omp_get_thread_num
  logical omp_in_parallel
#if defined (_OPENMP)
  if (omp_in_parallel()) then
    omp_mythread = omp_get_thread_num()
  else
    omp_mythread = 0
  end if
#else
  omp_mythread = 0
#endif
END FUNCTION

!> =============================================================================
SUBROUTINE init_omp
  implicit none
  integer omp_get_num_threads
!...............................................................................
#if defined (_OPENMP)
!$omp parallel
  omp_nthreads = omp_get_num_threads()
!$omp end parallel
#else
  omp_nthreads = 1
#endif
END SUBROUTINE init_omp

!> =============================================================================
END MODULE
