! $Id$
!> =============================================================================
MODULE omp_mod
  USE omp_lib
  implicit none
  private
  logical, save:: omp_trace=.false.
  type, public:: omp_t
    integer:: nthreads=1
    logical:: trace=.false.
    logical:: master=.true.
  contains
    procedure:: init
    procedure, nopass:: mythread => omp_mythread
    procedure, nopass:: numthreads => omp_numthreads
  end type
  public:: init_omp, omp_mythread, omp_numthreads
  logical, public, save:: omp_master
  integer, public, save:: omp_nthreads
  type(omp_t), public, save:: omp
  !$omp threadprivate (omp_master, omp)
CONTAINS

!> =============================================================================
SUBROUTINE init (self)
  class(omp_t):: self
  call init_omp
  !$omp single
  self%nthreads = omp_nthreads
  !$omp end single
  self%master = (omp_get_thread_num() == 0)
END SUBROUTINE init

!> =============================================================================
SUBROUTINE init_omp
!...............................................................................
#if defined (_OPENMP)
!$omp parallel
!$omp single
  omp_nthreads = omp_get_num_threads()
  omp_master = (omp_get_thread_num() == 0)
!$omp end single
!$omp end parallel
#else
  omp_nthreads = 1
  omp_master = .true.
#endif
END SUBROUTINE init_omp

#if ! defined (_OPENMP)
#if ! defined (__xlc__)
!> =============================================================================
LOGICAL FUNCTION omp_in_parallel()
  omp_in_parallel = .false.
END FUNCTION
#endif
#endif

!> =============================================================================
INTEGER FUNCTION omp_mythread()
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
INTEGER FUNCTION omp_numthreads()
#if defined (_OPENMP)
  if (omp_in_parallel()) then
    omp_numthreads = omp_get_num_threads()
  else
    omp_numthreads = 1
  end if
#else
  omp_numthreads = 1
#endif
END FUNCTION

!> =============================================================================
END MODULE
