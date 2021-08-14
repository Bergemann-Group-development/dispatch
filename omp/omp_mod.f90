! $Id$
!===============================================================================
MODULE omp_mod
  USE omp_lib
  USE io_unit_mod
  implicit none
  private
  public:: omp_get_thread_num, omp_get_num_threads, omp_in_parallel
  logical, save:: omp_trace=.false.
  type, public:: omp_t
    integer:: thread=0
    integer:: thread_global=0, nthreads_global=0
    integer:: ncores=1
    integer:: nthreads=1
    integer:: nsockets=1
    integer:: thread_task=1
    logical:: trace=.false.
    logical:: master=.true.
  contains
    procedure:: init
    procedure, nopass:: mythread
    procedure, nopass:: numthreads => omp_numthreads
    procedure, nopass:: set_stacksize
    procedure, nopass:: delay
  end type
  logical, public, save:: omp_master
  integer, public, save:: omp_nthreads
  integer, public, save:: omp_mythread
  public:: omp_numthreads
  type(omp_t), public, save:: omp
  !$omp threadprivate (omp_master, omp, omp_mythread)
CONTAINS

!===============================================================================
SUBROUTINE init (self)
  class(omp_t):: self
  !-----------------------------------------------------------------------------
  !print*,'OMP threads:'
  !$omp parallel
#if defined (_OPENMP)
  omp_master = (omp_get_thread_num() == 0)
  omp_mythread = omp_get_thread_num()
  !$omp single
  omp_nthreads = omp_get_num_threads()
  !$omp end single
#else
  omp_nthreads = 1
  omp_master = .true.
  omp_mythread = 0
#endif
  omp%nthreads = omp_nthreads
  omp%thread = omp_mythread
  omp%master = omp_master
  io_unit%master = io_unit%master .and. omp_master
  !$omp end parallel
END SUBROUTINE init

#if ! defined (_OPENMP)
#if ! defined (__xlc__)
#if ! defined (__PGI)
!===============================================================================
LOGICAL FUNCTION omp_in_parallel()
  omp_in_parallel = .false.
END FUNCTION
#endif
#endif
#endif

!===============================================================================
INTEGER FUNCTION mythread()
#if defined (_OPENMP)
  if (omp_in_parallel()) then
    mythread = omp_get_thread_num()
  else
    mythread = 0
  end if
#else
  mythread = 0
#endif
END FUNCTION

!===============================================================================
INTEGER FUNCTION omp_numthreads()
#if defined (_OPENMP)
  if (omp_in_parallel()) then
    omp_numthreads = omp_get_num_threads()
  else
    omp_numthreads = omp%nthreads
  end if
#else
  omp_numthreads = 1
#endif
END FUNCTION

!===============================================================================
!> Adjust OMP_STACKSIZE according to number of passive scalars, with baseline
!> being 4M for MHD
!===============================================================================
SUBROUTINE set_stacksize (nv)
#ifdef __INTEL_COMPILER
#ifdef _OPENMP
  USE omp_lib, only : kmp_set_stacksize_s, kmp_get_stacksize_s, kmp_size_t_kind
  integer(kind=kmp_size_t_kind) :: stacksize, new_stacksize
  real:: r
#endif
#endif
  integer:: nv, recommended
  character(len=16):: envvar
  logical:: kmp, omp
  logical, save:: first_time=.true., printed=.false.
  !-----------------------------------------------------------------------------
  if (io_unit%do_validate) return
  recommended = (4*1024**2*nv) / 8
#ifdef _OPENMP
#ifdef __INTEL_COMPILER
  call getenv ('KMP_STACKSIZE', envvar)
  kmp = envvar /= ''
  call getenv ('OMP_STACKSIZE', envvar)
  omp = envvar /= ''
  if (omp .and. .not.kmp .and. io_unit%master .and. first_time) then
    if (.not.printed) print '(/2x,a)', &
      '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print '(2x,a)', &
      '! WARNING: Use KMP_STACKSIZE rather than OMP_STACKSIZE with Intel-MPI!   !'
    print '(2x,a,/)', &
      '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    printed = .true.
  end if
  stacksize = kmp_get_stacksize_s()
  new_stacksize = recommended
  if (kind(r)==kind(1.0E0_8)) new_stacksize = new_stacksize * 2
  if (new_stacksize > stacksize) then
    call kmp_set_stacksize_s(new_stacksize)
    if (io_unit%master .and. first_time) then
      if (.not.printed) print '(1x,a)', &
        '*************************************************************************************'
      print '(1x,a,f5.1,a)', &
      '* WARNING! OpenMP stacksize has been reset to ', new_stacksize/1024.**2, &
      ' MB to avoid segfault'
      printed = .true.
    end if
  else if (first_time.and.io_unit%master) then
    print '(1x,a,f5.1,a,i3,a)', 'OpenMP stacksize default:', stacksize/1024.**2, ' MB for', &
      nv, ' variables'

    first_time = .false.
  endif
#else
  call getenv ('OMP_STACKSIZE', envvar)
  if (envvar == '' .and. io_unit%master .and. first_time) then
    print '(/,2(1x,a,/),1x,a,i9)', &
    '*************************************************************************************', &
    '* omp_mod.f90:', &
    '* WARNING: environment variable OMP_STACKSIZE not set, recommended value =',recommended
    printed = .true.
    call getenv ('KMP_STACKSIZE', envvar)
    if (envvar /= '' .and. io_unit%master .and. first_time) then
      if (.not.printed) print *, &
        '*************************************************************************************'
      print '(2(1x,a,/),1x,a)', &
        '* omp_mod.f90:', &
        '* WARNING: environment variable KMP_STACKSIZE is set, this assumes Intel MPI'
      printed = .true.
    end if
  end if
#endif
#endif
  if (printed .and. io_unit%master .and. first_time) then
    first_time = .false.
    print '(1x,a)', &
    '*************************************************************************************'
  end if
END SUBROUTINE set_stacksize

!===============================================================================
!> Active spin delay in milliseconds
!===============================================================================
SUBROUTINE delay (ms)
  real:: ms
#if defined (_OPENMP)
  real(8):: wc
  !.............................................................................
  wc = omp_get_wtime()
  do while (omp_get_wtime()-wc < ms*1d-3)
  end do
#endif
END SUBROUTINE delay

!===============================================================================
END MODULE omp_mod
