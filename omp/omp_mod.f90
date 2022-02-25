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
    integer:: place=0
    integer:: thread_global=0, nthreads_global=0
    integer:: ncores=1
    integer:: nthreads=1
    integer:: nplaces=0
    integer:: nsockets=1
    integer:: thread_task=0
    logical:: trace=.false.
    logical:: master=.true.
    real(8):: wait_time=0d0
  contains
    procedure:: init
    procedure:: info
    procedure, nopass:: mythread
    procedure, nopass:: numthreads => omp_numthreads
    procedure, nopass:: set_stacksize
    procedure, nopass:: delay
  end type
  public:: omp_numthreads
  type(omp_t), public, save:: omp
  !$omp threadprivate (omp)
CONTAINS

!===============================================================================
SUBROUTINE init (self)
  class(omp_t):: self
  integer:: thread, nthreads, place, nplaces
  logical:: master
  !-----------------------------------------------------------------------------
  !print*,'OMP threads:'
  !$omp parallel private(master,thread,place)
#if defined (_OPENMP)
  !$omp single
  !call omp_set_nested(.true.)
  call omp_set_max_active_levels(2)
  nthreads = omp_get_num_threads()
  nplaces = omp_get_num_places()
  !$omp end single
  thread = omp_get_thread_num()
  master = (thread == 0)
  if (nplaces > 0) then
    place = omp_get_place_num()
  else
    place = thread
  end if
#else
  nthreads = 1
  master = .true.
  thread = 0
  place = 0
  nplaces = 0
#endif
  omp%nthreads = nthreads
  omp%thread = thread
  omp%place = place
  omp%nplaces = nplaces
  omp%master = master
  io_unit%master = io_unit%master .and. master
  !$omp end parallel
  omp%nthreads = nthreads
  omp%nplaces = nplaces
END SUBROUTINE init

!===============================================================================
!> Write information about threads and places
!===============================================================================
SUBROUTINE info (self)
  class(omp_t):: self
  integer:: thread
  !-----------------------------------------------------------------------------
  !$omp single
  write (io_unit%mpi,'(a,2i4)') &
    ' omp%init: nthreads,nplaces=', omp%nthreads, omp%nplaces
  !$omp end single
  do thread=0,self%nthreads-1
    !$omp barrier
    if (thread==omp%thread) then
      write (io_unit%mpi,'(a,2i4)') &
        '  omp%init: thread,place=', omp%thread, omp%place
    end if
  end do
END SUBROUTINE info

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
SUBROUTINE set_stacksize (nw, verbose)
#ifdef __INTEL_COMPILER
#ifdef _OPENMP
  USE omp_lib, only : kmp_set_stacksize_s, kmp_get_stacksize_s, kmp_size_t_kind
  integer(kind=kmp_size_t_kind) :: stacksize, new_stacksize
  real:: r
#endif
#endif
  integer:: nw, recommended
  logical:: verbose
  character(len=16):: envvar
  logical:: kmp, omp
  logical, save:: first_time=.true., printed=.false.
  !-----------------------------------------------------------------------------
  if (io_unit%do_validate) return
  recommended = 6*(4*nw)
#ifdef _OPENMP
#ifdef __INTEL_COMPILER
  call getenv ('KMP_STACKSIZE', envvar)
  kmp = envvar /= ''
  call getenv ('OMP_STACKSIZE', envvar)
  omp = envvar /= ''
  if (omp .and. .not.kmp .and. omp_get_num_threads()==0 .and. first_time) then
    if (verbose .and. .not.printed) print '(/2x,a)', &
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
    if (first_time .and. omp_get_num_threads()==0) then
      if (verbose .and. .not.printed) print '(1x,a)', &
        '*************************************************************************************'
      print '(1x,a,i5,a)', &
      '* WARNING! OpenMP stacksize has been reset to ', &
      nint(new_stacksize/1024.**2), ' MB to avoid segfaults'
      printed = .true.
    end if
  else if (first_time .and. omp_get_num_threads()==0) then
      print '(1x,a,f5.1,a,i8,a)', 'OpenMP stacksize (default):', &
      stacksize/1024.**2, ' MB for', nw, ' mem-words'
    first_time = .false.
  endif
#else
  call getenv ('OMP_STACKSIZE', envvar)
  if (verbose .and. envvar == '' .and. omp_get_num_threads()==0 .and. first_time) then
    print '(/,2(1x,a,/),1x,a,i9)', &
    '*************************************************************************************', &
    '* omp_mod.f90:', &
    '* WARNING: environment variable OMP_STACKSIZE not set, recommended value =',recommended
    call getenv ('KMP_STACKSIZE', envvar)
    if (envvar /= '' .and. omp_get_num_threads()==0 .and. first_time) then
      if (.not.printed) print *, &
        '*************************************************************************************'
      print '(2(1x,a,/),1x,a)', &
        '* omp_mod.f90:', &
        '* WARNING: environment variable KMP_STACKSIZE is set, this assumes Intel MPI'
    end if
    printed = .true.
  end if
#endif
#endif
  if (omp_get_num_threads()==0 .and. first_time) then
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
