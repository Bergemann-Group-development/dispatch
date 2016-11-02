! $Id$
!> =============================================================================
MODULE omp_mod
  USE omp_lib
  USE io_unit_mod
  implicit none
  private
  public:: omp_get_thread_num, omp_get_num_threads, omp_in_parallel
  logical, save:: omp_trace=.false.
  type, public:: omp_t
    integer:: thread=0
    integer:: ncores=1
    integer:: nthreads=1
    integer:: nsockets=1
    logical:: trace=.false.
    logical:: master=.true.
  contains
    procedure:: init
    procedure, nopass:: mythread
    procedure, nopass:: numthreads => omp_numthreads
    procedure, nopass:: set_stacksize
  end type
  logical, public, save:: omp_master
  integer, public, save:: omp_nthreads
  integer, public, save:: omp_mythread
  public:: omp_numthreads
  type(omp_t), public, save:: omp
  !$omp threadprivate (omp_master, omp, omp_mythread)
CONTAINS

!> =============================================================================
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
  !print *, omp_mythread, omp%thread, omp_master, omp_nthreads
  !$omp end parallel
  call detect_cores
END SUBROUTINE init

#if ! defined (_OPENMP)
#if ! defined (__xlc__)
!> =============================================================================
LOGICAL FUNCTION omp_in_parallel()
  omp_in_parallel = .false.
END FUNCTION
#endif
#endif

!> =============================================================================
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
!> Adjust OMP_STACKSIZE according to number of passive scalars, with baseline
!> being 4M for MHD
!> =============================================================================
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
  logical, save:: first_time=.true., printed=.false.
  !-----------------------------------------------------------------------------
  recommended = (4*1024**2*nv) / 8
#ifdef _OPENMP
#ifdef __INTEL_COMPILER
  call getenv ('OMP_STACKSIZE', envvar)
  if (envvar /= '' .and. io_unit%master .and. first_time) then
    if (.not.printed) print '(1x,a)', &
      '*************************************************************************************'
    print '(1x,a)', &
      '* WARNING: environment variable OMP_STACKSIZE is set, this assumes non-Intel MPI'
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
  else if (first_time) then
    print '(1x,a,f5.1,a)', 'OpenMP stacksize default:', stacksize/1024.**2, ' MB'
    first_time = .false.
  endif
#else
  call getenv ('KMP_STACKSIZE', envvar)
  if (envvar /= '' .and. io_unit%master .and. first_time) then
    if (.not.printed) print '(1x,a)', &
      '*************************************************************************************'
    print '(1x,a)', &
      '* WARNING: environment variable KMP_STACKSIZE is set, this assumes Intel MPI'
    printed = .true.
  end if
  call getenv ('OMP_STACKSIZE', envvar)
  if (envvar == '' .and. io_unit%master .and. first_time) then
    print '(/,1x,a,/,1x,a,i9)', &
    '*************************************************************************************', &
    '* WARNING: environment variable OMP_STACKSIZE not set, recommended value =',recommended
    printed = .true.
  end if
#endif
#endif
  if (printed .and. io_unit%master .and. first_time) then
    first_time = .false.
    print '(1x,a)', &
    '*************************************************************************************'
  end if
END SUBROUTINE set_stacksize

!> =============================================================================
!> Attempt to determine the actual number of cores. If /proc/cpuinfo exists
!> one has to count the number of sockets, but observing hom many different 
!> "physical id" there are.  The number of cores per socket is in the fild
!> "cpu cores".
!> =============================================================================
SUBROUTINE detect_cores
  character(len=32):: file="/proc/cpuinfo"
  character(len=32):: line
  logical:: exist
  integer:: i, iostat, id, oid, n_socket, n_core
  type:: socket_t
    type(socket_t), pointer:: next => null()
    integer:: id
  end type
  type(socket_t), pointer:: p, o, head
  !.............................................................................
  nullify(head)
  n_socket = 0
  oid = -1
  inquire (file=file, exist=exist)
  if (exist) then
    open (file=file, unit=io_unit%os, form='formatted', status='old')
    do while (.true.)
      read (io_unit%os,'(a32)',iostat=iostat) line
      if (iostat/=0) then
        exit
      else
        i = index(line,':')
        if (line(1:8)=='cpu core') then
          read (line(i+1:),*) n_core
        else if (line(1:8)=='physical') then
          read (line(i+1:),*) id
          if (new(id)) then
            n_socket = n_socket + 1
          end if
        end if
      end if
    end do
    close (io_unit%os)
  end if
  if (n_socket>0) n_core = n_core*n_socket
  !$omp parallel
  omp%ncores = n_core
  !$omp end parallel
!-----------------------------------------------------------------------------
! For good measure, although totally insignificant, deallocate the list
!-----------------------------------------------------------------------------
  p => head
  do while (associated(p))
    o => p
    p => p%next
    deallocate (o)
  end do
CONTAINS
  logical function new(id)
    integer:: id
    p => head
    o => p
!---------------------------------------------------------------------------
! If the id already exists, return false
!---------------------------------------------------------------------------
    do while (associated(p))
      o => p
      if (p%id == id) then
        new = .false.
        return
      end if
      p => p%next
    end do
!---------------------------------------------------------------------------
! If the id did not exist, add it, and return true
!---------------------------------------------------------------------------
    allocate (p)
    p%id = id
    if (associated(o)) then
      o%next => p
    else
      head => p
    end if
    new = .true.
  end function
END SUBROUTINE detect_cores

!> =============================================================================
END MODULE omp_mod
