!===============================================================================
!> Basic MPI calls; initialization, end, barriers, asserts, wallclock, ...
!>
!> To make use of this module, make sure it is compiled together with your code,
!> using either the Makefile and the 'make' command, or compiling manually, with
!>
!> mpi = ../../../../../mpi             # or wherever
!> mpifort -c $mpi/mpi_mod.f90          # compiler mpi_mod.f90
!> mpifort -c $mpi/mpi_...              # add any other module your need
!> mpifort -c your_code.f90             # compile your code
!> mpifort *.o -o your_code.x           # link together into your_code.x
!>
!> In your_code.f90 you add lines such as:
!>
!> USE mpi_mod                          ! makes the mpi% object available
!> ...
!> call mpi%init                        ! initializes the mpi%object
!> ...
!> if (mpi%rank==0) print *, 'this will be printed only on MPI process 0'
!> ...
!> call mpi%end                         ! closes MPI
!>
!> Too see which variables and procedures are inside mpi%, just look below!
!>
!===============================================================================
MODULE mpi_mod
  USE omp_mod
  USE omp_timer_mod
  USE io_unit_mod
  USE trace_mod
  USE shared_mod
  USE mpi_lib
  implicit none
  private
  integer:: mpi_err             ! MPI error status
  logical:: mpi_master          ! true on MPI master (rank 0)
  logical:: mpi_trace=.false.   ! debug trace option
  logical:: master=.true.       ! true on OMP master thread on MPI master process
  integer:: mpi_size            ! number of MPI processes
  integer:: mpi_rank            ! process number (zero-based)
  integer:: output=2
  logical:: first_time = .true.
  !-----------------------------------------------------------------------------
  ! Object that support atomic incremental counting across MPI ranks.  The
  ! counter value %v reflects the actual counter value, but only at the very
  ! moment when it was upated, not at the moment when it happens to be accessed.
  !-----------------------------------------------------------------------------
  type, public:: mpi_counter_t
    integer:: window=0          ! MPI window
    real(8):: r=0               ! shared window value
    integer(8):: i=0            ! shared window value
    integer(8):: value=0        ! actual counter value, after update
  contains
    procedure, private:: init4
    procedure, private:: init8
    procedure, private:: init8d
    generic:: init => init4, init8, init8d
    procedure, private:: reset4
    procedure, private:: reset8
    procedure, private:: reset8d
    generic:: reset => reset4, reset8, reset8d
    procedure, private:: update4
    procedure, private:: update8
    procedure, private:: update8d
    generic:: update => update4, update8, update8d
  end type
  !-----------------------------------------------------------------------------
  ! Object that contains the most important MPI parameters
  !-----------------------------------------------------------------------------
  type, public:: mpi_t
    integer:: size=1, osize=1, rank=0, comm, verbose=0
    integer:: mode, max_tag
    logical:: master
    logical:: ok = .false.
    type(mpi_counter_t):: id, file_size
    type(mpi_counter_t):: end_counter
  contains
    procedure:: init
    procedure:: info
    procedure:: barrier => barrier_
    procedure:: print => print_
    procedure:: assert
    procedure:: abort => abort_mpi
    procedure:: end => end_mpi
  end type
  type(mpi_t), public:: mpi
CONTAINS

!===============================================================================
!> Initialize an object with the most important pieces of information about MPI:
!> mpi%size = the number of MPI ranks
!> mpi%rank = the rank of the current MPI process (0-indexed)
!> mpi%coord = the 3-D position of the rank, if a Cartesian arrangement is used
!===============================================================================
SUBROUTINE init (self, mpi_dims, dims)
  class(mpi_t):: self
  integer, optional, dimension(3):: mpi_dims, dims
  real(8):: offset
  !-----------------------------------------------------------------------------
  if (self%ok) return
  call trace%begin ('mpi_t%init')
  if (first_time) then
    first_time = .false.
    call init_mpi (self)
    call detect_cores
#ifdef MPI
    if (mpi%mode == MPI_THREAD_MULTIPLE) then
      call MPI_Bcast (omp%ncores, 1, MPI_INTEGER, 0, self%comm, mpi_err)
      call omp_timer%get (offset)
      call MPI_Bcast (offset, 1, MPI_REAL8, 0, self%comm, mpi_err)
      call omp_timer%set (offset)
    else
      !$omp critical (mpi_cr)
      call MPI_Bcast (omp%ncores, 1, MPI_INTEGER, 0, self%comm, mpi_err)
      call omp_timer%get (offset)
      call MPI_Bcast (offset, 1, MPI_REAL8, 0, self%comm, mpi_err)
      call omp_timer%set (offset)
      !$omp end critical (mpi_cr)
     end if
#endif
  end if
#ifdef MPI
  mpi%size   = mpi_size
  mpi%rank   = mpi_rank
  mpi%master = mpi_rank==0
#else
  mpi%size   = 1
  mpi%rank   = 0
  mpi%comm   = 0
  mpi%master = (mpi_rank==0)
#endif
  shared%mpi_rank   = mpi%rank
  shared%mpi_size   = mpi%size
  shared%mpi_master = mpi%master
  mpi%ok     = .true.
  !-----------------------------------------------------------------------------
  ! For some reason, this must be done element by element.  Turning it around,
  ! and doing mpi = self (which is syntactically correct), values are screwed
  !-----------------------------------------------------------------------------
  self%size   = mpi%size
  self%rank   = mpi%rank
  self%comm   = mpi%comm
  self%master = mpi%master
  self%ok     = mpi%ok
  !-----------------------------------------------------------------------------
  call mpi%id%init (0)
  call mpi%file_size%init (0)
  call mpi%end_counter%init (0)
  !-----------------------------------------------------------------------------
  ! Number all cores
  !-----------------------------------------------------------------------------
  !$omp parallel
  omp%nthreads_global = mpi%size*omp%nthreads
  omp%thread_global = omp%thread + mpi%rank*omp%nthreads
  !$omp end parallel
  !-----------------------------------------------------------------------------
  ! Make sure output from startup is finished
  !-----------------------------------------------------------------------------
  flush (io_unit%output)
  call self%barrier (delay=0.1)
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE print_ (self)
  class(mpi_t):: self
  if (master) print *, 'mpi%size  =', mpi%size
END SUBROUTINE
!===============================================================================
!> MPI barrier with optional label and delay -- often necessary to arrange serial
!> printout from all ranks, one at a time
!===============================================================================
SUBROUTINE barrier_ (self, label, delay)
  class(mpi_t):: self
  character(len=*), optional:: label
  real, optional:: delay
  integer:: i
  real(8):: wc
  !.............................................................................
  if (.not.self%ok) return
  call trace%begin ('mpi_t%barrier')
  if (present(label)) then
    call barrier_mpi (self, label, delay=delay)
  else
    call barrier_mpi (self, 'mpi%barrier', delay=delay)
  end if
  call trace%end ()
END SUBROUTINE barrier_

!===============================================================================
SUBROUTINE init_mpi (self)
  implicit none
  class(mpi_t):: self
  character(len=mch):: filename
  integer:: mpi_provided
!-------------------------------------------------------------------------------
!  Start up MPI and get number of nodes and our node rank
!-------------------------------------------------------------------------------
  call trace%begin ('mpi_t%init_mpi')
#ifdef MPI
  call MPI_INIT_THREAD (MPI_THREAD_MULTIPLE, mpi_provided, mpi_err)
  self%mode = mpi_provided
  call MPI_Comm_dup (mpi_comm_world, self%comm, mpi_err)
  call MPI_Comm_size (self%comm, mpi_size, mpi_err)
  call MPI_Comm_rank (self%comm, mpi_rank, mpi_err)
  mpi_master = (mpi_rank == 0)
  master = mpi_master
  if (self%verbose > 1) then
    write (filename,'(a,i4.4,a)') 'rank_', mpi_rank, '.dat'
    open (output, file=filename, form='formatted', status='unknown')
  end if
#else
  mpi_rank = 0
  mpi_master = (mpi_rank == 0)
  master = mpi_master
  print*,'MPI is not enabled!'
#endif
  call trace%end ()
END SUBROUTINE init_mpi

!===============================================================================
!> Print information about MPI library support
!===============================================================================
SUBROUTINE info (self)
  class(mpi_t):: self
#ifdef MPI
  integer(8):: tag_ub
  logical:: flag
  !-----------------------------------------------------------------------------
  call MPI_Comm_get_attr (MPI_COMM_WORLD, MPI_TAG_UB, tag_ub, flag, mpi_err)
  if (.not. flag) then
    tag_ub = 2**15
  end if
  self%max_tag = tag_ub
  if (mpi%rank==0) then
    select case (self%mode)
    case(MPI_THREAD_MULTIPLE)
      print 1, ' MPI implementation has support for simultaneous MPI calls by different threads'
    1 format(a,i7)    
    case(MPI_THREAD_SINGLE)
      print 1, ' MPI implementation only has support for MPI calls in one thread at a time'
    case(MPI_THREAD_FUNNELED)
      print 1, ' MPI implementation only has support for MPI calls in master regions'
    case(MPI_THREAD_SERIALIZED)
      print 1, ' MPI implementation only has support for serial MPI calls by any thread'
    case default
      print 1, ' unknown MPI implementation'
    end select
    print '(" MPI tag max =",1p,g10.3," = 2**",i2)', &
      real(tag_ub), nint(log(real(tag_ub))/log(2.))
  end if
#endif
END SUBROUTINE info

!===============================================================================
SUBROUTINE end_mpi (self, label)
  class(mpi_t):: self
  character(len=*), optional:: label
!...............................................................................
  call trace%begin ('mpi_t%end_mpi')
  !call finalize_cart
  if (mpi_rank==0) then
    if (present(label)) then
      print*,'MPI finalized with message: "'//trim(label)//'"'
    end if
  end if
#ifdef MPI
  write (io_unit%mpi,*) self%rank, wallclock(), 'calling MPI_Finalize'
  flush (io_unit%mpi)
  if (mpi%mode == MPI_THREAD_MULTIPLE) then
    call MPI_Finalize (mpi_err)
  else
    !$omp critical (mpi_cr)
    call MPI_Finalize (mpi_err)
    !$omp end critical (mpi_cr)
  end if
  write (io_unit%mpi,*) self%rank, wallclock(), 'returned from MPI_Finalize'
  flush (io_unit%mpi)
#endif
  call exit        ! avoid "stop" which may produce one output line per process
  call trace%end ()
END SUBROUTINE

!===============================================================================
SUBROUTINE abort_mpi (self, reason)
  class(mpi_t):: self
  character(len=*), optional:: reason
  integer:: i
  !-----------------------------------------------------------------------------
  ! Give other ranks a chance to issue messages as well, then abort
  !-----------------------------------------------------------------------------
  flush (io_unit%output)
  call message (stderr)
  call message (stdlog)
#ifdef MPI
  call omp%delay (3e3)
  if (mpi%mode == MPI_THREAD_MULTIPLE) then
    call MPI_Abort (mpi_comm_world, 127, mpi_err)
  else
    !$omp critical (mpi_cr)
    call MPI_Abort (mpi_comm_world, 127, mpi_err)
    !$omp end critical (mpi_cr)
  end if
#endif
  call exit
contains
  !-----------------------------------------------------------------------------
  subroutine message (unit)
  integer:: unit
  integer, parameter:: n=80
  character(len=n):: text
  !.............................................................................
  flush (unit)
  write (unit,*) ''
  write (unit,'(80a)') ('*',i=1,n)
  write (text,'("  ABORT: rank",i5,6x,"thread",i4)') mpi_rank, omp%thread
  write (unit,'(80a)') '**',trim(text),(' ',i=1,n-len(trim(text))-4),'**'
  if (present(reason)) then
    text = ' REASON: '//reason
    write (unit,'(80a)') '**',trim(text),(' ',i=1,n-len(trim(text))-4),'**'
  end if
  write (unit,'(80a)') ('*',i=1,n)
  write (unit,*) ''
  flush (unit)
  end subroutine
END SUBROUTINE abort_mpi

!===============================================================================
SUBROUTINE clean_stop
!...............................................................................
#ifdef MPI
  call omp%delay (3e3)
  if (mpi%mode == MPI_THREAD_MULTIPLE) then
    call MPI_Abort (mpi_comm_world, 127, mpi_err)
  else
    !$omp critical (mpi_cr)
    call MPI_Abort (mpi_comm_world, 127, mpi_err)
    !$omp end critical (mpi_cr)
  end if
#endif
  call exit
END SUBROUTINE clean_stop

!===============================================================================
SUBROUTINE barrier_mpi (self, label, delay)
  class(mpi_t):: self
  real, optional:: delay
  character(len=*) label
!...............................................................................
  if (self%verbose > 2) then
    write (output,*) mpi_rank, 'barrier_mpi: '//label
    flush (output)
  end if
#ifdef MPI
  if (mpi%size > 1) then
    if (present(delay)) then
      call omp%delay (1e3*delay)
    end if
    !$omp barrier
    !$omp single
    call MPI_Barrier (mpi_comm_world, mpi_err)
    !$omp end single
  end if
#endif
END SUBROUTINE barrier_mpi

!===============================================================================
SUBROUTINE assert (self, message, code)
  class(mpi_t):: self
  character(len=*):: message
  integer code
  character (len=64):: codemsg
!..............................................................................
  if (code /= 0) then
    write (codemsg,'(2x,a," error: code =",i10)') message, code
    write (stderr,*) 'mpi_t%assert: rank', mpi_rank, 'mesg:', codemsg
    flush (stderr)
    write (stdlog,*) 'mpi_t%assert: rank', mpi_rank, 'mesg:', codemsg
    flush (stdlog)
  end if
END SUBROUTINE assert

!===============================================================================
!> Attempt to determine the actual number of cores. If /proc/cpuinfo exists
!> one has to count the number of sockets, but observing hom many different
!> "physical id" there are.  The number of cores per socket is in the field
!> "cpu cores".
!===============================================================================
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
  real(8):: wc
  !.............................................................................
  nullify(head)
  n_socket = 0
  n_core = 0
  oid = -1
  inquire (file=file, exist=exist)
  if (exist) then
    wc = wallclock()
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
    wc = wallclock()-wc
    if (wc > 5) print '(a)','############################################################'
    if (mpi%master) &
    print '(a,f6.1,a)',     'reading /proc/cpuinfo:', wc, ' sec'
    if (wc > 5) print '(a)','############################################################'
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

!===============================================================================
!> Create and initialize a counter to value = n. NOTE: This is COLLECTIVE call!
!===============================================================================
SUBROUTINE init4 (self, n)
  class(mpi_counter_t):: self
  integer:: n
  !.............................................................................
  call init8 (self, int(n,kind=8))
END SUBROUTINE init4

!===============================================================================
!> Create and initialize a counter to value = n. NOTE: This is COLLECTIVE call!
!===============================================================================
SUBROUTINE init8 (self, n)
  class(mpi_counter_t):: self
  integer(8):: n
  !.............................................................................
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND):: nbytes=8
  integer:: mpi_err
#endif
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_counter_t%init')
#ifdef MPI
  if (self%window==0) then
  !-----------------------------------------------------------------------------
  ! Do NOT make this depend on MPI_THREAD_MULTIPLE -- this really needs a
  ! critical region
  !-----------------------------------------------------------------------------
    !$omp critical (mpi_cr)
    if (self%window==0) then
      self%i = n
      call MPI_Win_create (self%i, nbytes, 8, MPI_INFO_NULL, MPI_COMM_WORLD, &
                           self%window, mpi_err)
      call mpi%assert ('MPI_Win_create', mpi_err)
    end if
    !$omp end critical (mpi_cr)
  end if
#else
  self%i = n
#endif
  call trace%end ()
END SUBROUTINE init8

!===============================================================================
!> Create and initialize a counter to value = n. NOTE: This is COLLECTIVE call!
!===============================================================================
SUBROUTINE init8d (self, r)
  class(mpi_counter_t):: self
  real(8):: r
  !.............................................................................
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND):: nbytes=8
  integer:: mpi_err
#endif
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_counter_t%init')
#ifdef MPI
  if (self%window==0) then
  !-----------------------------------------------------------------------------
  ! Do NOT make this depend on MPI_THREAD_MULTIPLE -- this really needs a
  ! critical region
  !-----------------------------------------------------------------------------
    !$omp critical (mpi_cr)
    if (self%window==0) then
      self%r = r
      call MPI_Win_create (self%r, nbytes, 8, MPI_INFO_NULL, MPI_COMM_WORLD, &
                           self%window, mpi_err)
      call mpi%assert ('MPI_Win_create', mpi_err)
    end if
    !$omp end critical (mpi_cr)
  end if
#else
  self%r = r
#endif
  call trace%end ()
END SUBROUTINE init8d

!===============================================================================
!> Reset the counter to n
!===============================================================================
SUBROUTINE reset4 (self, n)
  class(mpi_counter_t):: self
  integer:: n, i
  !.............................................................................
  call reset8 (self, int(n,kind=8))
END SUBROUTINE reset4

!===============================================================================
!> Reset the counter to n
!===============================================================================
SUBROUTINE reset8 (self, n)
  class(mpi_counter_t):: self
  integer(8):: n, i
  !.............................................................................
  call trace%begin ('mpi_counter_t%reset')
  i = self%update(0)
  i = self%update(n-i)
  if (io_unit%verbose > 0) then
    write (stderr,*) 'mpi_counter_t%reset: count =', i
  endif
  call trace%end ()
END SUBROUTINE reset8

!===============================================================================
!> Reset the counter to n
!===============================================================================
SUBROUTINE reset8d (self, n)
  class(mpi_counter_t):: self
  real(8):: n, i
  !.............................................................................
  call trace%begin ('mpi_counter_t%reset')
  i = self%update(0.d0)
  i = self%update(n-i)
  if (io_unit%verbose > 0) then
    write (io_unit%log,*) 'mpi_counter_t%reset: count =', i
  endif
  call trace%end ()
END SUBROUTINE reset8d

!===============================================================================
!> Add a value i to the counter on master.  If the value reaches 0, return 0,
!> while resetting the counter to n.  Allow both 4-byte and 8-byte arguments.
!===============================================================================
FUNCTION update4 (self, i) RESULT (j)
  class(mpi_counter_t):: self
  integer:: i, j
  !.............................................................................
  j = update8 (self, int(i,kind=8))
END FUNCTION update4

!===============================================================================
FUNCTION update8 (self, i) RESULT (j)
  class(mpi_counter_t):: self
  integer(8):: i, j, k
  integer, save :: itimer=0
  !.............................................................................
#ifdef MPI
  integer:: master=0, mpi_err
  integer(kind=MPI_ADDRESS_KIND):: offset=0
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_counter_t%update8', 2, itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Do NOT make this depend on MPI_THREAD_MULTIPLE -- this really needs a
  ! critical region and an MPI_LOCK_EXCLUSIVE window lock
  !-----------------------------------------------------------------------------
  !$omp critical (mpi_upd8_cr)
  call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, master, 0, self%window, mpi_err)
  call mpi%assert ('MPI_Win_lock', mpi_err)
  call MPI_Get_accumulate (i, 1, MPI_INTEGER8, &
                           k, 1, MPI_INTEGER8, master, offset, &
                              1, MPI_INTEGER8, MPI_SUM, self%window, mpi_err)
  call mpi%assert ('MPI_Get_accumulate', mpi_err)
  call MPI_Win_unlock (master, self%window, mpi_err)
  call mpi%assert ('MPI_Win_unlock', mpi_err)
  j = k+i
  !$omp end critical (mpi_upd8_cr)
#else
  call trace%begin ('mpi_counter_t%update8', itimer=itimer)
  !$omp atomic capture
  self%i = self%i + i
  j = self%i
  !$omp end atomic
#endif
  if (io_unit%verbose > 0) then
    write (io_unit%log,'(a,2(2x,3i4))') &
      'mpi_counter_t%update: incr, sum, rank, size =', &
      i, j, mpi%rank, mpi%size
  end if
  call trace%end (itimer)
END FUNCTION update8

!===============================================================================
FUNCTION update8d (self, i) RESULT (j)
  class(mpi_counter_t):: self
  real(8):: i, j, k
  integer, save :: itimer=0
  !.............................................................................
#ifdef MPI
  integer:: master=0, mpi_err
  integer(kind=MPI_ADDRESS_KIND):: offset=0
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_counter_t%update8d', 2, itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Do NOT make this depend on MPI_THREAD_MULTIPLE -- this really needs a
  ! critical region
  !-----------------------------------------------------------------------------
  !$omp critical (mpi_upd8d_cr)
  call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, master, 0, self%window, mpi_err)
  call mpi%assert ('MPI_Win_lock', mpi_err)
  call MPI_Get_accumulate (i, 1, MPI_REAL8, &
                           k, 1, MPI_REAL8, master, offset, &
                              1, MPI_REAL8, MPI_SUM, self%window, mpi_err)
  call mpi%assert ('MPI_Get_accumulate', mpi_err)
  call MPI_Win_unlock (master, self%window, mpi_err)
  call mpi%assert ('MPI_Win_unlock', mpi_err)
  j = k+i
  !$omp end critical (mpi_upd8d_cr)
#else
  call trace%begin ('mpi_counter_t%update8d', itimer=itimer)
  !$omp atomic capture
  self%r = self%r + i
  j = self%r
  !$omp end atomic
#endif
  if (io_unit%verbose > 0) then
    write (io_unit%log,'(a,2(2x,3i4))') &
      'mpi_counter_t%update: incr, sum, rank, size =', &
      i, j, mpi%rank, mpi%size
  end if
  call trace%end (itimer)
END FUNCTION update8d

END MODULE mpi_mod

!===============================================================================
!> Enable copying between any 4-byte arrays, without type checking restrictions
!===============================================================================
SUBROUTINE anonymous_copy (n, a, b)
  integer:: n
  real(4):: a(n), b(n)
  b = a
END SUBROUTINE anonymous_copy
