!===============================================================================
!> $Id$
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
  implicit none
#ifdef MPI
  include 'mpif.h'              ! some systems require this instead of 'use mpi'
#endif
  integer:: mpi_err             ! MPI error status
  logical:: mpi_master          ! true on MPI master (rank 0)
  logical:: mpi_trace=.false.   ! debug trace option
  logical:: master=.true.       ! true on OMP master thread on MPI master process
  integer:: mpi_size            ! number of MPI processes
  integer:: mpi_rank            ! process number (zero-based)
  integer:: output=2
  logical:: first_time = .true.
  !-----------------------------------------------------------------------------
  ! Object that returns the most important MPI parameters
  !-----------------------------------------------------------------------------
  type mpi_t
    integer:: size=1, rank=0, comm, verbose=0
    integer:: mode
    logical:: master
    logical:: ok = .false.
  contains
    procedure:: init
    procedure:: barrier => barrier_
    procedure:: print => print_
    procedure:: assert
    procedure:: abort => abort_mpi
    procedure:: end => end_mpi
    procedure, nopass:: delay
  end type
  type(mpi_t):: mpi
PRIVATE
PUBLIC mpi_t, mpi
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
  character(len=120):: ids = &
  '$Id$ mpi/mpi_mod.f90'
  !-----------------------------------------------------------------------------
  if (self%ok) return
  if (first_time) then
    first_time = .false.
    call init_mpi (self)
    call omp%init
    call detect_cores
#ifdef MPI
    call MPI_Bcast (omp%ncores, 1, MPI_INTEGER, 0, self%comm, mpi_err)
    call omp_timer%get (offset)
    call MPI_Bcast (offset, 1, MPI_REAL8, 0, self%comm, mpi_err)
    call omp_timer%set (offset)
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
  mpi%ok     = .true.
  if (mpi%master .and. first_time) &
    print'(1x,a)', trim(ids)
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
  ! Make sure output from startup is finished
  !-----------------------------------------------------------------------------
  flush (io_unit%output)
  call self%barrier (delay=0.1)
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
  if (present(delay)) then
    wc = wallclock()
    do while (wallclock()-wc < delay)
    end do
  end if
  if (present(label)) then
    call barrier_mpi (self, label)
  else
    call barrier_mpi (self, 'mpi%barrier')
  end if
END SUBROUTINE barrier_

!===============================================================================
SUBROUTINE init_mpi (self)
  implicit none
  class(mpi_t):: self
  character(len=mch):: id='mpi.f90 $Id$'
  character(len=mch):: filename
  integer:: mpi_provided
!-------------------------------------------------------------------------------
!  Start up MPI and get number of nodes and our node rank
!-------------------------------------------------------------------------------
#ifdef MPI
#ifdef __GFORTRAN__
  call MPI_INIT_THREAD (MPI_THREAD_SINGLE, mpi_provided, mpi_err)
#else
  call MPI_INIT_THREAD (MPI_THREAD_MULTIPLE, mpi_provided, mpi_err)
#endif
  self%mode = mpi_provided
  call MPI_Comm_dup (mpi_comm_world, self%comm, mpi_err)
  call MPI_Comm_size (self%comm, mpi_size, mpi_err)
  call MPI_Comm_rank (self%comm, mpi_rank, mpi_err)
  mpi_master = (mpi_rank == 0)
  master = mpi_master
  if (master) then
    print 1, '--------------------------------------------------------------------------------'
    !!!print 1, ' '//id
    print 1, ' MPI initialized, mpi%size =', mpi_size
  1 format(a,i7)
    select case (mpi_provided)
    case(MPI_THREAD_MULTIPLE)
      print 1, ' MPI implementation provides support for simultaneous MPI calls by different threads'
    case(MPI_THREAD_SINGLE)
      print 1, ' MPI implementation only provides support for MPI calls in one thread at a time'
    case(MPI_THREAD_FUNNELED)
      print 1, ' MPI implementation only provides support for MPI calls in master regions'
    case(MPI_THREAD_SERIALIZED)
      print 1, ' MPI implementation only provides support for serial MPI calls by any thread'
    case default
      print 1, ' unknown MPI implementation'
    end select
    print 1, '--------------------------------------------------------------------------------'
  end if
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
END SUBROUTINE init_mpi

!===============================================================================
SUBROUTINE end_mpi (self, label)
  class(mpi_t):: self
  character(len=*), optional:: label
!..............................................................................
  !call finalize_cart
  if (mpi_rank==0) then
    if (present(label)) then
      print*,'MPI finalized with message: "'//trim(label)//'"'
    end if
  end if
#ifdef MPI
  write (io_unit%mpi,*) 'calling MPI_Finalize'
  flush (io_unit%mpi)
  call MPI_Finalize (mpi_err)
  write (io_unit%mpi,*) 'returned from MPI_Finalize'
  flush (io_unit%mpi)
#endif
  !call exit        ! avoid "stop" which may produce one output line per process
END SUBROUTINE

!===============================================================================
SUBROUTINE abort_mpi (self, reason)
  class(mpi_t):: self
  character(len=*), optional:: reason
!...............................................................................
#ifdef MPI
  !-----------------------------------------------------------------------------
  ! Give other ranks a chance to issue messages as well, then abort
  !-----------------------------------------------------------------------------
  flush (io_unit%output)
  call message (stderr)
  call message (io_unit%log)
  call delay (3e3)
  call MPI_Abort (mpi_comm_world, 127, mpi_err)
#endif
  call exit
contains
  subroutine message (unit)
  integer:: unit
  flush (unit)
  if (present(reason)) then
    write (unit,'("ABORT:  rank",i5,4x,"thread",i4,4x,"reason: ",a)') &
      mpi_rank, omp%thread, trim(reason)
  else
    write (unit,'("ABORT:  rank",i5,4x,"thread",i4,4x,"generic")') &
      mpi_rank, omp%thread
  end if
  flush (unit)
  end subroutine
END SUBROUTINE abort_mpi

!===============================================================================
SUBROUTINE clean_stop
!...............................................................................
#ifdef MPI
  call MPI_Abort (mpi_comm_world, 127, mpi_err)
#endif
  call exit
END SUBROUTINE clean_stop

!===============================================================================
SUBROUTINE barrier_mpi (self, label)
  class(mpi_t):: self
  character(len=*) label
!...............................................................................
  if (mpi_trace) print *,  mpi_rank, 'barrier_mpi: '//label
  if (self%verbose > 2) then
    write (output,*) mpi_rank, 'barrier_mpi: '//label
    flush (output)
  end if
#ifdef MPI
  call MPI_Barrier (mpi_comm_world, mpi_err)
#endif
END SUBROUTINE barrier_mpi

!===============================================================================
SUBROUTINE assert (self, message, code)
  class(mpi_t):: self
  character(len=*):: message
  integer code
  character (len=16):: codemsg
!..............................................................................
  if (code/=0) then
    write (codemsg,'("code =",i10)') code
    write(*,*)      mpi_rank,' error: '//trim(message)
    if (self%verbose > 1) then
      write(output,*) mpi_rank,' error: '//trim(message)
      flush (output)
    end if
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
!> Delay for about 1 ms, to encourage incoming messages to complete
!===============================================================================
SUBROUTINE delay (ms)
  real:: ms
  real(8):: wc
  !-----------------------------------------------------------------------------
  wc = wallclock()
  do while ((wallclock()-wc) < 1e-3*ms)
  end do
END SUBROUTINE delay

END MODULE mpi_mod

!===============================================================================
!> Enable copying between any 4-byte arrays, without type checking restrictions
!===============================================================================
SUBROUTINE anonymous_copy (n, a, b)
  integer:: n
  real(4):: a(n), b(n)
  b = a
END SUBROUTINE anonymous_copy
