!> =============================================================================
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
!> =============================================================================
MODULE mpi_mod
  implicit none
  include 'mpif.h'              ! some systems require this instead of 'use mpi'
  integer:: mpi_err             ! MPI error status
  logical:: mpi_master          ! true on MPI master (rank 0)
  logical:: mpi_trace=.false.   ! debug trace option
  logical:: master=.true.       ! true on OMP master thread on MPI master process
  integer:: mpi_comm            ! short hand for mpi_comm_world
  integer:: mpi_size            ! number of MPI processes
  integer:: mpi_rank            ! process number (zero-based)
  integer:: output=2
  integer, parameter:: mch=120
  logical:: first_time = .true.
  !-----------------------------------------------------------------------------
  ! Object that returns the most important MPI parameters
  !-----------------------------------------------------------------------------
  type mpi_t
    integer:: size, rank, comm, verbose=0
    logical:: master
    logical:: ok = .false.
  contains
    procedure:: init
    procedure:: barrier => barrier_
    procedure:: print => print_
    procedure:: assert
    procedure:: abort => abort_mpi
    procedure:: end => end_mpi
  end type
  type(mpi_t):: mpi
PRIVATE
PUBLIC mpi_t, mpi
CONTAINS

!> =============================================================================
!> Initialize an object with the most important pieces of information about MPI:
!> mpi%size = the number of MPI ranks
!> mpi%rank = the rank of the current MPI process (0-indexed)
!> mpi%coord = the 3-D position of the rank, if a Cartesian arrangement is used
!> =============================================================================
SUBROUTINE init (mp, mpi_dims, dims)
  class(mpi_t):: mp
  integer, optional, dimension(3):: mpi_dims, dims
  integer:: rank
  !.............................................................................
  if (first_time) then
    first_time = .false.
    call init_mpi (mp)
  end if
  mpi%size   = mpi_size
  mpi%rank   = mpi_rank
  mpi%comm   = mpi_comm_world
  mpi%master = mpi_rank==0
  mpi%ok     = .true.
  !-----------------------------------------------------------------------------
  ! For some reason, this must be done element by element.  Turning it around,
  ! and doing mpi = mp (which is syntactically correct), values are screwed
  !-----------------------------------------------------------------------------
  mp%size   = mpi%size
  mp%rank   = mpi%rank
  mp%comm   = mpi%comm
  mp%master = mpi%master
  mp%ok     = mpi%ok
END SUBROUTINE init

!> =============================================================================
!> =============================================================================
SUBROUTINE print_ (mp)
  class(mpi_t):: mp
  if (master) print *, 'mpi%size  =', mpi%size
END SUBROUTINE
!> =============================================================================
!> MPI barrier with optional label and delay -- often necessary to arrange serial
!> printout from all ranks, one at a time
!> =============================================================================
SUBROUTINE barrier_ (mp, label, delay)
  class(mpi_t):: mp
  character(len=*), optional:: label
  integer, optional:: delay
  integer:: i
  real:: a = 1.0
  !.............................................................................
  if (present(delay)) then
    !call sleep(delay)
    do i=1,delay*1000000
      a = a + sin(a)
    end do
  end if
  if (present(label)) then
    call barrier_mpi (mp, label)
  else
    call barrier_mpi (mp, 'mpi%barrier')
  end if
END SUBROUTINE barrier_

!> =============================================================================
SUBROUTINE init_mpi (mp)
  implicit none
  class(mpi_t):: mp
  character(len=mch):: id='mpi.f90 $Id$'
  character(len=mch):: filename
  integer:: mpi_provided
!-------------------------------------------------------------------------------
!  Start up MPI and get number of nodes and our node rank
!-------------------------------------------------------------------------------
  call MPI_INIT_THREAD (MPI_THREAD_MULTIPLE, mpi_provided, mpi_err)
  mpi_comm = mpi_comm_world
  call MPI_Comm_size (mpi_comm_world, mpi_size, mpi_err)
  call MPI_Comm_rank (mpi_comm_world, mpi_rank, mpi_err)
  mpi_master = (mpi_rank == 0)
  master = mpi_master
  if (master) then
    print 1, '--------------------------------------------------------------------------------'
    print 1, ' '//id
    print 1, ' MPI initialized, mpi%size =', mpi_size
  1 format(a,i7)
    select case (mpi_provided)
    case(MPI_THREAD_MULTIPLE)
      print 1, ' MPI implementation provides support for simultaneous MPI calls by different threads'
    case(MPI_THREAD_SINGLE)
      print 1, ' MPI implementation provides only support for MPI calls outside openmp regions'
    case(MPI_THREAD_FUNNELED)
      print 1, ' MPI implementation provides only support for MPI calls in master regions'
    case(MPI_THREAD_SERIALIZED)
      print 1, ' MPI implementation provides only support for serial MPI calls by any thread'
    case default
      print 1, ' unknown MPI implementation'
    end select
    print 1, '--------------------------------------------------------------------------------'
  end if
  if (mp%verbose > 1) then
    write (filename,'(a,i4.4,a)') 'rank_', mpi_rank, '.dat'
    open (output, file=filename, form='formatted', status='unknown')
  end if
END SUBROUTINE init_mpi

!> =============================================================================
SUBROUTINE end_mpi (mp, label)
  class(mpi_t):: mp
  character(len=*), optional:: label
!..............................................................................
  !call finalize_cart
  if (mpi_rank==0) then
    if (present(label)) then
      print*,'MPI finalized with message: "'//trim(label)//'"'
    end if
  end if
  call MPI_Finalize (mpi_err)
  call exit        ! avoid "stop" which may produce one output line per process
END SUBROUTINE

!> =============================================================================
SUBROUTINE abort_mpi (mp, reason)
  class(mpi_t):: mp
  character(len=*):: reason
!...............................................................................
  print*,'mpi_abort: '//trim(reason)
  if (mp%verbose > 2) then
    write(output,'(3(a,2x,i6))') 'mpi_abort: rank =',mpi_rank,'thread',trim(reason)
    call flush (output)
  end if
  call MPI_Abort (mpi_comm_world, 127, mpi_err)
  call exit
END SUBROUTINE abort_mpi

!> =============================================================================
SUBROUTINE clean_stop
!...............................................................................
  call MPI_Abort (mpi_comm_world, 127, mpi_err)
  call exit
END SUBROUTINE clean_stop

!> =============================================================================
SUBROUTINE barrier_mpi (mp, label)
  class(mpi_t):: mp
  character(len=*) label
!...............................................................................
  if (mpi_trace) print *,  mpi_rank, 'barrier_mpi: '//label
  if (mp%verbose > 2) then
    write (output,*) mpi_rank, 'barrier_mpi: '//label
    call flush (output)
  end if
  call MPI_Barrier (mpi_comm_world, mpi_err)
END SUBROUTINE barrier_mpi

!> =============================================================================
SUBROUTINE assert (mp, message, code)
  class(mpi_t):: mp
  character(len=*):: message
  integer code
  character (len=16):: codemsg
!..............................................................................
  if (code/=0) then
    write (codemsg,'("code =",i10)') code
    write(*,*)      mpi_rank,' error: '//trim(message)
    if (mp%verbose > 1) then
      write(output,*) mpi_rank,' error: '//trim(message)
      call flush (output)
    end if
  end if
END SUBROUTINE assert

END MODULE mpi_mod
