!> =============================================================================
!> $Id$
!> Basic MPI calls; initialization, end, barriers, asserts, wallclock, ...
!> =============================================================================
MODULE mpi_basic
  USE omp_mod                   ! OMP stuff
 !USE mpi                       ! Most systems accept this, but ..
  implicit none
  include 'mpif.h'              ! .. some systems may require this instead
  integer:: mpi_err             ! MPI error status
  logical:: mpi_master          ! true on MPI master (rank 0)
  logical:: mpi_trace=.false.   ! debug trace option
  logical:: master=.true.       ! true on OMP master thread on MPI master process
  integer:: mpi_comm            ! short hand for mpi_comm_world
  integer:: mpi_size            ! number of MPI processes
  integer:: mpi_rank            ! process number (zero-based)
  integer, parameter:: mch=120  ! max character string size
  public mpi_comm_world         ! main MPI communicator
  public MPI_INTEGER            ! MPI magic number for integers
  public MPI_INT                ! ditto
  public MPI_REAL               ! for (4-byte) real's
  public MPI_STATUS_SIZE
  integer:: input=1, output=2, data_unit=3, log_unit=4
  logical:: do_debug = .false.
  logical:: do_trace = .true.
  integer:: verbose = 1
  !-----------------------------------------------------------------------------
  ! Object that returns the most important MPI parameters
  !-----------------------------------------------------------------------------
  type mpi_t
    integer:: size, rank
  contains
    procedure:: init => init_mpi_info
  end type
  type(mpi_t):: info
PUBLIC
CONTAINS

!> =============================================================================
!> Initialize an object with the most important pieces of information about MPI:
!> mpi%size = the number of MPI ranks
!> mpi%rank = the rank of the current MPI process (0-indexed)
!> mpi%coord = the 3-D position of the rank, if a Cartesian arrangement is used
!> =============================================================================
SUBROUTINE init_mpi_info (mpi)
  class(mpi_t):: mpi
  !.............................................................................
  mpi%size = mpi_size
  mpi%rank = mpi_rank
END SUBROUTINE init_mpi_info

!> =============================================================================
SUBROUTINE end_mpi
!..............................................................................
  !call finalize_cart
  call MPI_Finalize (mpi_err)
  call exit        ! avoid "stop" which may produce one output line per process
END SUBROUTINE

!> =============================================================================
SUBROUTINE abort_mpi (reason)
  character(len=*):: reason
!...............................................................................
  print*,'mpi_abort: '//trim(reason)
  write(output,'(3(a,2x,i6))') 'mpi_abort: rank =',mpi_rank,'thread',omp_mythread(), trim(reason)
  call flush (output)
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
SUBROUTINE barrier_mpi (label)
  character(len=*) label
!...............................................................................
  if (mpi_trace) print *,  mpi_rank, 'barrier_mpi: '//label
  write (output,*) mpi_rank, 'barrier_mpi: '//label
  call flush (output)
  call MPI_Barrier (mpi_comm_world, mpi_err)
END SUBROUTINE barrier_mpi

!> =============================================================================
SUBROUTINE assert (message, flag)
  logical:: flag
  character(len=*):: message
!..............................................................................
  if (.not.flag) then
    write(*,*)      mpi_rank,' error: '//trim(message)
    write(output,*) mpi_rank,' error: '//trim(message)
    call flush (output)
  end if
END SUBROUTINE assert

!> =============================================================================
SUBROUTINE assert_mpi (message, code)
  character(len=*):: message
  integer code
  character (len=16):: codemsg
!..............................................................................
  if (code/=0) &
    print*,mpi_rank,' error: '//trim(message)//', code =', code
  write (codemsg,'("code =",i10)') code
  call assert (trim(message)//codemsg, code==0)
END SUBROUTINE assert_mpi

END MODULE mpi_basic

!> =============================================================================
SUBROUTINE init_mpi
  USE mpi_basic
  implicit none
  character(len=mch):: id='mpi.f90 $Id$'
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
  if (master) print*,id
  call init_omp
  if (master) then
    print 1, '***********************************************************************************'
    select case (mpi_provided)
    case(MPI_THREAD_MULTIPLE)
      print 1, 'MPI implementation provides support for simultaneous MPI calls by different threads'
      1 format(1x,a,i7)
    case(MPI_THREAD_SINGLE)
      print 1, 'MPI implementation provides only support for MPI calls outside openmp regions'
    case(MPI_THREAD_FUNNELED)
      print 1, 'MPI implementation provides only support for MPI calls in master regions'
    case(MPI_THREAD_SERIALIZED)
      print 1, 'MPI implementation provides only support for serial MPI calls by any thread'
    case default
      print 1, 'unknown MPI implementation'
    end select
    print 1, '***********************************************************************************'
  end if
  !
  if (master) then
    print 1,'number of MPI processes:', mpi_size
    print 1,'OMP threads per MPI process:', omp_nthreads
  end if
END SUBROUTINE init_mpi
