!> =============================================================================
!> $Id: fa7d9b42d0f7e0aec02fc33bcaebb2ccbdbf57ff $
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
  private:: print_id
  integer:: input=1, output=2, data_unit=3, log_unit=4
PUBLIC
CONTAINS

!> =============================================================================
SUBROUTINE init_mpi
  character(len=mch):: id='mpi.f90 $Id: fa7d9b42d0f7e0aec02fc33bcaebb2ccbdbf57ff $'
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
  call print_id (id)
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
FUNCTION wallclock() result (delta)
  real(kind=8):: time, delta
  real(kind=8), save:: offset=-1.
!...............................................................................
#ifdef _OPENMP
  real(kind=8):: OMP_Get_Wtime, OMP_Get_Wtick
  logical:: OMP_in_parallel
  if (OMP_in_parallel()) then
    time = OMP_Get_Wtime()
  else
    !$omp parallel
    time = OMP_Get_Wtime()
    !$omp end parallel
  end if
  if (offset == -1.) then
    !print *,'tick:', OMP_Get_Wtick()
#else
  time = MPI_Wtime()
  if (offset == -1.) then
    !print *,'tick:', MPI_Wtick()
#endif
    offset=time
  end if
  delta = time-offset
END FUNCTION wallclock

!> =============================================================================
SUBROUTINE print_id (id)
  character(len=mch) id
  character(len=mch), save:: hl= &
    '--------------------------------------------------------------------------------'
!..............................................................................
  if (id .ne. ' ' .or. mpi_trace) then
    if (master) then
      print'(a)',trim(hl)
      print'(1x,a,f12.2)', trim(id), wallclock()
    end if
    if (.not.mpi_trace) id = ' '
  end if
END SUBROUTINE print_id

!> =============================================================================
SUBROUTINE trace_mpi (message)
  character(len=*):: message
!..............................................................................
  if (mpi_trace) print*,mpi_rank,' trace: '//trim(message)
END SUBROUTINE trace_mpi

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
