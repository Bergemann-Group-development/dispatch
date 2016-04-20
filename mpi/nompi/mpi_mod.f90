!> =============================================================================
!> Basic MPI calls; initialization, end, barriers, asserts, wallclock for use
!> when MPI is *NOT* activated.
!>
!> =============================================================================
MODULE mpi_mod
  implicit none
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
  mpi%size   = 1
  mpi%rank   = 0
  mpi%comm   = 0
  mpi%master = (mpi_rank==0)
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
    call barrier_mpi (mp, 'nompi%barrier')
  end if
END SUBROUTINE barrier_

!> =============================================================================
SUBROUTINE init_mpi (mp)
  implicit none
  class(mpi_t):: mp
  character(len=mch):: filename
  integer:: mpi_provided
!-------------------------------------------------------------------------------
!  Start up MPI and get number of nodes and our node rank
!-------------------------------------------------------------------------------
  print*,'MPI is not enabled! This subroutine does nothing!'
END SUBROUTINE init_mpi

!> =============================================================================
SUBROUTINE end_mpi (mp, label)
  class(mpi_t):: mp
  character(len=*), optional:: label
!..............................................................................
  !call finalize_cart
  if (mpi_rank==0) then
    if (present(label)) then
      print*,'NOMPI finalized with message: "'//trim(label)//'"'
    else
      print*,'NOMPI finalized'
    end if
  end if
  call exit        ! avoid "stop" which may produce one output line per process
END SUBROUTINE

!> =============================================================================
SUBROUTINE abort_mpi (mp, reason)
  class(mpi_t):: mp
  character(len=*):: reason
!...............................................................................
  print*,'nompi_abort: '//trim(reason)
  if (mp%verbose > 2) then
    write(output,'(3(a,2x,i6))') 'nompi_abort: rank =',mpi_rank,'thread',trim(reason)
    call flush (output)
  end if
  call exit
END SUBROUTINE abort_mpi

!> =============================================================================
SUBROUTINE clean_stop
!...............................................................................
  call exit
END SUBROUTINE clean_stop

!> =============================================================================
SUBROUTINE barrier_mpi (mp, label)
  class(mpi_t):: mp
  character(len=*) label
!...............................................................................
  if (mpi_trace) print *,  mpi_rank, 'barrier_nompi: '//label
  if (mp%verbose > 2) then
    write (output,*) mpi_rank, 'barrier_nompi: '//label
    call flush (output)
  end if
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
