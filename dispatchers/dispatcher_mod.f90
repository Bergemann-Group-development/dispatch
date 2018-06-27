!> =============================================================================
!> Do not use a dispatcher, but call task_list%execute, which relies on threads
!> handling a ready queue.
!> =============================================================================
MODULE dispatcher_mod
  USE io_mod
  USE trace_mod
  USE task_list_mod
  USE dispatcher0_mod
  USE dispatcher1_mod
  USE dispatcher2_mod
  USE mpi_mesg_mod
  USE data_io_mod
  implicit none
  private
  type, public:: dispatcher_t
    procedure(method0), nopass, pointer:: method => null()
  contains
    procedure:: init
    procedure, nopass:: method0
    procedure, nopass:: method1
    procedure, nopass:: method2
  end type
  integer, save:: verbose=0
  integer, save:: method=0
  logical, save:: test=.false.
  real(8):: test_seconds=30.
  type(dispatcher_t), public:: dispatcher
CONTAINS

!> =============================================================================
!> Choose and initialize the dispatcher method.  For full flexibility, default
!> parameters for mpi_mesg_mod are first set, and may then be overwritten in the
!> mpi_mesg%init call.
!> =============================================================================
SUBROUTINE init (self)
  class(dispatcher_t):: self
  integer:: iostat
  namelist /dispatcher_params/ method, verbose, test, test_seconds
  !-----------------------------------------------------------------------------
  call trace%begin ('dispatcher_t%init')
  rewind (io%input)
  read (io%input, dispatcher_params, iostat=iostat)
  if (io%master) write (io%output, dispatcher_params)
  select case (method)
  case (0)
    self%method => method0
    call dispatcher0%init
    !---------------------------------------------------------------------------
    ! By default this method uses the most optimal send and receive procedures,
    ! which means each thread takes care of its own sent_list, and each thread
    ! checks a dedicated set of tasks for incoming messages.
    !---------------------------------------------------------------------------
    mpi_mesg%recv_active  = .true.
    mpi_mesg%queue_unpack = .true.
    mpi_mesg%recv_priv    = .true.
    mpi_mesg%send_priv    = .true.
  case (1)
    self%method => method1
    dispatcher1%verbose = verbose
    !---------------------------------------------------------------------------
    ! This method uses send_priv=.true. to make each thread takes care sending.
    ! It uses the check_active() procedure in task_mesg_mod to actively check
    ! for incoming messages.  Since this is done by only the master thread,
    ! there is no point in making an unpack queue.
    !---------------------------------------------------------------------------
    mpi_mesg%send_priv    = .true.
    mpi_mesg%recv_active  = .true.
    mpi_mesg%recv_priv    = .false.
    mpi_mesg%queue_unpack = .false.
  case (2)
    self%method => method2
    !---------------------------------------------------------------------------
    ! This method uses send_priv=.true. to make each thread takes care sending.
    ! It uses the check_active() procedure in task_mesg_mod to actively check
    ! for incoming messages.  Since this is done by only the master thread,
    ! there is no point in making an unpack queue.
    !---------------------------------------------------------------------------
    mpi_mesg%send_priv    = .true.
    mpi_mesg%recv_active  = .true.
    mpi_mesg%recv_priv    = .false.
    mpi_mesg%queue_unpack = .false.
  case default
    call io%abort ('unknown method in dispatcher_mod')
  end select
  call mpi_mesg%init
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> I/O params for tests
!> =============================================================================
SUBROUTINE set_io
  io%do_output = .false.
  io%print_time = 50.
  io%end_time = 1e10
  io%job_seconds = test_seconds
END SUBROUTINE set_io

!> =============================================================================
!> Execute a task list
!> =============================================================================
SUBROUTINE method0 (task_list)
  type(task_list_t), pointer:: task_list
  !-----------------------------------------------------------------------------
  task_list%method = method
  io%ntask = task_list%na
  if (test) call set_io
  call dispatcher0%execute (task_list, test)
END SUBROUTINE method0

!> =============================================================================
!> Execute a task list.  Since only the master thread is call check_mpi(), we
!> need to make sure that it handles all send/recv, hence send_priv=.false. and
!> recv_priv=.false., and we can use the simplest receive mechanism (recv_active
!> =.true.), with no buffering needed (queue_unpack=.false.).  These values are
!> imposed after the input namelists are read, to make sure the choices are
!> correct.
!> =============================================================================
SUBROUTINE method1 (task_list)
  type(task_list_t), pointer:: task_list
  !-----------------------------------------------------------------------------
  task_list%method = method
  io%ntask = task_list%na
  dispatcher1%verbose = verbose
  if (test) call set_io
  call dispatcher1%execute (task_list, test)
END SUBROUTINE method1

!> =============================================================================
!> Execute a task list
!> =============================================================================
SUBROUTINE method2 (task_list)
  type(task_list_t), pointer:: task_list
  !-----------------------------------------------------------------------------
  task_list%method = method
  io%ntask = task_list%na
  if (test) call set_io
  call dispatcher2%execute (task_list, test)
END SUBROUTINE method2

END MODULE dispatcher_mod
