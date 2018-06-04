!> =============================================================================
!> Do not use a dispatcher, but call task_list%execute, which relies on threads
!> handling a ready queue.
!> =============================================================================
MODULE dispatcher_mod
  USE io_mod
  USE trace_mod
  USE task_list_mod
  USE dispatcher1_mod
  USE mpi_mesg_mod
  USE data_io_mod
  !USE task_mod
  !USE patch_mod
  implicit none
  private
  type, public:: dispatcher_t
    procedure(method1), nopass, pointer:: method => null()
  contains
    procedure:: init
    procedure, nopass:: method1
  end type
  integer, save:: verbose=0
  integer, save:: method=0
  type(dispatcher_t), public:: dispatcher
CONTAINS

!> =============================================================================
!> Execute a task list
!> =============================================================================
SUBROUTINE init (self)
  class(dispatcher_t):: self
  integer:: iostat
  namelist /dispatcher_params/ method, verbose
  !-----------------------------------------------------------------------------
  call trace%begin ('dispatcher_t%init')
  rewind (io%input)
  read (io%input, dispatcher_params, iostat=iostat)
  if (io%master) write (io%output, dispatcher_params)
  select case (method)
  case (1)
    self%method => method1
    dispatcher1%verbose = verbose
  case default
    call io%abort ('unknown method in dispatcher_mod')
  end select
  call trace%end
END SUBROUTINE init

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
  mpi_mesg%recv_active = .true.
  mpi_mesg%queue_unpack = .false.
  mpi_mesg%recv_priv = .false.
  mpi_mesg%send_priv = .false.
  dispatcher1%verbose = verbose
  call dispatcher1%execute (task_list)
END SUBROUTINE method1

END MODULE dispatcher_mod
