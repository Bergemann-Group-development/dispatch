!> =============================================================================
!> Do not use a dispatcher, but call task_list%execute, which relies on threads
!> handling a ready queue.
!> =============================================================================
MODULE dispatcher_mod
  USE io_mod
  USE task_list_mod
  USE dispatcher_method1_mod
  USE dispatcher_method2_mod
  USE mpi_mesg_mod
  implicit none
  private
  type, public:: dispatcher_t
    procedure(none), nopass, pointer:: method => null()
  contains
    procedure:: init
    procedure, nopass:: none
    procedure, nopass:: method1
    procedure, nopass:: method2
  end type
  integer:: method=0
  type(dispatcher_t), public:: dispatcher
CONTAINS

!> =============================================================================
!> Execute a task list
!> =============================================================================
SUBROUTINE init (self)
  class(dispatcher_t):: self
  integer:: iostat
  integer, save:: verbose=0
  namelist /dispatcher_params/ method, verbose
  !-----------------------------------------------------------------------------
  rewind (io%input)
  read (io%input, dispatcher_params, iostat=iostat)
  write (io%output, dispatcher_params)
  select case (method)
  case (0)
    self%method => none
  case (1)
    self%method => method1
    dispatcher_method1%verbose = verbose
  case (2)
    self%method => method2
  case default
    call io%abort ('unknown method in dispatcher_mod')
  end select
END SUBROUTINE init

!> =============================================================================
!> Execute a task list
!> =============================================================================
SUBROUTINE none (task_list)
  type(task_list_t), pointer:: task_list
  !-----------------------------------------------------------------------------
  call task_list%execute
END SUBROUTINE none

!> =============================================================================
!> Execute a task list
!> =============================================================================
SUBROUTINE method1 (task_list)
  type(task_list_t), pointer:: task_list
  !-----------------------------------------------------------------------------
  task_list%recv_active = .true.
  mpi_mesg%recv_priv = .false.
  call dispatcher_method1%execute (task_list)
END SUBROUTINE method1

!> =============================================================================
!> Execute a task list
!> =============================================================================
SUBROUTINE method2 (task_list)
  type(task_list_t), pointer:: task_list
  !-----------------------------------------------------------------------------
  task_list%recv_active = .true.
  mpi_mesg%recv_priv = .false.
  call dispatcher_method2%execute (task_list)
END SUBROUTINE method2

END MODULE dispatcher_mod
