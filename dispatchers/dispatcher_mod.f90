!===============================================================================
!> Do not use a dispatcher, but call task_list%execute, which relies on threads
!> handling a ready queue.
!===============================================================================
MODULE dispatcher_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE task_list_mod
  USE dispatcher0_mod
  USE dispatcher1_mod
  USE dispatcher2_mod
  USE dispatcher3_mod
  USE dispatcher4_mod
  USE dispatcher5_mod
  USE dispatcher6_mod
  USE mpi_mesg_mod
  USE data_io_mod
  USE experiment_mod
  USE list_mod
  USE task_mesg_mod
  implicit none
  private
  type, public:: dispatcher_t
    procedure(method0), nopass, pointer:: method => null()
  contains
    procedure:: init
    procedure:: execute
    procedure, nopass:: method0
    procedure, nopass:: method1
    procedure, nopass:: method2
    procedure, nopass:: method3
    procedure, nopass:: method4
    procedure, nopass:: method5
  end type
  integer, save:: verbose=0
  integer, save:: method=0
  logical, save:: test=.false.
  real(8):: test_seconds=30.
  type(dispatcher_t), public:: dispatcher
CONTAINS

!===============================================================================
!> Choose and initialize the dispatcher method.  For full flexibility, default
!> parameters for mpi_mesg_mod are first set, and may then be overwritten in the
!> mpi_mesg%init call.
!===============================================================================
SUBROUTINE init (self)
  class(dispatcher_t):: self
  integer:: iostat
  namelist /dispatcher_params/ method, verbose, test, test_seconds
  !-----------------------------------------------------------------------------
  call trace%begin ('dispatcher_t%init')
  rewind (io%input)
  read (io%input, dispatcher_params, iostat=iostat)
  write (io%output, dispatcher_params)
  !-----------------------------------------------------------------------------
  select case (method)
  case (0)
    self%method => method0
    call dispatcher0%init
  case (1)
    self%method => method1
  case (2)
    self%method => method2
  case (3)
    self%method => method3
  case (4)
    self%method => method4
  case (5)
    self%method => method5
  case (6)
    self%method => method6
  case default
    call io%abort ('unknown method in dispatcher_mod')
  end select
  !-----------------------------------------------------------------------------
  ! Pass the method choice on to task_method_t, in case it has an impact
  !-----------------------------------------------------------------------------
  task_mesg%method = method
  call mpi_mesg%init
  call trace%end()
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE execute (self, task_list)
  class(dispatcher_t):: self
  type(task_list_t), pointer:: task_list
  class(list_t), pointer:: list
  class(link_t), pointer:: link
  !-----------------------------------------------------------------------------
  call trace%begin ('dispatcher_t%execute')
  task_list%method = method
  task_list%dispatcher = .true.
  if (test) call set_io
  !-----------------------------------------------------------------------------
  ! Reset task_list status to re-establish nbor relations, checking consistency
  !-----------------------------------------------------------------------------
  call task_list%reset_status (check=.true.)
  call task_list%init_task_list_pointers (task_list)
  call task_list%info
  !-----------------------------------------------------------------------------
  ! Call the selected dispatcher method
  !-----------------------------------------------------------------------------
  call self%method (task_list)
  call trace%end()
END SUBROUTINE execute

!===============================================================================
!> I/O params for tests
!===============================================================================
SUBROUTINE set_io
  io%do_output = .false.
  io%print_time = 50.
  io%end_time = 1e10
  io%job_seconds = test_seconds
END SUBROUTINE set_io

!===============================================================================
!> Execute a task list, using the selected dispatcher method
!===============================================================================
SUBROUTINE method0 (task_list)
  type(task_list_t), pointer:: task_list
  !-----------------------------------------------------------------------------
  task_list%method = method
  if (test) call set_io
  call dispatcher0%execute (task_list, test)
END SUBROUTINE method0

!===============================================================================
SUBROUTINE method1 (task_list)
  type(task_list_t), pointer:: task_list
  if (test) call set_io
  call dispatcher1%execute (task_list, test)
END SUBROUTINE method1

!===============================================================================
SUBROUTINE method2 (task_list)
  type(task_list_t), pointer:: task_list
  if (test) call set_io
  call dispatcher2%execute (task_list, test)
END SUBROUTINE method2

!===============================================================================
SUBROUTINE method3 (task_list)
  type(task_list_t), pointer:: task_list
  if (test) call set_io
  call dispatcher3%execute (task_list, test)
END SUBROUTINE method3

!===============================================================================
SUBROUTINE method4 (task_list)
  type(task_list_t), pointer:: task_list
  if (test) call set_io
  call dispatcher4%execute (task_list, test)
END SUBROUTINE method4

!===============================================================================
SUBROUTINE method5 (task_list)
  type(task_list_t), pointer:: task_list
  if (test) call set_io
  call dispatcher5%execute (task_list, test)
END SUBROUTINE method5

!===============================================================================
SUBROUTINE method6 (task_list)
  type(task_list_t), pointer:: task_list
  if (test) call set_io
  call dispatcher6%execute (task_list, test)
END SUBROUTINE method6

END MODULE dispatcher_mod
