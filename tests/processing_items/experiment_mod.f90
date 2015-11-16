!> =============================================================================
!> $Id$
!> Prototype experiment, showing the steps and components needed:
!> 1) a "create" procedure, which creates tasks and adds them to a task_list
!> 2) an "init" proc, which will be executed by the task_list%init
!> 3) an "update" procedure, which will be executed by task_queue%update
!> 3) an "output" procedure, which does whatever is desired
!> =============================================================================

MODULE experiment_mod
  USE item_mod
  USE task_mod
  USE task_list_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Override the init and output procedures, but keep the update and print
  !-----------------------------------------------------------------------------
  type, public, extends(tasks_t):: experiment_t
  contains
    procedure:: init
    procedure:: update
    procedure:: output
  end type
  !-----------------------------------------------------------------------------
  ! Override the init procedure, and keep everything else
  !-----------------------------------------------------------------------------
  type, public, extends(task_list_t):: experiment_list_t
  contains
    procedure:: init => init_experiment
  end type
  integer:: id=0
contains
!> =============================================================================
subroutine init (self)
  class(experiment_t):: self
  !.............................................................................
  call self%item_t%init                       ! generic item initialization
  print *,'experiment%init', self%id          ! just acknowledge, for now
end subroutine init
!> =============================================================================
subroutine output (self)
  class(experiment_t):: self
  !.............................................................................
  call self%print                             ! just a generic print, for now
end subroutine
!> =============================================================================
subroutine update (self)
  class(experiment_t):: self
  !.............................................................................
  call self%task_t%update
  print*, 'exeriment%update:', self%id, self%time
end subroutine
!> =============================================================================
SUBROUTINE init_experiment (list, name)
  class(experiment_list_t):: list
  character(len=*):: name
  type(experiment_t), pointer:: task
  class(item_t), pointer:: item
  integer:: i, n=5
  !.............................................................................
  call list%list_t%init (name)                ! generic initialization
  do i=1,n
    allocate (task)
    call task%init
    item => task
    call list%append (item)
  end do
  call list%task_list_t%init (name)           ! task list specific init
  call list%print
  print*,'XX0'
END SUBROUTINE
END MODULE
