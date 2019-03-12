!===============================================================================
!> $Id$
!> Prototype experiment, showing the steps and components needed:
!> 1) a "create" procedure, which creates tasks and adds them to a task_list
!> 2) an "init" proc, which will be executed by the task_list%init
!> 3) an "update" procedure, which will be executed by task_queue%update
!> 3) an "output" procedure, which does whatever is desired
!===============================================================================

MODULE experiment_mod
  USE task_mod
  USE task_list_mod
  implicit none
  private
  type, public, extends(task_t):: experiment_t
  contains
    procedure:: init
    procedure:: output
  end type
  integer:: id=0
contains

subroutine init (self)
  class(experiment_t):: self
  id = id+1
  self%id = id
  print *,'experiment%init', self%id
end subroutine init

subroutine output (self)
  class(experiment_t):: self
  call self%print
end subroutine

END MODULE experiment_mod

PROGRAM test
  USE item_mod
  USE task_list_mod
  USE experiment_mod
  implicit none
  type(task_list_t):: task_list
  type(experiment_t), pointer:: task
  class(item_t), pointer:: item
  integer:: i, n=5
  !.............................................................................
  call task_list%init ('test')
  do i=1,n
    allocate (task)
    call task%init
    item => task
    call task_list%append (item)
  end do
  call task_list%print
END PROGRAM test
