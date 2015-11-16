MODULE experiment_mod
  USE task_mod
  USE task_list_mod
  implicit none
  private
  type, public, extends(tasks_t):: experiment_t
  contains
    procedure:: init
    procedure:: print
  end type
  integer:: id=0
contains

subroutine init (self)
  class(experiment_t):: self
  id = id+1
  self%id = id
  print *,'experiment%init', self%id
end subroutine init

subroutine print (self)
  class(experiment_t):: self
  print*,'experiment%print:', self%id, self%time
end subroutine
END MODULE experiment_mod

PROGRAM test
  USE task_list_mod
  USE task_mod
  USE item_mod
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
