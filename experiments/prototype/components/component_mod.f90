!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Generic component module.
!>
!> A component has a task list, mass, and position, velocity, and time offsets.
!> All offsets are to be added on top of higher order components, which they are
!> members of.
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE component_mod
  USE task_list_mod
  USE task_mod
  implicit none
  private
  type, public:: component_t
    type(task_list_t), pointer:: task_list
    character(len=120):: name
    real(8):: size(3)=1.0, position(3)=0.0, velocity(3)=0.0, time=0.0, mass=1.0
    integer:: ntask
  contains
    procedure:: init, add
  end type
contains
!===============================================================================
!> Initialize the task list and its tasks
!===============================================================================
SUBROUTINE init (self, name)
  class(component_t):: self
  character(len=*):: name
  class(task_t), pointer:: task
  integer:: i
  !.............................................................................
  self%name = name
  allocate (self%task_list)
  call self%task_list%list_t%init ('component task list')
  print *,'component '//trim(self%name)//'  has mean density ', self%mass/self%size**3
  do i=1,self%ntask
    allocate (task)
    call task%init
    call self%add (task)
  end do
END SUBROUTINE init
!===============================================================================
!> Add a task to the task list.  This is generic, but the task specified by the
!> caller is specific.  From some context, a specific component is being built,
!> by adding specific tasks to this task list.
!===============================================================================
SUBROUTINE add (self, task)
  class(component_t):: self
  class(task_t), pointer:: task
  !.............................................................................
  task%position = task%position + self%position
  task%velocity = task%velocity + self%velocity
  task%time     = task%time     + self%time
  call self%task_list%append (task)
END SUBROUTINE add
END MODULE component_mod
