!===============================================================================
!> Module to collect a block of tasks together before updating them in parallel
!===============================================================================
MODULE task_block_mod
  USE io_mod
  USE trace_mod
  USE link_mod
  implicit none
  private
  type links_t
    class(link_t), pointer:: link => null()
  end type
  type, public:: task_block_t
    integer:: tasks_per_block = 32
    integer:: tasks = 0
    type(links_t), pointer:: links(:)
  contains
    procedure:: init
    procedure:: update
    procedure:: run
  end type
!  type(task_block_t), public:: task_block
CONTAINS

!===============================================================================
!> Initialize a task block
!===============================================================================
SUBROUTINE init (self)
  class(task_block_t):: self
  !.............................................................................
  call trace%begin ('task_block_t%init')
  allocate (self%links(self%tasks_per_block))
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Fill up a task block
!===============================================================================
SUBROUTINE update (self, link)
  class(task_block_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  call trace%begin ('task_block_t%update')
  self%tasks = self%tasks + 1
  self%links(self%tasks)%link => link
  if (self%tasks == self%tasks_per_block) then
    call self%run
    self%tasks = 0
  end if
  call trace%end ()
END SUBROUTINE update

!===============================================================================
!> Run the tasks in a task block
!===============================================================================
SUBROUTINE run (self)
  class(task_block_t):: self
  !.............................................................................
  integer:: i
  class(link_t), pointer:: link
  !-----------------------------------------------------------------------------
  call trace%begin ('task_block_t%run')
  !$omp distribute
  do i=1,self%tasks_per_block
    link => self%links(i)%link
    call link%task%run
  end do
  call trace%end ()
END SUBROUTINE run

END MODULE task_block_mod