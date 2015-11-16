!> =============================================================================
!> $Id$
!> A task_list_t has the distinction that it must be able to call an initial
!> condition routine, which only exists in experiment_t tasks
!> =============================================================================
MODULE task_list_mod
  USE list_mod
  USE item_mod
  USE task_mod
  USE bits_mod
  implicit none
  private
  type, public, extends(list_t):: task_list_t
    real(8):: time_end
    type(task_list_t), pointer:: ready_queue
  contains
    procedure:: init
    procedure:: check
    procedure:: process
    procedure:: apply
    procedure:: print
  end type
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (list, name)
  class(task_list_t):: list
  character(len=*):: name
  !.............................................................................
  type(task_list_t), pointer:: out
  type(link_t), pointer:: link
  class(item_t), pointer:: item
  !.............................................................................
  print*, 'task_list%init:', list%name
  allocate (list%ready_queue)
  call list%ready_queue%list_t%init ('ready queue')
  link => list%head
  do while (associated(link))
    item => link%item
    call list%ready_queue%append (item)
    link => link%next
  end do
  print*,'MK1'
  call list%ready_queue%print
END SUBROUTINE

!> =============================================================================
!> =============================================================================
SUBROUTINE check (list, item)
  class(task_list_t):: list
  !.............................................................................
  type(link_t), pointer:: link
  class(item_t), pointer:: item
  !.............................................................................
  print*, 'task_list%check:', list%name, item%id, item%time
  if (item%time > self%time_end) then
    call self%remove (item)
  end if
END SUBROUTINE

!> =============================================================================
!> =============================================================================
SUBROUTINE print (list)
  class(task_list_t):: list
  !.............................................................................
  type(task_list_t), pointer:: out
  type(link_t), pointer:: link
  class(item_t), pointer:: item
  !.............................................................................
  print*, 'task_list%print:', list%name
  link => list%head
  do while (associated(link))
    item => link%item
    call item%print
    link => link%next
  end do
END SUBROUTINE

!> =============================================================================
!> Process the task list, computing initial values where appropriate (which
!> can take significant time).  The list process is best done by a single
!> thread, while the other ones can help with the actual computation.
!> =============================================================================
SUBROUTINE process (list, queue)
  class(task_list_t):: list
  type(task_list_t), pointer:: queue
  class(link_t), pointer:: link
  class(item_t), pointer:: item
  type(tasks_t), pointer:: task
  !.............................................................................
  !$omp single
  allocate (queue)
  link => list%head
  do while (associated(link))
    item => link%item
    !$omp task
    select type (item)
    class is (task_t)
      if (item%is_set(bits%internal+bits%virtual)) then
        call item%init
      end if
    end select
    if (item%is_set(bits%internal)) then
      call queue%append (item)
    end if
    !$omp end task
    link => link%next
  end do
!  link => list%head
!  do while (associated(link))
!    item => link%item
!    select type (item)
!    type is (tasks_t)
!      nbor_list => task%nbor_list
!      call nbor_list%init (task_list, task)
!    end select
!    link => link%next
!  end do
  !$omp end single
END SUBROUTINE
!> =============================================================================
!> Process the tasks in a task list through a procedure, which optionally adds
!> tasks to an out list
!> =============================================================================
SUBROUTINE apply (self, procedure, out)
  class(task_list_t):: self
  type(task_list_t), pointer:: out
  type(link_t), pointer:: link
  class(item_t), pointer:: task
  !.............................................................................
  !$omp single
  link => self%head
  do while (associated(link))
    task => link%item
    !$omp task
    call procedure (task, out)
    !$omp end task
    link => link%next
  end do
  !$omp end single
END SUBROUTINE

END MODULE
