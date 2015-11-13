!> =============================================================================
!> $Id$
!> A task_list_t has the distinction that it must be able to call an initial
!> condition routine, which only exists in experiment_t tasks
!> =============================================================================
MODULE task_list_mod
  USE list_mod
  USE item_mod
  USE patch_mod
  implicit none
  type, extends(list_t):: task_list_t
  contains
    procedure:: process
  end type
CONTAINS

!> =============================================================================
!> Process the task list, computing initial values where appropriate (which
!> can take significant time).  The list process is best done by a single
!> thread, while the other ones can help with the actual computation.
!> =============================================================================
SUBROUTINE process (list, ready_queue)
  class(task_list_t):: list
  type(task_list_t), pointer:: ready_queue
  class(link_t), pointer:: link
  class(item_t), pointer:: item
  !.............................................................................
  !$omp single
  allocate (ready_queue)
  link => list%head
  do while (associated(link))
    item => link%item
    !$omp task
    select type (item)
    class is (patch_t)
      call item%init
    end select
    call ready_queue%append (item)
    !$omp end task
    link => link%next
  end do
  !$omp end single
END SUBROUTINE

END MODULE
