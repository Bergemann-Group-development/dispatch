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
  USE trace_mod
  USE mpi_mod
  USE io_mod
  implicit none
  private
  type, public, extends(list_t):: task_list_t
    real(8):: time_end
    type(task_list_t), pointer:: ready_queue
  contains
    procedure:: init
    procedure:: update
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
  type(task_list_t), pointer:: out
  type(link_t), pointer:: link
  class(item_t), pointer:: item
  real(8):: time_end
  namelist /task_list_params/ time_end
  !.............................................................................
  call trace_begin('task_list_t%init '//list%name)
  rewind (io%input)
  read (io%input, task_list_params)
  if (io%master) write (*, task_list_params)
  list%time_end = time_end
  !
  allocate (list%ready_queue)
  call list%ready_queue%list_t%init ('ready queue')
  link => list%head
  do while (associated(link))
    item => link%item
    call list%ready_queue%append (item)
    link => link%next
  end do
  call list%ready_queue%print
  call trace_end
END SUBROUTINE

!> =============================================================================
!> =============================================================================
SUBROUTINE update (list)
  class(task_list_t):: list
  !.............................................................................
  type(link_t), pointer:: link
  class(item_t), pointer:: item
  !.............................................................................
  call trace_begin('task_list_t%update '//list%name)
  link => list%ready_queue%head
  item => link%item
  call list%ready_queue%remove_link (link)
  call item%update
  call list%check (item)
  call trace_end
END SUBROUTINE

!> =============================================================================
!> =============================================================================
SUBROUTINE check (list, item)
  class(task_list_t):: list
  !.............................................................................
  type(link_t), pointer:: link
  class(item_t), pointer:: item
  !.............................................................................
  call trace_begin('task_list_t%check '//list%name)
  print*, 'task_list%check:', trim(list%name), item%id, item%time, list%time_end
  if (item%time > list%time_end) then
    call list%remove (item)
    call list%ready_queue%remove (item)
  end if
  call list%ready_queue%add_by_time (item)
  call trace_end
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
  call trace_begin('task_list_t%process '//list%name)
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
  call trace_end
END SUBROUTINE
!> =============================================================================
!> Process the tasks in a task list through a procedure, which optionally adds
!> tasks to an out list
!> =============================================================================
SUBROUTINE apply (list, procedure, out)
  class(task_list_t):: list
  type(task_list_t), pointer:: out
  type(link_t), pointer:: link
  class(item_t), pointer:: task
  !.............................................................................
  call trace_begin('task_list_t%apply '//list%name)
  !$omp single
  link => list%head
  do while (associated(link))
    task => link%item
    !$omp task
    call procedure (task, out)
    !$omp end task
    link => link%next
  end do
  !$omp end single
  call trace_end
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
END MODULE
