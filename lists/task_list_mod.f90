!> ==================-===========================================================
!> $Id$
!> A task_list_t has the distinction that it must be able to call an initial
!> condition routine, which only exists in experiment_t tasks.  Procedures in
!> the module use these foreign procedures:
!>
!>    task_t%init         : to initialize a task
!>    task_t%update       : to update a task
!>    task_t%send_to      : send data to foreign task
!>    list_t%remove_task  : remove task from list
!>    list_t%remove_head  : remove head from list
!>
!> =============================================================================
MODULE task_list_mod
  USE list_mod
  USE task_mod
  USE tasks_mod
  USE bits_mod
  USE trace_mod
  USE mpi_mod
  USE io_mod
  implicit none
  private
  type, public, extends(list_t):: task_list_t
    type(task_list_t), pointer:: ready_queue => null()
  contains
    procedure:: init
    procedure:: update
    procedure:: check_mpi
    procedure:: check_nbors
    procedure:: check_ready
!    procedure:: print
  end type
CONTAINS

!> =============================================================================
!> Initialize the task list, by allocating a ready queue and initializing the
!> task list tasks
!> =============================================================================
SUBROUTINE init (self, name)
  class(task_list_t):: self
  character(len=*):: name
  type(link_t), pointer:: link
  class(task_t), pointer:: task
  !.............................................................................
  call trace_begin('task_list_t%init '//name)
  allocate (self%ready_queue)
  call self%ready_queue%list_t%init ('ready queue')
  link => self%head
  do while (associated(link))
    task => link%task
    call task%init
    link => link%next
  end do
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Update the state of the task list, taking the steps necessary to update the
!> head task on the ready_queue, and check for consequences
!> =============================================================================
SUBROUTINE update (self)
  class(task_list_t):: self
  !.............................................................................
  type(link_t), pointer:: head
  class(task_t), pointer:: task
  !.............................................................................
  call trace_begin('task_list_t%update '//self%name)
  head => self%ready_queue%head                        ! head link
  task => head%task                                    ! head task
  print*,self%ready_queue%name
  if (associated(task)) then
    call self%check_mpi                                 ! check incoming MPI
    call self%ready_queue%remove_head                   ! remove the head link
    call task%update                                    ! update the task
    if (task%time > task%end_time) then                 ! finished:
      call self%remove_task (task)                      ! also remove from self
    else                                                ! not finished:
      call trace('9')
      call self%check_nbors (head)                      ! any nbor ready?
    end if
    deallocate (head)
  end if
  call trace_end
END SUBROUTINE update

!> =============================================================================
!> Check for MPI messages
!> =============================================================================
SUBROUTINE check_mpi (self)
  class(task_list_t):: self
  !.............................................................................
  call trace_begin ('check_mpi')
  call trace_end
END SUBROUTINE check_mpi

!> =============================================================================
!> Among a task and its neighbor tasks, move local tasks to ready_queue if they
!> are ready, and send task data to non-locals.
!> =============================================================================
SUBROUTINE check_nbors (self, link)
  class(task_list_t):: self
  type(link_t), pointer:: link, nbor
  class(task_t), pointer:: task, nbtask
  !.............................................................................
  call trace_begin('task_list_t%check_nbors '//self%name)
  task => link%task                                   ! main task
  nbor => link%nbor                                   ! first nbor
  do while (associated (nbor))                        ! keep going until end
    nbtask => nbor%task                               ! neighbor task
    if (nbtask%rank==mpi%rank) then                   ! local rank:
      call self%ready_queue%check_ready (nbor)        ! nbor ready?
    else                                              ! foreign rank:
      call task%send_to (nbtask)                      ! send task data
    end if
    nbor => nbor%next                                 ! next nbor
  end do
  call self%ready_queue%check_ready (link)            ! finally check link task
  call trace_end
END SUBROUTINE

!> =============================================================================
!> Check if the task of link is ready, and if so, add it to self
!> =============================================================================
SUBROUTINE check_ready (self, link)
  class(task_list_t):: self
  class(link_t):: link
  class(task_t), pointer:: task, nbtask
  type(link_t), pointer:: nbor
  !.............................................................................
  call trace_begin('task_list_t%check_ready '//self%name)
  task => link%task
  nbor => link%nbor                                   ! start on nbor list
  do while (associated (nbor))                        ! keep going until end
    nbtask => nbor%task
    if (task%ahead_of(nbtask)) return                 ! give up -- do nothing
    nbor => nbor%next                                 ! next nbor
  end do
  call self%add_by_time (task)                        ! add task to list
  call trace_end
END SUBROUTINE check_ready

END MODULE
