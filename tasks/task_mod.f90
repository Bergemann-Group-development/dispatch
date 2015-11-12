!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for task
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE task_mod
  USE item_mod
  USE list_mod
  implicit none
  integer, parameter:: ntime=7
  !-----------------------------------------------------------------------------
  ! Basic task info
  !-----------------------------------------------------------------------------
  type, abstract, extends(item_t):: task_t
    integer:: id = 0
    integer:: rank = 0
    real(8):: dtime = 0d0
    real(8):: position(3) = 0d0
    real(8):: velocity(3) = 0d0
    real(8):: size(3) = 0d0
    integer:: status = 0
    type(list_t), pointer:: nbor_list
  contains
    procedure                :: execute => execute_task
    procedure(one), deferred :: check
    procedure(two), deferred :: send_to    ! deferred to specific task
  end type
  abstract interface
    subroutine one (self)
    import task_t
    class(task_t):: self
    end subroutine
    subroutine two (self, rank)
    import task_t
    class(task_t):: self
    integer:: rank
    end subroutine
  end interface
  !-----------------------------------------------------------------------------
  ! Timeslot info added
  !-----------------------------------------------------------------------------
  type, extends(task_t):: tasks_t
    real(8):: t(ntime)
    real(8):: dt(ntime)
    integer:: iit(ntime)
    integer:: nt = ntime
    integer:: it = 0
    integer:: verbose = 0
  contains
    procedure:: init => void
    procedure:: check => void
    procedure:: update => void
    procedure:: send_to => void2
    procedure:: allocate => void
    procedure:: deallocate => void
    procedure:: print => void
  end type
PRIVATE
PUBLIC ntime, task_t, tasks_t
CONTAINS

SUBROUTINE void (self)
  class(tasks_t):: self
END SUBROUTINE

SUBROUTINE void2 (self, rank)
  class(tasks_t):: self
  integer:: rank
END SUBROUTINE

!> =============================================================================
!> This method specifies the steps to perform on a task ready for execution:
!> After taking a time step it broadcasts the new state to neighbors, checks if
!> any of them can be moved to the ready_queue, and finally checks itself.
!> =============================================================================
SUBROUTINE execute_task (self)
  class(task_t):: self
  class(item_t), pointer:: item
  class(link_t), pointer:: nbor
  class(task_t), pointer:: nbtask
  !.............................................................................
  call self%update                        ! take a timestep
  nbor => self%nbor_list%head             ! start on nbor_list
  item => nbor%item
  select type (item)
  class is (task_t)
    nbtask => item                        ! neighbor task
  end select
  do while (associated(nbor))             ! continue to the end
    if (nbtask%rank == self%rank) then    ! same rank task
      call nbtask%check                   ! nbor task to ready_queue?
    else                                  ! neighbor on different rank
      call nbtask%send_to (nbtask%rank)   ! send to neighbor rank
    end if
    nbor => nbor%next                     ! next nbor
  end do
  call self%check                         ! self to ready_queue?
END SUBROUTINE

END MODULE
