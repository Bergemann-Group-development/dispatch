 !> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for tasks.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE task_mod
  USE item_mod
  USE list_mod
  USE bits_mod
  implicit none
  integer, parameter:: ntime=7
  !-----------------------------------------------------------------------------
  ! Basic task info
  !-----------------------------------------------------------------------------
  type, extends(item_t):: task_t
    integer:: rank = 0
    real(8):: position(3) = 0d0
    real(8):: velocity(3) = 0d0
    real(8):: size(3) = 0d0
    integer:: n_check = 0
    type(list_t), pointer:: nbor_list
  contains
  end type
  !-----------------------------------------------------------------------------
  ! Time slot info and procedures added in tasks_t
  !-----------------------------------------------------------------------------
  type, extends(task_t):: tasks_t
    real(8):: t(ntime)
    real(8):: dt(ntime)
    integer:: iit(ntime)
    integer:: nt = ntime
    integer:: it = 0
    integer:: verbose = 0
  contains
    procedure:: execute
    procedure:: check
    procedure:: set_virtual
    procedure:: send_to
    procedure:: overlaps
  end type
PRIVATE
PUBLIC ntime, task_t, tasks_t
CONTAINS

!> =============================================================================
!> =============================================================================
FUNCTION overlaps (self, item)
  class(tasks_t):: self
  class(item_t):: item
  logical:: overlaps
  overlaps = .true.
END FUNCTION overlaps

SUBROUTINE set_virtual (self)
  class(tasks_t):: self
END SUBROUTINE

SUBROUTINE send_to (self, rank)
  class(tasks_t):: self
  integer:: rank
END SUBROUTINE

!> =============================================================================
!> This method specifies the steps to perform on a task ready for execution:
!> After taking a time step it broadcasts the new state to neighbors, checks if
!> any of them can be moved to the ready_queue, and finally checks itself.
!> =============================================================================
SUBROUTINE execute(self)
  class(tasks_t):: self
  class(item_t), pointer:: item
  class(link_t), pointer:: nbor
  class(task_t), pointer:: nbtask
  !.............................................................................
  print*,'XX1'
  call self%update                        ! take a timestep
  nbor => self%nbor_list%head             ! start on nbor_list
  item => nbor%item
  select type (item)
  class is (task_t)
    nbtask => item                        ! neighbor task
  end select
  do while (associated(nbor))             ! continue to the end
    if (nbtask%rank == self%rank) then    ! same rank task
      call nbtask%check                  ! nbor task to ready_queue?
    else                                  ! neighbor on different rank
      select type (nbtask)
      type is (tasks_t)
        call nbtask%send_to (nbtask%rank) ! send to neighbor rank
      end select
    end if
    nbor => nbor%next                     ! next nbor
  end do
  call self%check                         ! self to ready_queue?
END SUBROUTINE

!> =============================================================================
!> Simple check function
!> =============================================================================
SUBROUTINE check (self)
  class(tasks_t):: self
  class(list_t), pointer:: list
  class(link_t), pointer:: nbor, link
  class(item_t), pointer:: item
  !.............................................................................
  nbor => self%nbor_list%head             ! start on nbor_list
  do while (associated(nbor))             ! continue to the end
    item => nbor%item
    if (item%time < self%time) then       ! fail if neighbor behind
      select type (item)
      type is (tasks_t)
        item%n_check = item%n_check + 1   ! count failures
      end select
      return
    end if
    nbor => nbor%next                     ! next nbor
  end do
!allocate (link)
!  link%item => self
!  call list%add_by_time (link)            ! add to list
END SUBROUTINE

END MODULE
