 !> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for tasks.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE task_mod
  USE item_mod
  USE list_mod
  USE bits_mod
  USE trace_mod
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
    type(list_t), pointer:: nbor_list => null()
  contains
  end type
  !-----------------------------------------------------------------------------
  ! Time slot info and procedures added in tasks_t
  !-----------------------------------------------------------------------------
  type, extends(task_t):: tasks_t
    real(8):: t(ntime) = 0d0
    real(8):: dt(ntime) = 1d0
    integer:: iit(ntime) = 1
    integer:: nt = ntime
    integer:: it = 1
    integer:: verbose = 1
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
  call trace_begin ('task_t%overlap')
  overlaps = .true.
  call trace_end
END FUNCTION overlaps

SUBROUTINE set_virtual (self)
  class(tasks_t):: self
  call trace_begin ('task_t%set_virtual')
  call trace_end
END SUBROUTINE

SUBROUTINE send_to (self, rank)
  class(tasks_t):: self
  integer:: rank
  call trace_begin ('task_t%send_to')
  call trace_end
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
  call trace_begin ('task_t%execute')
  call self%update                        ! take a timestep
  nbor => self%nbor_list%head             ! start on nbor_list
  do while (associated(nbor))             ! continue to the end
    item => nbor%item
    select type (item)
    class is (task_t)                     ! select task_t items
      nbtask => item                      ! neighbor task
    end select
    if (nbtask%rank == self%rank) then    ! same rank task
      call nbtask%check                   ! nbor task to ready_queue?
    else                                  ! neighbor on different rank
      select type (nbtask)
      type is (tasks_t)
        call nbtask%send_to (nbtask%rank) ! send to neighbor rank
      end select
    end if
    nbor => nbor%next                     ! next nbor
  end do
  call self%check                         ! self to ready_queue?
  call trace_end
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
  call trace_begin ('task_t%check')
  nbor => self%nbor_list%head             ! start on nbor_list
  do while (associated(nbor))             ! continue to the end
    item => nbor%item
    if (item%time < self%time) then       ! fail if neighbor behind
      select type (item)
      type is (tasks_t)
        item%n_check = item%n_check + 1   ! count failures
      end select
      call trace_end
      return
    end if
    nbor => nbor%next                     ! next nbor
  end do
!allocate (link)
!  link%item => self
!  call list%add_by_time (link)            ! add to list
  call trace_end
END SUBROUTINE

END MODULE
