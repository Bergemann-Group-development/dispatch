MODULE refine_mod
  USE list_mod
  USE link_mod
  USE task_mod

  type:: refine_t
    logical:: on                              ! allow external access
    integer:: levelmin=0                      ! Min level of refinment
    integer:: levelmax=0                      ! Max level of refinment
  contains
    procedure :: check
  end type
  type(refine_t), public:: refine

contains
  !===============================================================================
  !> Check if refinement is desired and needed for the current task; if so, push
  !> new tasks onto the (ready) queue with the same time as the parent, i.e., to
  !> the head of the queue.  If the task is a leaf task, and the level of
  !> refinement does not appear to be needed, request derefinement.
  !===============================================================================
  SUBROUTINE check (self, tasklist, link, was_refined, was_derefined)
    class(refine_t):: self
    class(list_t), target:: tasklist
    class(link_t), pointer:: link
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: task
    logical, intent(out):: was_refined, was_derefined
    !.............................................................................
    was_refined = .false.
    was_derefined = .false.
  END SUBROUTINE check
END MODULE refine_mod