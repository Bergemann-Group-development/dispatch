!> =============================================================================
!> $Id$
!> Skeleton experiment, which does nothing, but has the necessary procedures
!> =============================================================================
MODULE experiment_mod
  USE list_mod
  USE mhd_mod
  type experiment_t
  contains
    procedure, nopass:: init
  end type
CONTAINS

!> =============================================================================
!> Populate a task_list (complete list of all tasks relevant to this rank --
!> local tasks and virtual task = copies of tasks from neighbor ranks), and a
!> ready_list (list of tasks lokal to this rank, which are ready to be executed.
!> =============================================================================
SUBROUTINE init (task_list, ready_queue)
  class(list_t):: task_list, ready_queue
END SUBROUTINE

!> =============================================================================
!> Update a local task that is ready to be updated; all guard zone values have
!> been loaded, and all other preparations have been made.
!> =============================================================================
SUBROUTINE update (task)
  class(mhd_t):: task
END SUBROUTINE

END MODULE experiment_mod
