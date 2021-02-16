!===============================================================================
!> Each experiment must have a tasks_mod.f90, which is capable of allocating
!> all types of tasks used in that experiment, given only its task%type.  The
!> allocated task will be %pre_init()'ed and will be added to the task_list.
!>
!> A generic tasks_mod.f90 is available in the tasks/ directory, which covers
!> the case where only experiment_t task types are needed.
!===============================================================================
MODULE tasks_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE experiment_mod
  type, public:: tasks_t
  contains
    procedure:: allocate
  end type
  type(tasks_t), public:: tasks
CONTAINS

!===============================================================================
!> Allocate a task of a specific type
!===============================================================================
SUBROUTINE allocate (self, type, task)
  class(tasks_t)               :: self
  character(16)                :: type
  class(task_t), pointer       :: task
  !-----------------------------------------------------------------------------
  call trace%begin ('tasks_t%pre_init')
  if (type=='experiment_t') then
    allocate (experiment_t:: task)
  else if (type=='task_t') then
    allocate (experiment_t:: task)
  end if
  call trace%end ()
END SUBROUTINE allocate

END MODULE tasks_mod