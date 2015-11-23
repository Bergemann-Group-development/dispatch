!> =============================================================================
!> The duty of the experiment_mod module is to define a generic experiment_t
!> object, which has an init procedure and a task list.  The tasks can be of
!> any type belonging to the task class.
!> =============================================================================
MODULE experiment_mod
  USE envelope_mod
  USE planets_mod
  USE pebbles_mod
  USE star_mod
  USE task_list_mod
  USE task_mod
  USE trace_mod
  implicit none
  private
  type, public:: experiment_t
    type(task_list_t), pointer:: task_list
  contains
    procedure:: init
  end type
contains

!> =============================================================================
!> The init procedure just splices together the relevant task lists.  The
!> specific behavior of each component is described in its respective module
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  type(star_t), pointer:: star
  type(pebbles_t), pointer:: pebbles
  type(planets_t), pointer:: planets
  type(envelope_t), pointer:: envelope
  class(task_t), pointer:: task
  class(task_list_t), pointer:: task_list
  !-----------------------------------------------------------------------------
  ! Start with the envelope%task_list
  !-----------------------------------------------------------------------------
  call trace_begin ('experiment%init')
  allocate (envelope)
  call envelope%init ('env')
  self%task_list => envelope%task_list
  !-----------------------------------------------------------------------------
  ! Create a star object, initialize it, and add its task to the task list
  !-----------------------------------------------------------------------------
  allocate (star)
  call star%init
  task => star
  call self%task_list%append (task)                     ! a star is a task
  !-----------------------------------------------------------------------------
  ! Create a planets component, initialize it, and append its task list to the
  ! existing task list
  !-----------------------------------------------------------------------------
  allocate (planets)
  call planets%init
  task_list => planets%task_list
  call self%task_list%append_list (task_list)   ! planets is a component
  !-----------------------------------------------------------------------------
  ! Create a pebble component, initialize it, and append
  !-----------------------------------------------------------------------------
  allocate (pebbles)
  call pebbles%init
  task_list => pebbles%task_list
  call self%task_list%append_list (task_list)      ! pebble is a component
  call self%task_list%init
  call trace_end
END SUBROUTINE init

END MODULE experiment_mod
