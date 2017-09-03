!> =============================================================================
!> The duty of the `experiment_mod` module is to define a generic (or specific)
!> `experiment_t` object, which has an initialisation procedure and a task list.
!> The tasks can be of any type belonging to the task class.
!> =============================================================================
!> This particular `experiment_mod` sets up a rigid cylinder and a shock that
!> subsequently passes, and is deflected by, the rigid cylinder.  The Newton-Euler
!> eqns. are used to describe the motion of the cylinder, such that the shock
!> can actually impulsively move the cylinder.
!> =============================================================================
MODULE setup_mod
  USE io_mod
  USE trace_mod
  USE list_mod
  USE task_list_mod
  USE experiment_mod
  USE components_mod
  USE scaling_mod
  implicit none
  private

  real(8), parameter :: pi = 3.14159265358979323846

  ! components of the simulation
  ! an `setup_t` object is required to tie everything together
  type, public :: setup_t
    type(task_list_t), pointer :: task_list
    character(len=64) :: experiment_name = 'shocked_cylinder'
  contains
    procedure :: init => init_setup
  end type
  type(setup_t), public :: setup

CONTAINS

!> =============================================================================
!> The init procedure just splices together the relevant task lists.  The
!> specific behavior of each component is described in its respective module
!> =============================================================================
SUBROUTINE init_setup (self)
  class(setup_t) :: self
  class(experiment_t), pointer :: task
  class(list_t), pointer :: task_list
  type(box_t), pointer :: box
  type(solid_cylinder_t), pointer :: cylinder
  type(envelope_annulus_t), pointer :: envelope
  !-----------------------------------------------------------------------------
  call trace_begin ('setup_t%init')

  call io%init (self%experiment_name)
  call scaling%init()

  ! Create an empty task list.
  allocate(self%task_list)
  !-----------------------------------------------------------------------------
  ! Start with a box for the computational domain, possible made of many patches
  !-----------------------------------------------------------------------------
  allocate (box)
  call box%init()
  task_list => box%task_list
  call self%task_list%append_list (task_list)

  !-----------------------------------------------------------------------------
  ! Next, the patch represented a solid cylinder; only one patch and therefore
  ! only one task
  !-----------------------------------------------------------------------------
  allocate (cylinder)
  call cylinder%init()
  task_list => cylinder%task_list
  call self%task_list%append_list (task_list)

  !-----------------------------------------------------------------------------
  ! Finally, the "envelope" component, which surrounds the solid cylinder, and
  ! may be composed of several patches.
  !-----------------------------------------------------------------------------
  allocate (envelope)
  call envelope%init(cylinder)
  task_list => envelope%task_list
  call self%task_list%append_list (task_list)

  !-----------------------------------------------------------------------------
  ! Done. Perform final list initialisations.
  !-----------------------------------------------------------------------------
  call self%task_list%init(trim(self%experiment_name)//' task list')
  call trace_end
END SUBROUTINE init_setup

END MODULE setup_mod
