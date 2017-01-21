!> =============================================================================
!> The duty of the `experiment_mod` module is to define a generic (or specific)
!> `experiment_t` object, which has an initialisation procedure and a task list.
!> The tasks can be of any type belonging to the task class.
!> =============================================================================
!> This particular `experiment_mod` sets up a rigid sphere within a cylindrical
!> annulus that is rotating at the Keplerian speed of its radial centre. The
!> sphere also moves with the local Keplerian speed, but the gas is slightly
!> sub-Keplerian, resulting in a headwind.
!> =============================================================================
MODULE setup_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE list_mod
  USE task_list_mod
  USE experiment_mod
  USE components_mod
  USE scaling_mod
  USE planetesimal_mod
  USE parameters_mod
  USE rubiks_mod
  implicit none
  private

  real(8), parameter :: pi = asin(1.0d0) * 2.0d0

  ! components of the simulation
  ! a `setup_t` object is required to tie everything together
  type, public :: setup_t
    type(task_list_t), pointer :: task_list => null()
    character(len=64) :: experiment_name = 'planetesimal annulus'
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
  class(task_t), pointer :: rootpatch
  class(list_t), pointer :: task_list
  type(annulus_t), pointer :: annulus
  type(planetesimal_t), pointer :: planetesimal
  type(rubiks_t), pointer :: nested_grids, ring_grids
  logical, save :: first_time = .true.
  !-----------------------------------------------------------------------------
  call trace_begin ('setup_t%init')

  call io%init (self%experiment_name)
  call scaling%init()

  ! read input params
  if (first_time) then
    first_time = .false.
    call experiment_params%init()
  end if

  ! Create an empty task list.
  allocate(self%task_list)
  
  !-----------------------------------------------------------------------------
  ! Start with the the planetesimal; for the moment, it is a global data structure
  ! rather than a task.
  !-----------------------------------------------------------------------------
  allocate (planetesimal)
  planetesimal%mstar = experiment_params%mstar
  call planetesimal%init()
  ! point to the planetesimal task so that its parameters are globally accessible
  ! FIXME: This should be a list for multiple planets...
  experiment_params%planetesimal => planetesimal

  !-----------------------------------------------------------------------------
  ! Next, a cylindrical annulus for the computational domain, possibly made of
  ! many patches.
  !-----------------------------------------------------------------------------
  allocate (annulus)
  call annulus%init()
  if (associated(annulus%task_list)) then
    task_list => annulus%task_list
    call self%task_list%append_list (task_list)
  end if

  !-----------------------------------------------------------------------------
  ! Set-up a ring of patches with centred on the planetesimal's orbital radius
  !-----------------------------------------------------------------------------
  rootpatch => annulus%task_list%head%task
  allocate (ring_grids)
  ring_grids%rotating = annulus%rotating
  call ring_grids%init_ring(rootpatch%level)
  if (associated(ring_grids%task_list)) then
    task_list => ring_grids%task_list
    call self%task_list%append_list(task_list)
  end if

  !-----------------------------------------------------------------------------
  ! Set-up a series of "Rubik's cube"-like wedges centred on the planetesimal
  ! position
  !-----------------------------------------------------------------------------
  if (associated(ring_grids%task_list)) then
    rootpatch => ring_grids%task_list%head%task
  end if
  allocate (nested_grids)
  nested_grids%rotating = annulus%rotating
  call nested_grids%init_hierarchy(1,rootlvl=rootpatch%level)
  if (associated(nested_grids%task_list)) then
    task_list => nested_grids%task_list
    call self%task_list%append_list(task_list)
  end if

  !-----------------------------------------------------------------------------
  ! Done. Perform final list initialisations.
  !-----------------------------------------------------------------------------
  call self%task_list%init(trim(self%experiment_name)//' task list')
  call trace_end
END SUBROUTINE init_setup

END MODULE setup_mod
