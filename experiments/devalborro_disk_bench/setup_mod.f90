!===============================================================================
!> The duty of the `experiment_mod` module is to define a generic (or specific)
!> `experiment_t` object, which has an initialisation procedure and a task list.
!> The tasks can be of any type belonging to the task class.
!===============================================================================
MODULE setup_mod
  USE io_mod
  USE trace_mod
  USE list_mod
  USE task_list_mod
  USE task_mod
  USE experiment_mod
  USE components_mod
  USE parameters_mod
  USE scaling_mod
  implicit none
  private

  real(8), parameter :: pi = 3.14159265358979323846

  ! components of the simulation
  ! an `setup_t` object is required to tie everything together
  type, public :: setup_t
    type(task_list_t), pointer :: task_list
    character(len=64) :: experiment_name = 'de Val-Borro benchmark'
  contains
    procedure :: init => init_setup
  end type
  type(setup_t), public :: setup

CONTAINS

!===============================================================================
!> The init procedure just splices together the relevant task lists.  The
!> specific behavior of each component is described in its respective module
!===============================================================================
SUBROUTINE init_setup (self)
  class(setup_t) :: self
  class(experiment_t), pointer :: task
  class(list_t), pointer :: task_list
  type(annulus_t), pointer :: annulus
  class(task_t), pointer :: rootpatch
  logical, save :: first_time = .true.
  logical :: ramp_mass
  real(8) :: mtotal, aplanet, mplanet, gpp_soften, aspect_rat, phiplanet
  real(4) :: kill_rmin, kill_rmax
  namelist / IC_params / mtotal, aplanet, mplanet, gpp_soften, kill_rmin, &
                              kill_rmax, aspect_rat, ramp_mass, phiplanet
  !-----------------------------------------------------------------------------
  call trace_begin ('setup_t%init')
  call io%header('begin setup_t%init')

  ! read input params
  if (first_time) then
    first_time = .false.
    ! default values
    mtotal  = experiment_params%mtotal
    aplanet = experiment_params%aplanet
    mplanet = experiment_params%mplanet
    gpp_soften = experiment_params%gpp_soften
    kill_rmin  = experiment_params%kill_rmin
    kill_rmax  = experiment_params%kill_rmax
    aspect_rat = experiment_params%aspect_rat
    ramp_mass = experiment_params%ramp_mass
    phiplanet = experiment_params%phiplanet
    rewind (io%input); read (io%input, initial_params)
    if (io%master) then
      write (*, initial_params)
      write (io_unit%nml, initial_params)
      flush (io_unit%nml)
    end if
    experiment_params%mtotal = mtotal
    experiment_params%aplanet = aplanet
    experiment_params%mplanet = mplanet
    experiment_params%gpp_soften = gpp_soften
    experiment_params%mstar = mtotal - mplanet
    experiment_params%kill_rmin = kill_rmin
    experiment_params%kill_rmax = kill_rmax
    experiment_params%aspect_rat = aspect_rat
    experiment_params%ramp_mass = ramp_mass
    experiment_params%phiplanet = phiplanet
  end if

  ! Create an empty task list.
  allocate(self%task_list)
  !-----------------------------------------------------------------------------
  ! Start with a box for the computational domain, possibly made of many patches
  !-----------------------------------------------------------------------------
  allocate (annulus)
  call annulus%init(experiment_params)
  task_list => annulus%task_list
  call self%task_list%append_list (task_list)

  !-----------------------------------------------------------------------------
  ! at the moment, you can only have only one planet.
  ! the planet is not a task, but rather a global data structure.
  !-----------------------------------------------------------------------------
  call planet%init()

  !-----------------------------------------------------------------------------
  ! Done. Perform final list initialisations.
  !-----------------------------------------------------------------------------
  call self%task_list%init(trim(self%experiment_name)//' task list')
  call self%task_list%reset_status
  call self%task_list%count_status
  call self%task_list%init_bdries
  call self%task_list%init_all_nbors
  call self%task_list%print_tasks

  call trace_end
END SUBROUTINE init_setup

END MODULE setup_mod
