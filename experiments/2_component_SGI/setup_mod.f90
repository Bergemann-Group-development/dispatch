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
  USE rubiks_mod
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
  type(rubiks_t), pointer:: nested_grids
  class(task_t), pointer :: rootpatch
  logical, save :: first_time = .true.
  logical :: ramp_mass, wave_kill
  real(8) :: sigma0, gamma_c, mstar, dust2gas, T0, mu, gamma_T, amp, lambda, &
             nu, a, d_int
  real(4) :: kill_rmin, kill_rmax, roffset, R, Rc, RT
  integer :: nv, nc
  namelist / IC_params / R, sigma0, Rc, mstar, kill_rmin, kill_rmax,roffset, &
                              gamma_c, dust2gas, nv, nc, T0, RT, mu, gamma_T, &
                              wave_kill, amp, lambda, nu, a, d_int
  !-----------------------------------------------------------------------------
  call trace_begin ('setup_t%init')

  call io%init (self%experiment_name)

  ! read input params
  if (first_time) then
    ! default values
    R = experiment_params%R
    sigma0 = experiment_params%sigma0
    Rc = experiment_params%Rc
    gamma_c = experiment_params%gamma_c
    mstar  = experiment_params%mstar
    kill_rmin  = experiment_params%kill_rmin
    kill_rmax  = experiment_params%kill_rmax
    roffset = experiment_params%roffset
    dust2gas  = experiment_params%dust2gas
    nv = experiment_params%nv
    nc = experiment_params%nc
    T0 = experiment_params%T0
    RT = experiment_params%RT
    mu = experiment_params%mu
    gamma_T = experiment_params%gamma_T
    nu = experiment_params%nu
    a = experiment_params%a
    d_int = experiment_params%d_int
    wave_kill = experiment_params%wave_kill
    amp = experiment_params%amp
    lambda = experiment_params%lambda
    rewind (io%input); read (io%input, initial_params)
    if (io%master) write (*, initial_params)
    experiment_params%R = R
    experiment_params%mstar = mstar
    experiment_params%mu = mu
  end if
  
  call scaling%init()
  
  if (first_time) then
    first_time = .false.
    rewind (io%input); read (io%input, initial_params)
    if (io%master) write (*, initial_params)
    experiment_params%sigma0 = sigma0
    experiment_params%Rc = Rc
    experiment_params%gamma_c = gamma_c
    experiment_params%kill_rmin = kill_rmin
    experiment_params%kill_rmax = kill_rmax
    experiment_params%roffset = roffset
    experiment_params%dust2gas  = dust2gas
    experiment_params%T0 = T0
    experiment_params%RT = RT
    experiment_params%gamma_T = gamma_T
    experiment_params%nv = nv
    experiment_params%nc = nc
    experiment_params%wave_kill = wave_kill
    experiment_params%nu = nu
    experiment_params%a = a
    experiment_params%d_int = d_int
    experiment_params%amp = amp
    experiment_params%lambda = lambda
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
  ! Set-up a series of "Rubik's cube"-like wedges centred on the planetesimal
  ! position
  !-----------------------------------------------------------------------------
  rootpatch => annulus%task_list%head%task
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
  call self%task_list%init_all_nbors
  call trace_end
END SUBROUTINE init_setup

END MODULE setup_mod
