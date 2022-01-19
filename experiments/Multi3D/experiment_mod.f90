!===============================================================================
!> Top level MHD task definition -- can either rely entirely on lower levels,
!> including extras_t%initial_conditions, or can optionally do such things here
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE solver_mod
  USE units_mod, only: cgs
  USE scaling_mod, only: scaling
  USE kinds_mod
  USE composition_mod, only: composition_type
  USE atom_var, only: atom_type
  USE geometry_var, only:geometry_type
  USE spectrum_mod, only:spectrum_type
  USE m3d_patch_mod, only: m3d_t, halo_m3d_t
  USE opac_patch_mod, only: opacity_t
  USE rad_patch_mod, only: radiation_t
  USE atom_task_mod, only: atom_task_t
  USE rate_patch_mod, only: rate_patch_t
  USE list_mod
  USE link_mod
  USE task_mod
  USE math_mod
  USE timer_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
    type(composition_type) :: composition
    class(spectrum_type), pointer :: sp
    type(geometry_type) :: geometry
    class(opacity_t), pointer, dimension(:) :: back_opac
    class(opacity_t), pointer, dimension(:,:) :: opac
    class(radiation_t), pointer, dimension(:,:) :: rt_patches
    class(atom_task_t), pointer :: atom_task
    real(8), allocatable, dimension(:,:,:) :: ne, temp, rho, vx, vy, vz
    real(8), allocatable, dimension(:,:,:,:) :: nh
    integer :: idx_pe=6, n_mu, n_phi, n_nu, iter=0, maxiter=1
  contains
    procedure pre_init
    procedure init
    procedure update
    procedure init_level
    procedure bndries
    procedure pre_init_composition
    procedure init_composition
    procedure append_to
    procedure link_sub_task
    procedure link_and_init
    procedure init_sub_patches
    procedure freeze_mhd
  end type
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE pre_init (self)
  class(experiment_t), target :: self
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('experiment_mod::pre_init', itimer=itimer)
  call self%gpatch_t%pre_init ()
  !call self%pre_init_atom()
  call self%pre_init_composition()
  call trace_end (itimer=itimer)
END SUBROUTINE pre_init

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer, save:: itimer=0
  real(8), save:: timestep=0.01d0
  !----------------------------------------------------------------------------
  call trace_begin('experiment_mod::init', itimer=itimer)
!  call self%solver_t%init
  call self%gpatch_t%init                                           ! keep

  allocate(self%nh(self%n(1), self%n(2), self%n(3), 6))

  self%dtime = timestep

  call trace_end (itimer=itimer)
END SUBROUTINE init

!===============================================================================
!> Setup initial states in the patches
!===============================================================================
SUBROUTINE bndries (self)
  class(experiment_t):: self
  real(var_k), dimension(:,:,:,:,:), pointer :: mem
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%bndries')
  mem => self%mem(:,:,:,:,self%it,:)
  if (any(self%boundaries%topbot)) call self%boundaries%boundary(mem)
  call trace%end()
END SUBROUTINE bndries

!===============================================================================
!> This update only concerns the MHD task; the RT tasks are called directly
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  integer :: it_l, in_l, in_lte_l
  !----------------------------------------------------------------------------
  !call timer%tic ()
  !call trace_begin('experiment_t%update')
  !call self%solver_t%update

  if (self%iter == 0) then
    call self%init_composition()
    call self%freeze_mhd()
  endif

  self%iter = self%iter + 1
  !call trace_end
  !call timer%toc ('update took ', 1)
END SUBROUTINE update

!===============================================================================
SUBROUTINE init_level (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  if (any(self%n > 1)) then
    ! The level does not depend on the patch size in non-regular grids
    self%level = 1
  end if
END SUBROUTINE init_level

!===============================================================================
SUBROUTINE pre_init_composition (self)
  class(experiment_t) :: self
  character(len=80) :: absdat_file, abund_file
  type(composition_type), save :: composition
  logical, save :: first_time=.true.
  !-----------------------------------------------------------------------------
  if (first_time) then
    write (*,*) ' - Initializing the NEW opacity package'
    absdat_file = './input_multi3d/absdat'
    abund_file  = './input_multi3d/abund'
    call composition%init(absdat_file, abund_file)
    first_time = .false.
  endif

  ! copy composition with all the information from the abund
  ! and absdat file into dispackage to make it threadprivate
  self%composition = composition
  ! this is needed, otherwise all pointers will still point
  ! to ions and molecules in multi3d%composition
  call self%composition%reassign_pointer()
END SUBROUTINE pre_init_composition


!===============================================================================
SUBROUTINE init_composition (self)
  class(experiment_t) :: self
  !-----------------------------------------------------------------------------
  ! loading temp, ne and rho in to composition objekt
  ! ion and molecule equilibrium are solved immediately
  call self%composition%load_state_variables(self%temp, self%ne, self%rho)
  call self%composition%get_nhyd(self%nh)
END SUBROUTINE init_composition

!===============================================================================
SUBROUTINE append_to(self, task_list)
  class (experiment_t), target :: self
  class(list_t) :: task_list
  class(link_t),  pointer :: link
  class (task_t), pointer :: task
  !-----------------------------------------------------------------------------
  allocate (halo_m3d_t:: link)
  link%task => self
  self%link => link
  call task_list%append (link)
END SUBROUTINE append_to

!===============================================================================
SUBROUTINE init_sub_patches(self)
  class (experiment_t), target :: self
  logical, save :: first_time=.true.
  class(m3d_t), pointer :: patch
  class(atom_task_t), pointer :: atom_task
  class(opacity_t), pointer, dimension(:) :: back_opac
  class(opacity_t), pointer, dimension(:,:) :: sum_opac
  class(radiation_t), pointer, dimension(:,:) :: rt_patches
  integer :: iostat, n_omega, n
  real, dimension(:), allocatable:: mu, w_mu, mux, muy, muz
  integer, save :: n_nu=1, n_mu=1, n_phi=4
  integer :: i, dir, i_phi, i_mu, p, maxiter
  namelist /experiment_params/ n_nu, n_mu, n_phi, maxiter
  !-----------------------------------------------------------------------------
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, experiment_params)
    end if
    call io%namelist_warning ('experiment_params', iostat)
    self%maxiter = maxiter
    first_time = .false.
  endif

  !-------------------------------------------------------------------------
  ! Assign angles and integration weights
  !-------------------------------------------------------------------------
  n_omega = 2 * n_phi * n_mu ! only for comparison with 'set_a2'
  allocate (mu  (n_omega))
  allocate (w_mu(n_omega))
  !call radau%init (self%n_mu, mu, w_mu)
  mu = 0.57735026 ! only for comparison with 'set_a2'
  mux = [-0.57735026,  0.57735026,  0.57735026, -0.57735026]
  muy = [ 0.57735026,  0.57735026, -0.57735026, -0.57735026]
  w_mu = 1

  allocate(atom_task)
  allocate(atom_task%opac(0:n_omega, n_nu))
  allocate(atom_task%rate(n_omega, n_nu))
  self%atom_task => atom_task

  call self%link_sub_task(atom_task)
  call atom_task%pre_init()
  call atom_task%init()
  call atom_task%append_to_task_list(self%task_list)
  self%sp => atom_task%sp

  allocate(back_opac(n_nu))
  self%back_opac => back_opac

  allocate(sum_opac(n_omega, n_nu))
  self%opac => sum_opac

  allocate(rt_patches(0:n_omega, n_nu))
  self%rt_patches => rt_patches

  do n=1, n_nu
      back_opac(n)%kind = 'back_opac'
      call self%link_and_init(back_opac(n), n, n_nu)

      atom_task%opac(0,n)%kind = 'cont_opac'
      call self%link_and_init(atom_task%opac(0,n), n, n_nu)

      rt_patches(0,n)%kind = 'mean_rad'
      call self%link_and_init(rt_patches(0,n), n, n_nu)

  end do

  do p=1,4 ! do this for line_opac, angle_rad, sum_opac, line_rate
    i = 1
    do i_phi=1, n_phi
      do i_mu=1, n_mu
        do dir=-1,1,2
          do n=1, n_nu
            if (p==1) then
              patch => atom_task%opac(i,n)
              patch%kind = 'line_opac'
            elseif (p==2) then
              patch => atom_task%rate(i,n)
              patch%kind = 'rad_rate'
            elseif (p==3) then
              patch => sum_opac(i,n)
              patch%kind = 'sum_opac'
            elseif (p==4) then
              patch => rt_patches(i,n)
              patch%kind = 'angle_rad'
            endif

            patch%i_omega = i
            patch%i_phi = i_phi
            patch%i_mu = i_mu
            patch%dir = dir

            patch%wmu = w_mu(i) * math%pi4 / (2*n_phi)

            patch%vec(1) = mux(i_phi)
            patch%vec(2) = muy(i_phi)
            patch%vec(3) = dir*mu(i)

            patch%mu  = patch%vec(3) ! only for comparison with 'set_a2'
            patch%phi = patch%vec(1) ! only for comparison with 'set_a2'

            call self%link_and_init(patch, n, n_nu)
          end do
          i = i+1
        end do
      end do
    end do
  end do
END SUBROUTINE init_sub_patches

!===============================================================================
SUBROUTINE link_and_init (self, m3d, n, n_nu)
  class(experiment_t), target :: self
  class(m3d_t) :: m3d
  integer, intent(in) :: n    ! index  of frequency bin
  integer, intent(in) :: n_nu ! number of frequency bins
  !-----------------------------------------------------------------------------
  call self%link_sub_task(m3d)

  select type(m3d)
    class is (opacity_t)
      select case(m3d%kind)
      case ('back_opac')
        m3d%composition = self%composition
        call m3d%composition%reassign_pointer()
        m3d%mean_rad => self%rt_patches(0,n)
      case ('cont_opac')
        m3d%atom => self%atom_task%atom
      case ('line_opac')
        m3d%atom => self%atom_task%atom
        m3d%vproj = m3d%vec(1)*self%vx + m3d%vec(2)*self%vy + m3d%vec(3)*self%vz
      case ('sum_opac')
        m3d%back_opac => self%back_opac(n)
        m3d%cont_opac => self%atom_task%opac(0,n)
        m3d%line_opac => self%atom_task%opac(m3d%i_omega,n)
      end select

    class is (rate_patch_t)
      m3d%atom => self%atom_task%atom
      m3d%rad  => self%rt_patches(m3d%i_omega,n)
      m3d%cont_opac => self%atom_task%opac(0,n)
      m3d%line_opac => self%atom_task%opac(m3d%i_omega,n)
      m3d%vproj = m3d%vec(1)*self%vx + m3d%vec(2)*self%vy + m3d%vec(3)*self%vz

    class is (radiation_t)
      select case(m3d%kind)
      case ('mean_rad')
        m3d%sheep => self%rt_patches(1:, n)
      case ('angle_rad')
        m3d%opac => self%opac(m3d%i_omega,n)
      end select
  end select

  call m3d%pre_init()
  call m3d%set_freq_indices(n, self%sp%nnu, n_nu)
  call m3d%init()
  call m3d%append_to_task_list()

  m3d%boundaries%topbot = self%boundaries%topbot

END SUBROUTINE link_and_init

!===============================================================================
SUBROUTINE link_sub_task (self, task)
  class(experiment_t), target :: self
  class(task_t) :: task
  !-----------------------------------------------------------------------------
  select type (task)
  class is (m3d_t)
    task%patch => self
    task%temp  => self%temp
    task%rho   => self%rho 
    task%ne    => self%ne
    task%sp    => self%sp
    task%maxiter = self%maxiter
  class is (atom_task_t)
    task%patch => self
    task%temp  => self%temp
    task%rho   => self%rho 
    task%ne    => self%ne
    task%nh =>  self%nh
    task%composition => self%composition
  end select
END SUBROUTINE link_sub_task

!===============================================================================
SUBROUTINE freeze_mhd (self)
  class(experiment_t), target :: self
  integer :: i
  !-----------------------------------------------------------------------------
  call self%set(bits%frozen)
  self%time = io%end_time

  !do i=1, size(self%back_opac)
  !  !call self%back_opac(i)%set(bits%frozen)
  !  self%back_opac(i)%time = huge(1.0)
  !end do
END SUBROUTINE

END MODULE
