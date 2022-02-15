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
  USE spectrum_mod, only:spectrum_type
  USE m3d_patch_mod, only: m3d_t, halo_m3d_t, pcheck
  USE opac_patch_mod, only: opacity_t
  USE rad_patch_mod, only: radiation_t
  USE atom_patch_mod, only: atom_patch_t
  USE rate_patch_mod, only: rate_patch_t
  USE list_mod
  USE link_mod
  USE task_mod
  USE math_mod
  USE timer_mod
  USE ray_mod, only: set_carlson_quadratures
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
    type(composition_type) :: composition
    class(spectrum_type), pointer :: sp
    class(opacity_t), pointer, dimension(:) :: back_opac
    class(opacity_t), pointer, dimension(:,:) :: opac
    class(radiation_t), pointer, dimension(:) :: out_rad
    class(radiation_t), pointer, dimension(:,:) :: rt_patches
    class(atom_patch_t), pointer :: atom_patch
    real(8), allocatable, dimension(:,:,:) :: ne, temp, rho, vx, vy, vz
    real(8), allocatable, dimension(:,:,:,:) :: nh
    integer :: idx_pe=6, n_mu, n_phi, n_nu, iter=0, maxiter, ilambd
    logical :: nh_missing=.true.
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: update
    procedure:: post_update
    procedure:: init_level
    procedure:: bndries
    procedure:: pre_init_composition
    procedure:: init_composition
    procedure:: append_to
    procedure:: set_t_offset
    procedure:: set_maxiter
    procedure:: link_sub_task
    procedure:: link_and_init
    procedure:: init_sub_patches
    procedure:: nbor_relations => experiment_nbor_relations
  end type
  
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE pre_init (self)
  class(experiment_t), target :: self
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_mod::pre_init', itimer=itimer)
  call self%gpatch_t%pre_init ()
  !call self%pre_init_atom()
  call self%pre_init_composition()
  call trace%end (itimer=itimer)
END SUBROUTINE pre_init

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer, save:: itimer=0
  real(8), save:: timestep=1.d0
  !----------------------------------------------------------------------------
  call trace%begin('experiment_mod::init', itimer=itimer)
!  call self%solver_t%init
  call self%gpatch_t%init                                           ! keep

  allocate(self%nh(self%n(1), self%n(2), self%n(3), 6))

  self%time     = -1.0
  self%dtime    = timestep
  self%out_next = huge(1.0)

  self%kind = 'mhd_patch'

  call trace%end (itimer=itimer)
END SUBROUTINE init

!===============================================================================
!> Setup initial states in the patches
!===============================================================================
SUBROUTINE bndries (self)
  class(experiment_t):: self
  real(var_k), dimension(:,:,:,:,:), pointer :: mem
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%bndries', itimer=itimer)
  mem => self%mem(:,:,:,:,self%it,:)
  if (any(self%boundaries%topbot)) call self%boundaries%boundary(mem)
  call trace%end (itimer)
END SUBROUTINE bndries

!===============================================================================
!> This update only concerns the MHD task; the RT tasks are called directly
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  integer :: it_l, in_l, in_lte_l
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update', itimer=itimer)
  !call self%solver_t%update

  if (self%iter == 0) then
    call self%init_composition()
  endif

  call self%post_update()
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
SUBROUTINE post_update (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  if (self%verbose > 1) then
      write(*,'(F5.1, A12, I4, 12(ES15.7))') self%time, trim(self%kind), self%id, &
      self%mem(self%li(1)+pcheck(1)-1,self%li(2)+pcheck(2)-1,self%li(3)+pcheck(3)-1,:,1,1)
  endif

  if (self%iter == 0) then
    self%finished = .true.
    self%time = io%end_time
    call self%set(bits%frozen)
  endif
  
  self%iter = self%iter + 1
END SUBROUTINE post_update

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
  if (self%nh_missing) then
    call self%composition%get_nhyd(self%nh)
  endif
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
  class(atom_patch_t), pointer :: atom_patch
  class(opacity_t), pointer, dimension(:) :: back_opac
  class(opacity_t), pointer, dimension(:,:) :: sum_opac
  class(radiation_t), pointer, dimension(:) :: out_rad
  class(radiation_t), pointer, dimension(:,:) :: rt_patches
  integer :: iostat, n_omega, n, n_out, i_out, n_mu
  real(8), dimension(:), allocatable:: wmu, long_wmu
  real(8), dimension(:,:), allocatable:: vec, long_vec
  integer, save :: n_nu=1, verbose=0, maxiter=1, ilambd=3
  character(len=32), save :: quad_scheme='set_a4'
  character(len=32), save :: long_scheme='set_a2'
  integer :: i, dir, p
  namelist /experiment_params/ n_nu, maxiter, verbose, ilambd, quad_scheme, long_scheme
  !-----------------------------------------------------------------------------
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, experiment_params)
    end if
    call io%namelist_warning ('experiment_params', iostat)
    write (io%namelist, experiment_params)
    first_time = .false.
  endif

  self%ilambd  = ilambd
  self%maxiter = maxiter
  self%verbose = verbose

  !-------------------------------------------------------------------------
  ! Assign angles and integration weights
  !-------------------------------------------------------------------------
  call set_carlson_quadratures (quad_scheme, vec, wmu, up_and_down=.true.)
  call set_carlson_quadratures (long_scheme, long_vec, long_wmu, up_and_down=.false.)
  n_omega = size(vec, 2)
  n_out = size(long_vec, 2)
  n_mu = n_omega + n_out

  allocate(atom_patch)
  allocate(atom_patch%opac(0:n_mu, n_nu))
  allocate(atom_patch%rate(n_omega, n_nu))
  self%atom_patch => atom_patch
  atom_patch%n_omega = n_omega
  atom_patch%n_nu = n_nu

  call self%link_sub_task(atom_patch)
  call atom_patch%pre_init()
  call self%set_t_offset(atom_patch)
  call self%set_maxiter(atom_patch)
  call atom_patch%init()
  call atom_patch%append_to_task_list(self%task_list)
  self%sp => atom_patch%sp

  allocate(back_opac(n_nu))
  self%back_opac => back_opac

  allocate(sum_opac(n_mu, n_nu))
  self%opac => sum_opac

  allocate(rt_patches(0:n_mu, n_nu))
  self%rt_patches => rt_patches

  if (self%boundaries%topbot(1)) then
    allocate(out_rad(n_out))
    self%out_rad => out_rad
  endif

  do n=1, n_nu
      rt_patches(0,n)%kind = 'mean_rad'
      call self%link_and_init(rt_patches(0,n), n, n_nu)

      atom_patch%opac(0,n)%kind = 'cont_opac'
      call self%link_and_init(atom_patch%opac(0,n), n, n_nu)

      back_opac(n)%kind = 'back_opac'
      call self%link_and_init(back_opac(n), n, n_nu)
  end do

  do p=1,4 ! do this for line_opac, angle_rad, sum_opac, line_rate
    do i=1, n_mu
      do n=1, n_nu
        if (p==4) then
          patch => atom_patch%opac(i,n)
          patch%kind = 'line_opac'
        elseif (p==3) then
          patch => sum_opac(i,n)
          patch%kind = 'sum_opac'
        elseif (p==2) then
          patch => rt_patches(i,n)
          patch%kind = 'angle_rad'
        elseif (p==1) then
          if (i > n_omega) cycle
          patch => atom_patch%rate(i,n)
          patch%kind = 'rad_rate'
        endif

        patch%i_mu = i
        patch%i_nu = n

        if (i > n_omega) then
          patch%out_ray = .true.
          patch%i_out = i - n_omega 
          patch%vec = long_vec(:,patch%i_out)
          patch%wmu = long_wmu(patch%i_out)
          patch%dir = 1
        else
          patch%vec = vec(:,i)
          patch%wmu = wmu(i)
          patch%dir = merge(1, -1, vec(3,i) > 0)
        endif


        call self%link_and_init(patch, n, n_nu)
      end do
    end do
  end do

  if (associated(self%out_rad)) then
    do i=1, n_out
      out_rad(i)%kind = 'out_rad'
      out_rad(i)%i_mu  = i + n_omega
      out_rad(i)%i_out = i
      out_rad(i)%out_ray = .true.
      call self%link_and_init(out_rad(i), 1, 1)
    end do
  endif
END SUBROUTINE init_sub_patches

!===============================================================================
SUBROUTINE set_t_offset (self, task)
  class (experiment_t) :: self
  class(task_t) :: task
  real(8) :: t_offset
  !-----------------------------------------------------------------------------
  select type (task)
  class is (atom_patch_t)
    t_offset = 0.1
  class is (m3d_t)
    select case (task%kind)
    case ('back_opac', 'cont_opac', 'line_opac')
      t_offset = 0.2
    case ('sum_opac')
      t_offset = 0.4
    case ('angle_rad')
      t_offset = 0.6
    case ('mean_rad', 'out_rad')
      t_offset = 0.8
    case ('rad_rate')
      t_offset = 0.9
    end select
    
    if (task%out_ray) then
      t_offset = t_offset + self%ilambd + self%maxiter
    endif
  end select

  task%time = self%time + t_offset * self%dtime
END SUBROUTINE set_t_offset

!===============================================================================
SUBROUTINE set_maxiter (self, task)
  class (experiment_t) :: self
  class(task_t) :: task
  !-----------------------------------------------------------------------------
  select type (task)
  class is (atom_patch_t)
    task%maxiter = self%maxiter
    task%rate%maxiter = task%maxiter - 1
    task%opac(0,:)%maxiter = task%maxiter
    task%opac(1:task%n_omega,:)%maxiter = task%maxiter - 1
    task%opac(task%n_omega+1:,:)%maxiter = 0
  class is (m3d_t)
    select case (task%kind)
    case ('back_opac')
      ! LTE iteration counts as first lambda iteration
      task%maxiter = max(0, self%ilambd - 1)
    case ('mean_rad')
      task%maxiter = self%ilambd - 1
    case ('angle_rad', 'sum_opac')
      if (task%out_ray) then
        task%maxiter = 0
      else
        task%maxiter = self%ilambd + self%maxiter - 1
      endif
    case ('out_rad')
      task%maxiter = 0
    end select
  end select
END SUBROUTINE set_maxiter

!===============================================================================
SUBROUTINE link_sub_task (self, task)
  class(experiment_t), target :: self
  class(task_t) :: task
  !-----------------------------------------------------------------------------
  select type (task)
  class is (m3d_t)
    task%patch => self
    task%temp = self%temp
    task%rho  = self%rho 
    task%ne   = self%ne
    task%sp   = self%sp

  class is (atom_patch_t)
    task%patch => self
    task%temp  => self%temp
    task%rho   => self%rho 
    task%ne    => self%ne
    task%nh =>  self%nh
    task%composition => self%composition
    task%ilambd = self%ilambd
  end select

  task%verbose = self%verbose
END SUBROUTINE link_sub_task

!===============================================================================
SUBROUTINE link_and_init (self, m3d, n, n_nu)
  class(experiment_t), target :: self
  class(m3d_t), target :: m3d
  integer, intent(in) :: n    ! index  of frequency bin
  integer, intent(in) :: n_nu ! number of frequency bins
  integer :: i
  !-----------------------------------------------------------------------------
  call self%link_sub_task(m3d)

  select type(m3d)
    class is (opacity_t)
      select case(m3d%kind)
      case ('back_opac')
        m3d%composition = self%composition
        call m3d%composition%reassign_pointer()
        m3d%mean_rad => self%rt_patches(0, n)
      case ('cont_opac')
        m3d%atom => self%atom_patch%atom
      case ('line_opac')
        m3d%atom => self%atom_patch%atom
        !-----------------------------------------------------------------------------
        ! There appears to be an issue with angles in M3D. The vectors of the final
        ! intensity rays are not in line with the vectors used for the projected
        ! velocities. The x and y coordinates seem to be flipped.
        !-----------------------------------------------------------------------------
        if (m3d%out_ray) then
          m3d%vproj = -m3d%vec(1)*self%vx - m3d%vec(2)*self%vy + m3d%vec(3)*self%vz 
        else
          m3d%vproj = m3d%vec(1)*self%vx + m3d%vec(2)*self%vy + m3d%vec(3)*self%vz 
        endif
      case ('sum_opac')
        m3d%back_opac => self%back_opac(n)
        m3d%cont_opac => self%atom_patch%opac(0,n)
        m3d%line_opac => self%atom_patch%opac(m3d%i_mu,n)
        m3d%rad => self%rt_patches(m3d%i_mu,n)
      end select

    class is (rate_patch_t)
      m3d%atom => self%atom_patch%atom
      m3d%rad  => self%rt_patches(m3d%i_mu,n)
      m3d%cont_opac => self%atom_patch%opac(0,n)
      m3d%line_opac => self%atom_patch%opac(m3d%i_mu,n)
      m3d%vproj = m3d%vec(1)*self%vx + m3d%vec(2)*self%vy + m3d%vec(3)*self%vz

    class is (radiation_t)
      select case(m3d%kind)
      case ('mean_rad')
        m3d%sheep => self%rt_patches(1:, n)
      case ('out_rad')
        m3d%sheep => self%rt_patches(m3d%i_mu, :)
        do i=1, size(m3d%sheep)
          m3d%sheep(i)%out_rad => m3d
        end do
        m3d%vec = m3d%sheep(1)%vec
      case ('angle_rad')
        m3d%opac => self%opac(m3d%i_mu,n)
        m3d%mean_rad => self%rt_patches(0, n)
      end select
  end select

  call m3d%pre_init()
  call self%set_t_offset(m3d)
  call self%set_maxiter(m3d)
  call m3d%set_freq_indices(n, self%sp%nnu, n_nu)
  call m3d%init()
  call m3d%append_to_task_list()
END SUBROUTINE link_and_init

!===============================================================================
SUBROUTINE experiment_nbor_relations(self, another, needed, needs_me, download)
  class(experiment_t) :: self
  class(task_t), pointer :: another         ! another task in the list
  logical :: needed, needs_me, download    ! Is the other task: needed, 
                                           !    does it need me, should I download?
  
  call trace%begin ('atom_patch_t%nbor_relations')
  needed = .false.; needs_me = .false.; download = .false.

  if (another%kind == self%kind) then
    ! First check the main relations, based on location and level of refinement
    call self%gpatch_t%nbor_relations(another, needed, needs_me, download)
  endif

  select type(another)
  class is (solver_t)
    ! use default
  class is (atom_patch_t)
    if (another%patch%id == self%id) then
      needed = .false.; needs_me = .true.;  download = .false.
    endif
  class is (opacity_t)
    if (another%kind == 'back_opac' .and. another%patch%id == self%id) then
      needed = .false.; needs_me = .true.;  download = .false.
    endif
  class default
    needed = .false.; needs_me = .false.; download = .false.
  end select

  call trace%end()
END SUBROUTINE experiment_nbor_relations

END MODULE
