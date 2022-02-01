module opac_patch_mod
  !USE rt_integral_mod
  !USE radau_mod
  !USE boundaries_mod
  !USE eos_mod
  !USE list_mod
  USE halo_mod
  USE gpatch_mod
  USE solver_mod
  USE link_mod
  USE task_mod
  !USE index_mod
  !USE mesh_mod
  USE io_mod
  !USE aux_mod
  USE timer_mod
  !USE scaling_mod
  USE trace_mod
  !USE mpi_mod
  !USE omp_mod
  !USE omp_timer_mod
  !USE math_mod
  USE bits_mod
  !USE units_mod
  !USE scratch_mod
  USE atom_var, only: atom_type
  USE composition_mod, only:composition_type
  USE m3d_patch_mod, only:m3d_t
  implicit none
  private

  !-----------------------------------------------------------------------------
  type, public, extends(m3d_t) :: opacity_t
  integer :: i_chi=1, i_eta=2, i_snu=2
  class(atom_type), pointer :: atom => null()
  type(composition_type) :: composition
  real(8), dimension(:,:,:), allocatable :: vproj
  real(8), dimension(:,:,:,:), allocatable :: nh, wphi
  class(opacity_t), pointer :: back_opac => null()
  class(opacity_t), pointer :: cont_opac => null()
  class(opacity_t), pointer :: line_opac => null()
  class(m3d_t), pointer :: rad => null()
  class(m3d_t), pointer :: mean_rad => null()
  contains
    procedure :: pre_init => opacity_pre_init
    procedure :: init => opacity_init
    procedure :: update => opacity_update
    procedure :: nbor_relations => opacity_nbor_relations
    procedure :: init_composition
    procedure :: back_opac_update
    procedure :: cont_opac_update
    procedure :: line_opac_update
    procedure :: sum_opac_update
  end type opacity_t




  !-----------------------------------------------------------------------------
  integer, save:: verbose=0                           ! verbosity
    
contains

!===============================================================================
SUBROUTINE opacity_pre_init (self)
  class(opacity_t), target:: self
  ! ----------------------------------------------------------------------------
  call self%m3d_t%pre_init()

  self%nw = 2
  self%mv = self%nv
  !self%nv is assigned in experiment_mod : init_sub_patches
  self%nt = 1
  self%it = 1
  self%new = self%it

  if (self%kind == 'sum_opac') then
    self%ng = [3,3,1]
  endif
  
  ! ----------------------------------------------------------------------------
  ! Allocate and set unsigned and pervolume. Opacity, source function,
  ! intensity and populations are all positive definite
  ! Although quantities are generally per volume, but density is not stored
  ! locally, thus this flag must not be set.
  ! ----------------------------------------------------------------------------
  allocate (self%unsigned(self%nv), self%pervolume(self%nv))
  self%pervolume = .false.
  self%unsigned  = .true.
END SUBROUTINE opacity_pre_init

!===============================================================================
SUBROUTINE init_composition (self)
  class(opacity_t) :: self
  !-----------------------------------------------------------------------------
  ! loading temp, ne and rho in to composition objekt
  ! ion and molecule equilibrium are solved immediately
  call self%composition%load_state_variables(self%temp, self%ne, self%rho)

  allocate(self%nh(self%n(1), self%n(2), self%n(3), 6))
  call self%composition%get_nhyd(self%nh)
END SUBROUTINE init_composition

!===============================================================================
SUBROUTINE opacity_init (self)
  class(opacity_t) :: self
  ! ----------------------------------------------------------------------------
  call self%m3d_t%init()

  if (self%kind == 'sum_opac') then
    self%t_offset = 0.4d0
  else
    self%t_offset = 0.2d0
  endif

  self%time = self%patch%time + self%t_offset*self%dtime
END SUBROUTINE opacity_init

!===============================================================================
SUBROUTINE opacity_update (self)
  class(opacity_t) :: self
  real(4), dimension(:,:,:,:), pointer :: chi, eta
  ! ----------------------------------------------------------------------------
  chi => self%authority_mem(:,:,:,:, self%it, self%i_chi)
  eta => self%authority_mem(:,:,:,:, self%it, self%i_eta)

  select case (self%kind)
  case ('back_opac')
    call self%back_opac_update(chi, eta)
  case ('cont_opac')
    call self%cont_opac_update(chi, eta)
  case ('line_opac')
    call self%line_opac_update(chi, eta)
  case ('sum_opac')
    call self%sum_opac_update()
  end select

  call self%post_update()
END SUBROUTINE opacity_update

!===============================================================================
SUBROUTINE back_opac_update (self, chi, eta)
  class(opacity_t), target :: self
  real(4), dimension(:,:,:,:), intent(inout) :: chi, eta
  real(8), dimension(:), pointer :: nu
  real(4), dimension(:,:,:,:), pointer :: Jold
  character(len=64) :: text
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('opacity_t%back_opac_update', itimer=itimer)
  nu   => self%sp%nu(self%fs:self%fe)
  Jold => self%mean_rad%authority_mem(:,:,:,:,self%mean_rad%it,1)

  if (self%iter == 0) then
    chi = 0
    eta = 0
    call self%init_composition()
    call self%composition%absorption(nu, chi)
  endif

  call self%composition%scattering(nu, Jold, eta)

  call trace%end(itimer)
END SUBROUTINE back_opac_update

!===============================================================================
SUBROUTINE cont_opac_update (self, chi, eta)
  class(opacity_t) :: self
  real(4), dimension(:,:,:,:), intent(inout) :: chi, eta
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('opacity_t%cont_opac_update', itimer=itimer)
  chi = 0
  eta = 0
  call self%atom%add_continuum_opacities(self%fs, self%fe, self%sp%nu, self%temp, chi, eta)
  call trace%end(itimer)
END SUBROUTINE cont_opac_update

!===============================================================================
SUBROUTINE line_opac_update (self, chi, eta)
  class(opacity_t) :: self
  real(4), dimension(:,:,:,:), intent(inout) :: chi, eta
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('opacity_t%line_opac_update', itimer=itimer)
  chi = 0
  eta = 0
  if (self%iter == 0) then
    call self%atom%compute_wphi(self%fs, self%fe, self%sp%nu, self%temp, self%vproj, self%wphi)
  endif
  call self%atom%add_line_opacities(self%fs, self%fe, self%sp%nu, self%temp, self%vproj, chi, eta)
  call trace%end(itimer)
END SUBROUTINE line_opac_update

!===============================================================================
SUBROUTINE sum_opac_update (self)
  class(opacity_t) :: self
  real(4), dimension(:,:,:,:,:), pointer :: mem
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('opacity_t%sum_opac_update', itimer=itimer)
  mem => self%authority_mem(:,:,:,:,self%it,:)
  mem = self%back_opac%mem(:,:,:,:,self%back_opac%it,:)
  mem = self%cont_opac%mem(:,:,:,:,self%cont_opac%it,:) + mem
  mem = self%line_opac%mem(:,:,:,:,self%line_opac%it,:) + mem
  mem(:,:,:,:,self%i_snu) = mem(:,:,:,:,self%i_eta) / mem(:,:,:,:,self%i_chi)
  mem => self%mem(:,:,:,:,self%it,:)
  call self%boundaries%mirror(mem)
  call trace%end(itimer)
END SUBROUTINE sum_opac_update

!===============================================================================
SUBROUTINE opacity_nbor_relations(self, another, needed, needs_me, download)
  class(opacity_t) :: self
  class(task_t), pointer :: another        ! another task in the list
  logical :: needed, needs_me, download    ! Is the other task: needed, 
                                           ! does it need me, should I download?
  class(task_t), pointer :: other
  logical :: default
  ! ----------------------------------------------------------------------------
  call trace%begin ('opacity_t%nbor_relations')

  ! If this RT task has no RT on, it should have no nbor relations
  if (.not.(self%on)) then
    needed = .false.; needs_me = .false.; download = .false.
    call trace%end()
    return
  end if

  needed = .false.; needs_me = .false.; download = .false.

  if (self%kind == 'sum_opac') then
    if (another%id == self%back_opac%id) then
      needed = .true.; needs_me = .false.; download = .false.

    elseif (another%id == self%cont_opac%id) then
      needed = .true.; needs_me = .false.; download = .false.

    elseif (another%id == self%line_opac%id) then
      needed = .true.; needs_me = .false.; download = .false.

    elseif (another%id == self%rad%id) then
      needed = .false.; needs_me = .true.; download = .false.
    else
      select type (another)
      class is (opacity_t)
        if (another%i_nu == self%i_nu) then
          if (another%i_omega == self%i_omega) then
            if (another%kind == self%kind) then
              if (self%is_upstream(another)) then
                needed = .true.; needs_me = .false.; download = .true.
              endif
            endif
          endif
        endif
      end select
    endif
  else
    select type (another)
    class is (opacity_t)
      if (another%i_nu == self%i_nu) then
        if (another%kind == 'sum_opac') then
          if (another%patch%id == self%patch%id) then
            needed = .false.; needs_me = .true.; download = .false.
          endif
        endif
      endif
    end select
  endif

  select case (self%kind)
  case ('back_opac')
    if (another%id == self%patch%id) then
      needed = .true. ; needs_me = .false.; download = .false.
    elseif (another%id == self%mean_rad%id) then
      needed = .true. ; needs_me = .false.; download = .false.
    endif
  case ('line_opac', 'cont_opac')
    if (all(another%position == self%position)) then
      if (another%kind == trim(self%atom%id) // '_patch') then
        needed = .true. ; needs_me = .false.; download = .false.
      endif
    endif
  end select

  call trace%end()
END SUBROUTINE opacity_nbor_relations
  
  
      
  end module opac_patch_mod
