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
  class(opacity_t), pointer :: back_opac, cont_opac, line_opac
  class(m3d_t), pointer :: mean_rad
  contains
    procedure :: pre_init => opacity_pre_init
    procedure :: init => opacity_init
    procedure :: update => opacity_update
    procedure :: nbor_relations => opacity_nbor_relations
    procedure :: init_composition
  end type opacity_t




  !-----------------------------------------------------------------------------
  integer, save:: verbose=0                           ! verbosity
    
contains

!===============================================================================
SUBROUTINE opacity_pre_init (self)
  class(opacity_t), target:: self
  ! ----------------------------------------------------------------------------
  call self%m3d_t%pre_init()

  if (self%kind == 'sum_opac') then
    self%ng = [4,4,1]
  endif

  self%nv = 2
  self%mv = self%nv
  !self%nw is assigned in experiment_mod : init_opacity_patches
  self%nt = 1
  self%it = 1
  self%new = self%it
  
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
  class(atom_type), pointer :: atom
  class(composition_type), pointer :: composition
  integer :: fs, fe
  real(4), dimension(:,:,:,:), pointer :: chi, eta
  real(4), dimension(:,:,:,:,:), pointer :: mem
  real(4), dimension(:,:,:,:), pointer :: Jold
  ! ----------------------------------------------------------------------------
  !call timer%tic ()

  fs = self%fs
  fe = self%fe
  mem => self%authority_mem

  chi => mem(:,:,:, self%i_chi, :)
  eta => mem(:,:,:, self%i_eta, :)

  chi = 0
  eta = 0

  if (associated(self%atom)) atom => self%atom

  select case (self%kind)
  case ('back_opac')
    if (self%iter == 0) then
      call self%init_composition()
      !call self%set(bits%frozen)
    endif

    if (self%patch%is_set(bits%frozen)) self%time = huge(1.0)

    Jold => self%mean_rad%authority_mem(:,:,:,1,:)
    call self%composition%absorption(self%sp%nu(fs:fe), Jold, chi, eta)

  case ('cont_opac')
    call atom%add_continuum_opacities(fs, fe, self%sp%nu, self%temp, chi, eta)

  case ('line_opac')
    if (self%iter == 0) then
      call atom%compute_wphi(fs, fe, self%sp%nu, self%temp, self%vproj, self%wphi)
    endif
    call atom%add_line_opacities(fs, fe, self%sp%nu, self%temp, self%vproj, chi, eta)

  case ('sum_opac')
    mem = self%back_opac%mem(:,:,:,:,self%it,:)
    mem = self%cont_opac%mem(:,:,:,:,self%it,:) + mem
    mem = self%line_opac%mem(:,:,:,:,self%it,:) + mem
    mem(:,:,:,self%i_snu,:) = mem(:,:,:,self%i_eta,:) / mem(:,:,:,self%i_chi,:)
    call self%bndries
  end select

  self%iter = self%iter + 1

  !select case (self%dir)
  !case(0, 1)
  !  write(*,'(A10, 2(I4), 3(ES18.7))'), trim(self%kind), self%i_omega, self%i_nu, chi(1,1,1,1), eta(1,1,1,1)
  !case(-1)
  !  write(*,'(A10, 2(I4), 3(ES18.7))'), trim(self%kind), self%i_omega, self%i_nu, chi(1,1,self%n(3),1), eta(1,1,self%n(3),1)
  !end select

  !call timer%toc ('update took ', 1)
END SUBROUTINE opacity_update

!===============================================================================
SUBROUTINE opacity_nbor_relations(self, another, needed, needs_me, download)
  class(opacity_t) :: self
  class(task_t), pointer :: another        ! another task in the list
  logical :: needed, needs_me, download    ! Is the other task: needed, 
                                           ! does it need me, should I download?
  class(task_t), pointer :: other
  logical :: default
  
  call trace%begin ('opacity_t%nbor_relations')

  ! If this RT task has no RT on, it should have no nbor relations
  if (.not.(self%on)) then
    needed = .false.; needs_me = .false.; download = .false.
    call trace%end()
    return
  end if

  ! First check the main relations, based on location and level of refinement
  call self%m3d_t%nbor_relations(another, needed, needs_me, download)

  if (self%kind == 'sum_opac') then
    if (another%id == self%patch%id) then
      needed = .false.; needs_me = .false.; download = .false.

    elseif (another%id == self%back_opac%id) then
      needed = .true.; needs_me = .false.; download = .false.

    elseif (another%id == self%cont_opac%id) then
      needed = .true.; needs_me = .false.; download = .false.

    elseif (another%id == self%line_opac%id) then
      needed = .true.; needs_me = .false.; download = .false.
    else
      select type (another)
      class is (opacity_t)
        if (another%i_nu == self%i_nu .and. another%on) then
          if (another%kind == self%kind) then
            if (another%i_omega == self%i_omega) then
              if (self%is_upstream(another)) then
                needed = .true.; needs_me = .false.; download = .true.
              endif
            endif
          endif
        endif
      end select
    endif
  elseif (self%kind == 'back_opac') then
    if (another%id == self%mean_rad%id) then
      needed = .true.; needs_me = .false.; download = .false.
    endif
  endif

  call trace%end()
END SUBROUTINE opacity_nbor_relations
  
  
      
  end module opac_patch_mod