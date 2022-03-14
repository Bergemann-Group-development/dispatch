module opac_patch_mod
  !USE rt_integral_mod
  !USE radau_mod
  !USE boundaries_mod
  !USE eos_mod
  !USE list_mod
  USE os_mod
  USE halo_mod
  USE gpatch_mod
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
  USE os_mod
  !USE units_mod
  !USE scratch_mod
  USE atom_var, only: atom_type
  USE m3d_patch_mod, only:maxiter, m3d_patch_t
  USE ray_patch_mod, only:ray_patch_t
  USE composition_mod, only:composition_type
  implicit none
  private

  !-----------------------------------------------------------------------------
  type, public, extends(ray_patch_t) :: opacity_t
  integer :: i_abs=3, i_sca=4
  class(atom_type), pointer :: atom => null()
  class(composition_type), pointer :: composition
  real(8), dimension(:,:,:,:), pointer, contiguous :: wphi
  class(m3d_patch_t), pointer :: rad => null()
  class(m3d_patch_t), pointer :: mean_rad => null()
  contains
    procedure :: pre_init => opacity_pre_init
    procedure :: update => opacity_update
    procedure :: nbor_relations => opacity_nbor_relations
    procedure :: back_opac_update
    procedure :: cont_opac_update
    procedure :: line_opac_update
    procedure :: net_opac_update
  end type opacity_t

  !-----------------------------------------------------------------------------
  integer, save:: verbose=0                           ! verbosity
    
contains

!===============================================================================
SUBROUTINE opacity_pre_init (self)
  class(opacity_t), target:: self
  ! ----------------------------------------------------------------------------
  call trace%begin('opacity_t%opacity_pre_init')
  call self%ray_patch_t%pre_init()
  call self%append_to_task_list()

  if (self%kind == 'back_opac') then
    call self%set_nv_and_nw (self%fe-self%fs+1,4)
  else
    call self%set_nv_and_nw (self%fe-self%fs+1,2)
  endif

  if (self%kind == 'net_opac') then
    if (self%out_ray) then
      self%ng = [0,0,1]
    else
      if (self%ray%long) then
        self%ng = [3,3,1]
      else
        self%ng = [2,2,1]
      endif
    endif
  endif
  call trace%end()
END SUBROUTINE opacity_pre_init

!===============================================================================
SUBROUTINE opacity_update (self)
  class(opacity_t) :: self
  real(4), dimension(:,:,:,:), pointer :: chi, eta
  ! ----------------------------------------------------------------------------
  call self%pre_update()
  if (self%finished) then
    return
  endif
  
  chi => self%authority_mem(:,:,:,:, self%it, self%i_chi)
  eta => self%authority_mem(:,:,:,:, self%it, self%i_eta)

  if (self%kind /= 'back_opac') then
    chi = 0
  endif
  eta = 0

  select case (self%kind)
  case ('back_opac')
    call self%back_opac_update(chi, eta)
  case ('cont_opac')
    call self%cont_opac_update(chi, eta)
  case ('line_opac')
    call self%line_opac_update(chi, eta)
  case ('net_opac')
    call self%net_opac_update()
  end select

  call self%post_update()
END SUBROUTINE opacity_update

!===============================================================================
SUBROUTINE back_opac_update (self, chi, eta)
  class(opacity_t), target :: self
  real(4), dimension(:,:,:,:), intent(inout) :: chi, eta
  real(8), dimension(:), pointer :: nu
  real(4), dimension(:,:,:,:), pointer, contiguous :: Jold, abslte, sumsca
  character(len=64) :: text
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('opacity_t%back_opac_update', itimer=itimer)
  nu   => self%sp%nu(self%fs:self%fe)
  abslte => self%mem(:,:,:,:,self%it,self%i_abs)
  sumsca => self%mem(:,:,:,:,self%it,self%i_sca)

  if (self%iter == 0) then
    chi = 0
    !call self%init_composition()
    call self%composition%absorption(nu, chi, abslte, sumsca)
  endif

  if (self%iter > 0) then
    Jold => self%mean_rad%authority_mem(:,:,:,:,self%mean_rad%it,1)
    call self%composition%scattering(nu, Jold, abslte, sumsca, eta)
  else
    call os%allocate(Jold, shape(self%mem(:,:,:,:,1,1)), zero=.true.)
    call self%composition%scattering(nu, Jold, abslte, sumsca, eta)
    call os%deallocate(Jold)
  endif
  
  if (self%iter == self%maxiter-1) then
    self%time = self%time + (maxiter - self%iter - 1)*self%dtime
  endif

  call trace%end(itimer)
END SUBROUTINE back_opac_update

!===============================================================================
SUBROUTINE cont_opac_update (self, chi, eta)
  class(opacity_t) :: self
  real(4), dimension(:,:,:,:), intent(inout) :: chi, eta
  real(4), dimension(:,:,:), pointer, contiguous :: temp
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('opacity_t%cont_opac_update', itimer=itimer)
  temp => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_temp, self%it,1)
  call self%atom%add_continuum_opacities(self%fs, self%fe, self%sp%nu, temp, chi, eta)
  call trace%end(itimer)
END SUBROUTINE cont_opac_update

!===============================================================================
SUBROUTINE line_opac_update (self, chi, eta)
  class(opacity_t) :: self
  real(4), dimension(:,:,:,:), intent(inout) :: chi, eta
  real(4), dimension(:,:,:), pointer, contiguous :: temp
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('opacity_t%line_opac_update', itimer=itimer)
  chi = 0
  eta = 0
  temp => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_temp, self%it,1)

  if (self%iter == 0) then
    call self%comp_vproj(self%vproj)
    call self%atom%compute_wphi(self%fs, self%fe, self%sp%nu, temp, self%vproj, self%wphi)
  endif
  call self%atom%add_line_opacities(self%fs, self%fe, self%sp%nu, temp, self%vproj, chi, eta)
  call trace%end(itimer)
END SUBROUTINE line_opac_update

!===============================================================================
SUBROUTINE net_opac_update (self)
  class(opacity_t) :: self
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: nbtask
  real(4), dimension(:,:,:,:,:,:), pointer :: mem
  real(4), dimension(:,:,:,:,:), pointer :: mem1
  real(4), dimension(:,:,:,:,:), pointer :: chi, eta, Snu
  integer, save:: itimer=0
  integer :: iw,it,f
  !-----------------------------------------------------------------------------
  call trace%begin('opacity_t%net_opac_update', itimer=itimer)
  mem => self%authority_mem
  chi => self%authority_mem(:,:,:,:,:,self%i_chi)
  eta => self%authority_mem(:,:,:,:,:,self%i_eta)
  Snu => self%authority_mem(:,:,:,:,:,self%i_Snu)

  nbor => self%link%nbor
  do while (associated(nbor))
    if (nbor%needed .and. .not. nbor%download) then
      nbtask => nbor%task
      select type (nbtask)
      class is (opacity_t)
        do iw=1, self%nw
        do it=1, self%nt
        !$omp simd
        do f=1, self%nv
          mem(:,:,:,f,it,iw) = mem(:,:,:,f,it,iw) + nbtask%mem(:,:,:,f,it,iw)
        end do
        end do
        end do
      end select
    endif
    nbor => nbor%next
  end do

  do it=1, self%nt
  !$omp simd
  do f=1, self%nv
    Snu(:,:,:,f,it) = eta(:,:,:,f,it) / chi(:,:,:,f,it)
  end do
  end do
  
  mem1 => self%mem(:,:,:,:,self%it,:)
  call self%boundaries%mirror(mem1)
  call trace%end(itimer)
END SUBROUTINE net_opac_update

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
  call self%ray_patch_t%nbor_relations(another, needed, needs_me, download)

  if (needed .or. needs_me) then
    select type(another)
    class is(opacity_t)
      select case (self%kind)
      case ('net_opac')
        needs_me=.false.
      case ('line_opac', 'cont_opac', 'back_opac')
        if (.not. download) then
          needed=.false.
        endif
        select case (another%kind)
        case ('line_opac', 'cont_opac', 'back_opac')
          needs_me=.false.
        end select
      end select
    class is (m3d_patch_t)
      if (another%id == self%m3d_parent%id) then
        needs_me = .false.
        if (self%kind=='net_opac') then
          needed = .false.
        endif
      else
        if (self%kind=='net_opac' .and. another%kind=='angle_rad') then
          needed = .false.
        elseif(self%kind=='back_opac' .and. another%kind=='mean_rad') then
          needs_me = .false.
          self%mean_rad => another
        elseif(self%m3d_parent%kind=='atom_rad' .and. another%kind=='rate_rad') then
          needed = .false.
        else
          needed = .false.; needs_me = .false.
        endif
      endif
    class default
      needed = .false.; needs_me = .false.
    end select
  endif
  call trace%end()
END SUBROUTINE opacity_nbor_relations
  end module opac_patch_mod
