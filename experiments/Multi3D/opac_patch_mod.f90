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
  USE atom_mod, only: atom_type
  USE m3d_patch_mod, only: maxiter, m3d_patch_t, save_back_opac
  USE ray_patch_mod, only: ray_patch_t
  USE composition_mod, only: composition_type
  implicit none
  private

  !-----------------------------------------------------------------------------
  type, public, extends(ray_patch_t) :: opac_patch_t
  integer :: i_abs=3, i_sca=4
  class(atom_type), pointer :: atom => null()
  class(composition_type), pointer :: composition
  real(4), dimension(:,:,:,:), pointer, contiguous :: wphi
  class(m3d_patch_t), pointer :: rad => null()
  class(m3d_patch_t), pointer :: mean_rad => null()
  contains
    procedure :: pre_init => opacity_pre_init
    procedure :: update => opacity_update
    procedure :: nbor_relations => opacity_nbor_relations
    procedure :: back_opac_update
    procedure :: cont_opac_update
    !procedure :: line_opac_update
    !procedure :: net_opac_update
  end type opac_patch_t

  !-----------------------------------------------------------------------------
  integer, save:: verbose=0                           ! verbosity
    
contains

!===============================================================================
SUBROUTINE opacity_pre_init (self)
  class(opac_patch_t), target:: self
  integer :: nw, fs, fe
  ! ----------------------------------------------------------------------------
  call trace%begin('opac_patch_t%opacity_pre_init')
  call self%ray_patch_t%pre_init()

  self%ng = 0

  if (self%kind == 'cont_opac' .or. self%kind == 'line_opac') then
    if (self%kind == 'cont_opac') then
      fs = max(minval(self%atom%cont%ired),  self%fs)
      fe = min(maxval(self%atom%cont%iblue), self%fe)
    else
      fs = max(minval(self%atom%line%ired),  self%fs)
      fe = min(maxval(self%atom%line%iblue), self%fe)
    endif

    if (fs >= self%fe .or. fe <= self%fs) then
      ! do not add to task list
      call trace%end()
      return
    else
      self%fs = fs
      self%fe = fe
    endif
  endif

  if (self%kind == 'back_opac' .and. save_back_opac) then
    nw = 4
  else
    nw = 2
  endif

  call self%set_nv_and_nw (self%fe-self%fs+1,nw)
  call self%post_pre_init()
  call trace%end()
END SUBROUTINE opacity_pre_init

!===============================================================================
SUBROUTINE opacity_update (self)
  class(opac_patch_t) :: self
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
  !case ('line_opac')
  !  call self%line_opac_update(chi, eta)
  !case ('net_opac')
  !  call self%net_opac_update()
  end select

  call self%post_update()
END SUBROUTINE opacity_update

!===============================================================================
SUBROUTINE back_opac_update (self, chi, eta)
  class(opac_patch_t), target :: self
  real(4), dimension(:,:,:,:), intent(inout) :: chi, eta
  real(8), dimension(:), pointer :: nu
  real(4), dimension(:,:,:,:), pointer, contiguous :: Jold, abslte, sumsca
  character(len=64) :: text
  integer, save:: itimer=0
  integer :: m(6)
  !-----------------------------------------------------------------------------
  call trace%begin('opac_patch_t%back_opac_update', itimer=itimer)
  nu   => self%sp%nu(self%fs:self%fe)

  if (save_back_opac) then
    abslte => self%mem(:,:,:,:,self%it,self%i_abs)
    sumsca => self%mem(:,:,:,:,self%it,self%i_sca)
  else
    m = shape(self%mem)
    call os%allocate(abslte, m(1:4), label='abslte')
    call os%allocate(sumsca, m(1:4), label='sumsca')
  endif

  if (self%iter == 0 .or. .not. save_back_opac) then
    call self%composition%absorption(nu, chi, abslte, sumsca, include_absmet=self%out_ray==.false.)
  endif
  
  if (associated(self%mean_rad)) then
    call self%check_access(self%mean_rad)
    Jold => self%mean_rad%authority_mem(:,:,:,:,self%mean_rad%it,1)
    call self%composition%scattering(nu, abslte, eta, Jold, sumsca)
  else
    call self%composition%scattering(nu, abslte, eta)
  endif

  if (.not. save_back_opac) then
    call os%deallocate(abslte, label='abslte')
    call os%deallocate(sumsca, label='sumsca')
  endif

  call trace%end(itimer)
END SUBROUTINE back_opac_update

!===============================================================================
SUBROUTINE cont_opac_update (self, chi, eta)
  class(opac_patch_t) :: self
  real(4), dimension(:,:,:,:), intent(inout) :: chi, eta
  real(4), dimension(:,:,:), pointer, contiguous :: temp
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('opac_patch_t%cont_opac_update', itimer=itimer)
  temp => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_temp, self%it,1)
  call self%atom%add_continuum_opacities(self%fs, self%fe, self%sp%nu, temp, chi, eta)
  call trace%end(itimer)
END SUBROUTINE cont_opac_update

!===============================================================================
!SUBROUTINE line_opac_update (self, chi, eta)
!  class(opac_patch_t) :: self
!  real(4), dimension(:,:,:,:), intent(inout) :: chi, eta
!  real(4), dimension(:,:,:), pointer, contiguous :: temp
!  integer, save:: itimer=0
!  !-----------------------------------------------------------------------------
!  call trace%begin('opac_patch_t%line_opac_update', itimer=itimer)
!  chi = 0
!  eta = 0
!  temp => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_temp, self%it,1)
!
!  if (self%iter == 0) then
!    call self%comp_vproj(self%vproj)
!    if (.not. self%out_ray) then
!      call self%atom%compute_wphi(self%fs, self%fe, self%sp%nu, temp, self%vproj, self%wphi)
!    endif
!  endif
!  call self%atom%add_line_opacities(self%fs, self%fe, self%sp%nu, temp, self%vproj, chi, eta)
!  call trace%end(itimer)
!END SUBROUTINE line_opac_update

!===============================================================================
!SUBROUTINE net_opac_update (self)
!  class(opac_patch_t) :: self
!  class(nbor_t), pointer:: nbor
!  class(task_t), pointer:: nbtask
!  real(4), dimension(:,:,:,:,:,:), pointer :: mem
!  real(4), dimension(:,:,:,:,:), pointer :: mem1
!  real(4), dimension(:,:,:,:,:), pointer :: chi, eta, Snu
!  integer, save:: itimer=0
!  integer :: iw,it, f, fs, fe, fn
!  !-----------------------------------------------------------------------------
!  call trace%begin('opac_patch_t%net_opac_update', itimer=itimer)
!  mem => self%authority_mem
!  chi => self%authority_mem(:,:,:,:,:,self%i_chi)
!  eta => self%authority_mem(:,:,:,:,:,self%i_eta)
!  Snu => self%authority_mem(:,:,:,:,:,self%i_Snu)
!
!  nbor => self%link%nbor
!  do while (associated(nbor))
!    if (nbor%needed .and. .not. nbor%download) then
!      nbtask => nbor%task
!      select type (nbtask)
!      class is (opac_patch_t)
!      fs = nbtask%fs - self%fs + 1
!      fe = fs + nbtask%nv - 1
!        do iw=1, self%nw
!        do it=1, self%nt
!        !$omp simd
!        do f=fs,fe
!          fn = f - fs + 1
!          mem(:,:,:,f,it,iw) = mem(:,:,:,f,it,iw) + nbtask%mem(:,:,:,fn,it,iw)
!        end do
!        end do
!        end do
!      end select
!    endif
!    nbor => nbor%next
!  end do
!
!  do it=1, self%nt
!  !$omp simd
!  do f=1, self%nv
!    Snu(:,:,:,f,it) = eta(:,:,:,f,it) / chi(:,:,:,f,it)
!  end do
!  end do
!  
!  mem1 => self%mem(:,:,:,:,self%it,:)
!  call self%boundaries%mirror(mem1)
!  call trace%end(itimer)
!END SUBROUTINE net_opac_update

!===============================================================================
SUBROUTINE opacity_nbor_relations(self, another, needed, needs_me, download)
  class(opac_patch_t) :: self
  class(task_t), pointer :: another        ! another task in the list
  logical :: needed, needs_me, download    ! Is the other task: needed, 
                                           ! does it need me, should I download?
  class(task_t), pointer :: other
  logical :: default
  ! ----------------------------------------------------------------------------
  call trace%begin ('opac_patch_t%nbor_relations')
  call self%ray_patch_t%nbor_relations(another, needed, needs_me, download)

  if (needed .or. needs_me) then
    select type(another)
    class is(opac_patch_t)
      select case (self%kind)
      case ('line_opac', 'cont_opac', 'back_opac')
        if (another%kind/=self%kind) then
          needed=.false.
          needs_me=.false.
        endif
      end select
    class is (m3d_patch_t)
      if (another%id == self%m3d_parent%id) then
        needs_me = .false.
      elseif (another%kind == 'net_opac') then
        select case(self%kind)
          case ('line_opac', 'cont_opac', 'back_opac')
            needed=.false.
        end select
      else
        if(self%kind=='back_opac' .and. another%kind=='angle_rad' .and. .not. self%out_ray) then
          needed = .false.
        elseif(self%kind=='back_opac' .and. another%kind=='mean_rad') then
          needs_me = .false.
          self%mean_rad => another
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
