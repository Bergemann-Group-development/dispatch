MODULE atom_cont_patch_mod
    USE os_mod
    USE halo_mod
    USE gpatch_mod
    USE link_mod
    USE task_mod
    USE io_mod
    USE timer_mod
    USE trace_mod
    USE bits_mod
    USE os_mod
    USE handyfuncs_mod
    USE atom_mod, only: atom_type
    USE m3d_patch_mod, only: m3d_patch_t, comoving_frame
    USE ray_patch_mod, only: ray_patch_t
    implicit none
    private
  
    !-----------------------------------------------------------------------------
    type, public, extends(ray_patch_t) :: atom_cont_patch_t
    contains
      procedure :: pre_init => atom_cont_pre_init
      procedure :: execute
      procedure :: nbor_relations
      procedure :: init_nbors_hash
    end type atom_cont_patch_t
      
  contains
  
  !===============================================================================
  SUBROUTINE atom_cont_pre_init (self)
    class(atom_cont_patch_t), target:: self
    integer :: nw, fs, fe
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin('atom_cont_patch_t%opacity_pre_init', itimer=itimer)
    self%kind = 'atom_cont'
    call self%ray_patch_t%pre_init()
  
    self%ng = 0
    fs = max(minval(self%atom%cont%ired),  self%fs)
    fe = min(maxval(self%atom%cont%iblue), self%fe)
  
    if (fs >= self%fe .or. fe <= self%fs) then
      ! do not add to task list
      call trace%end(itimer)
      return
    else
      self%fs = fs
      self%fe = fe
    endif
  
    call self%set_nv_and_nw (self%fe-self%fs+1, 2)
    call self%post_pre_init()
    call trace%end(itimer)
  END SUBROUTINE atom_cont_pre_init
  
  !===============================================================================
  SUBROUTINE execute (self)
    class(atom_cont_patch_t) :: self
    real(4), dimension(:,:,:,:), pointer, contiguous :: chi, eta
    real(4), dimension(:,:,:), pointer, contiguous :: temp
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin (trim(self%kind)//'%update', itimer=itimer)
    chi => self%authority_mem(:,:,:,:, self%it, self%i_chi)
    eta => self%authority_mem(:,:,:,:, self%it, self%i_eta)
    temp => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_temp, self%it,1)
    chi = 0
    eta = 0
    call self%atom%add_continuum_opacities(self%fs, self%fe, self%sp%nu, temp, chi, eta)
    if (self%iter == 1) then
      self%dtime = self%dt_fixed
    endif
    call trace%end(itimer)
  END SUBROUTINE execute
  
  !===============================================================================
  SUBROUTINE nbor_relations(self, another, needed, needs_me, download)
    class(atom_cont_patch_t) :: self
    class(task_t), pointer :: another        ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                             ! does it need me, should I download?
    class(task_t), pointer :: other
    logical :: default
    ! ----------------------------------------------------------------------------
    call trace%begin ('atom_cont_patch_t%nbor_relations')
    call self%ray_patch_t%nbor_relations(another, needed, needs_me, download)
  
    if (needed .or. needs_me) then
      select type(another)
      class is (m3d_patch_t)
        if (another%id == self%m3d_parent%id) then
          needs_me = .false.
        elseif (another%kind == 'net_opac') then
          needed=.false.
        else
          needed = .false.; needs_me = .false.
        endif
      class default
        needed = .false.; needs_me = .false.
      end select
    endif
    call trace%end()
  END SUBROUTINE nbor_relations
  
  !===============================================================================
  !> Fast method to generate nbor lists in cases with very long task lists
  !===============================================================================
  SUBROUTINE init_nbors_hash (self, nbors)
    class(atom_cont_patch_t):: self
    class(nbor_t), pointer:: nbors
    !.............................................................................
    class(task_t), pointer:: task
    class(nbor_t), pointer:: nbor
    integer, save:: itimer=0
    integer :: key(4), i_mu, mu_start, mu_stop
    !-----------------------------------------------------------------------------
    call trace%begin('atom_cont_patch_t%init_nbors_hash', itimer=itimer)
    call self%m3d_patch_t%init_nbors_hash(nbors)
  
    if (comoving_frame .and. .not. self%out_ray) then
      mu_start = 0
      mu_stop  = 0
    else
      mu_start = merge(-self%n_out, 1, self%out_ray)
      mu_stop  = merge(-1, self%n_mu, self%out_ray)
    endif
  
    key(1:2) = self%key(1:2)
    do i_mu=mu_start, mu_stop
      key(3) = i_mu
      key(4) = hash_string('net_opac')
      call self%check_nbor(key, nbors)
    end do
    call trace%end (itimer)
  END SUBROUTINE init_nbors_hash
  
  END MODULE atom_cont_patch_mod
  