MODULE back_cont_patch_mod
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
    USE m3d_patch_mod, only: maxiter, m3d_patch_t, save_back_cont, hash_string, comoving_frame, ilambd
    USE ray_patch_mod, only: ray_patch_t
    USE composition_mod, only: composition_type
    implicit none
    private
  
    !-----------------------------------------------------------------------------
    type, public, extends(ray_patch_t) :: back_cont_patch_t
    integer :: i_abs=3, i_sca=4
    class(composition_type), pointer :: composition
    class(m3d_patch_t), pointer :: mean_rad => null()
    contains
      procedure :: pre_init => back_cont_pre_init
      procedure :: update => back_cont_update
      procedure :: nbor_relations
      procedure :: init_nbors_hash
    end type back_cont_patch_t
  
    !-----------------------------------------------------------------------------
    integer, save:: verbose=0                           ! verbosity
      
  contains
  
  !===============================================================================
  SUBROUTINE back_cont_pre_init (self)
    class(back_cont_patch_t), target:: self
    integer :: nw, fs, fe
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin('back_cont_patch_t%pre_init', itimer=itimer)
    self%kind = 'back_cont'
    call self%ray_patch_t%pre_init()

    self%static_mem = save_back_cont
  
    self%ng = 0
    nw = merge(4, 2, save_back_cont)
    call self%set_nv_and_nw (self%fe-self%fs+1,nw)
    call self%post_pre_init()
    call trace%end(itimer)
  END SUBROUTINE back_cont_pre_init

  !===============================================================================
  SUBROUTINE back_cont_update (self)
    class(back_cont_patch_t) :: self
    real(4), dimension(:,:,:,:), pointer :: chi, eta
    real(8), dimension(:), pointer :: nu
    real(4), dimension(:,:,:,:), pointer, contiguous :: Jold, abslte, sumsca
    real(4), dimension(:,:,:), pointer, contiguous :: vproj
    character(len=64) :: text
    integer, save:: itimer=0
    integer :: m(6)
    !-----------------------------------------------------------------------------
    if (self%is_finished()) return
    call trace%begin('back_cont_patch_t%update', itimer=itimer)
    chi => self%authority_mem(:,:,:,:, self%it, self%i_chi)
    eta => self%authority_mem(:,:,:,:, self%it, self%i_eta)
    nu  => self%sp%nu(self%fs:self%fe)
    eta = 0
  
    if (save_back_cont) then
      abslte => self%mem(:,:,:,:,self%it,self%i_abs)
      sumsca => self%mem(:,:,:,:,self%it,self%i_sca)
    else
      m = shape(self%mem)
      call os%allocate(abslte, m(1:4), label='abslte')
      call os%allocate(sumsca, m(1:4), label='sumsca')
    endif
  
    if (self%iter == 0 .or. .not. save_back_cont) then
      call self%composition%absorption(nu, chi, abslte, sumsca, include_absmet=self%out_ray==.false.)
    endif
    
    if (.not. self%out_ray .and. ilambd > 0 .and. self%iter > 0) then
      call self%check_access(self%mean_rad)
      Jold => self%mean_rad%authority_mem(:,:,:,:,self%mean_rad%it,1)
      call self%composition%scattering(nu, abslte, eta, Jold, sumsca)
    else
      call self%composition%scattering(nu, abslte, eta)
    endif

    if (self%composition%ll_initialised .and. comoving_frame) then
      call os%allocate(vproj, self%n)
      call self%comp_vproj(vproj)
      call self%composition%add_linelist_opacities(nu, eta, chi, vproj)
      call os%deallocate(vproj)
    endif
  
    if (.not. save_back_cont) then
      call os%deallocate(abslte, label='abslte')
      call os%deallocate(sumsca, label='sumsca')
    endif
  
    call trace%end(itimer)
    call self%post_update()
  END SUBROUTINE back_cont_update
  
  !===============================================================================
  SUBROUTINE nbor_relations(self, another, needed, needs_me, download)
    class(back_cont_patch_t) :: self
    class(task_t), pointer :: another        ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                             ! does it need me, should I download?
    class(task_t), pointer :: other
    logical :: default
    ! ----------------------------------------------------------------------------
    call trace%begin ('back_cont_patch_t%nbor_relations')
    call self%ray_patch_t%nbor_relations(another, needed, needs_me, download)
  
    if (needed .or. needs_me) then
      select type(another)
      class is (m3d_patch_t)
        if (another%id == self%m3d_parent%id) then
          needs_me = .false.
        elseif (another%kind == 'net_opac') then
          needed=.false.
        elseif(another%kind=='angle_rad' .and. .not. self%out_ray) then
          needed = .false.
        elseif(another%kind=='mean_rad') then
          if (.not. comoving_frame) needs_me = .false.
          if (ilambd == 0) needed = .false.
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
    class(back_cont_patch_t):: self
    class(nbor_t), pointer:: nbors
    !.............................................................................
    class(task_t), pointer:: task
    class(nbor_t), pointer:: nbor
    integer, save:: itimer=0
    integer :: key(4), i_mu, mu_start, mu_stop
    !-----------------------------------------------------------------------------
    call trace%begin('back_cont_patch_t%init_nbors_hash', itimer=itimer)
    call self%m3d_patch_t%init_nbors_hash(nbors)
  
    key(1:3) = self%key(1:3)
    mu_start = merge(self%n_omega+1, 1, self%out_ray)
    mu_stop  = merge(self%n_mu, self%n_omega, self%out_ray)

    if (comoving_frame .or. ilambd > 0) then
      key(4) = hash_string('mean_rad')
      call self%check_nbor(key, nbors)
    endif
  
    if (comoving_frame .and. .not. self%out_ray) then
      key(4) = hash_string('net_opac')
      call self%check_nbor(key, nbors)
    else
      do i_mu=mu_start, mu_stop
        key(3) = i_mu
        key(4) = hash_string('net_opac')
        call self%check_nbor(key, nbors)
      end do
    endif
    do i_mu=mu_start, mu_stop
      if (.not. self%out_ray) then
        key(3) = i_mu
        key(4) = hash_string('angle_rad')
        call self%check_nbor(key, nbors)
      endif
    end do
    call trace%end (itimer)
  END SUBROUTINE init_nbors_hash
  
  END MODULE back_cont_patch_mod
  