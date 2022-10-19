MODULE rate_patch_mod
    USE halo_mod
    USE gpatch_mod
    USE link_mod
    USE task_mod
    USE io_mod
    USE timer_mod
    USE scaling_mod
    USE trace_mod
    USE os_mod
    USE bits_mod
    USE m3d_patch_mod, only: hash_string, comoving_frame
    USE ray_patch_mod, only: ray_patch_t
    USE rad_patch_mod, only: rad_patch_t
    USE trans_mod, only: transition_type
    USE atom_mod, only: atom_type
    USE lines_patch_mod, only: lines_patch_t
    USE mean_rad_patch_mod, only:mean_rad_patch_t
    implicit none

    !-----------------------------------------------------------------------------

    type, public, extends(ray_patch_t) :: rate_patch_t
      class(atom_type), pointer :: atom => null()
      class(rad_patch_t), pointer :: rad => null()
      class(lines_patch_t), pointer :: lines => null()
      integer, dimension(:), pointer, contiguous :: idx_map, idx_remap
      logical :: update_skipped=.false.
      contains
        procedure :: pre_init => rate_patch_pre_init
        procedure :: update => rate_patch_update
        procedure :: nbor_relations => rate_patch_nbor_relations
        procedure :: integrate_rates
        procedure :: print_post_update_info
        procedure :: check_in_bounds
        procedure :: init_nbors_hash
    end type rate_patch_t
  
    !-----------------------------------------------------------------------------
  
  contains
    
  !===============================================================================
  SUBROUTINE rate_patch_pre_init (self)
    class(rate_patch_t), target:: self
    class(transition_type), pointer :: trans
    integer :: i, k, nrad, fs, fe
    logical :: skip
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin('rate_patch_t%rate_patch_pre_init', itimer=itimer)
    self%kind = 'rate_patch'
    call self%ray_patch_t%pre_init()

    call os%allocate(self%idx_map, self%atom%nrad, zero=.true., label='rate_patch:idx_map')
    nrad = 0
    do i=1, size(self%atom%line)
      trans => self%atom%line(i)
      call trans%check_range(self%fs,self%fe,fs,fe,skip)
      if (skip) cycle
      nrad = nrad+1
      self%idx_map(i) = nrad
    end do
    do i=1, size(self%atom%cont)
      trans => self%atom%cont(i)
      call trans%check_range(self%fs,self%fe,fs,fe,skip)
      if (skip) cycle
      nrad = nrad+1
      self%idx_map(i+self%atom%nline) = nrad
    end do
    call os%allocate(self%idx_remap, nrad, label='rate_patch:idx_remap')
    do i=1, self%atom%nrad
      k = self%idx_map(i)
      if (k > 0) self%idx_remap(k) = i
    end do

    call self%set_nv_and_nw(nrad, 2)

    call self%post_pre_init()
    call trace%end(itimer)
  END SUBROUTINE rate_patch_pre_init

  !===============================================================================
  SUBROUTINE rate_patch_update (self)
    class(rate_patch_t) :: self
    real(4), dimension(:,:,:,:), pointer :: gamma_ij, gamma_ji, Inu, Psi, Ieff
    real(4), dimension(:,:,:), pointer, contiguous :: temp
    integer, save:: itimer=0
    !-----------------------------------------------------------------------------
    if (self%is_finished()) return
    call trace%begin('rate_patch_t%rate_patch_update', itimer=itimer)

    gamma_ij => self%mem(:,:,:,:, self%it, self%i_Gij)
    gamma_ji => self%mem(:,:,:,:, self%it, self%i_Gji)

    if (.not. self%rad%stiff) then
      self%static_mem = .false.
    endif

    if (.not. self%static_mem) then
      Psi  => self%rad%Psi
      Ieff => self%rad%Ieff
      Inu  => self%rad%authority_mem(:,:,:,:,self%rad%it,self%rad%i_inu)
      temp => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_temp, self%it,1)

      ! Zero the rates
      gamma_ij = 0.0d0
      gamma_ji = 0.0d0

      call self%check_access(self%lines)
      call self%atom%full_preconditioning(self%fs, temp, Ieff, Psi, gamma_ij, gamma_ji, &
                                      self%idx_map, self%lines%line, self%lines%idx_map)
      self%update_skipped = .false.
    else
      self%update_skipped = .true.
    endif

    call self%integrate_rates(gamma_ij, gamma_ji)

    if (self%rad%stiff) then
      self%static_mem = .true.
      call self%rad%compress_mem(self)
    endif

    call trace%end(itimer)
    call self%post_update()
    if (.not. self%static_mem) then
      call self%deallocate_mem(self, print=.true.)
    endif
  END SUBROUTINE rate_patch_update

  !===============================================================================
  SUBROUTINE integrate_rates (self, gamma_ij, gamma_ji)
    class(rate_patch_t) :: self
    integer :: m(4), i, k
    real(4), dimension(:,:,:,:), pointer :: gamma_ij, gamma_ji
    !-----------------------------------------------------------------------------
    call self%m3d_parent%mem_lock(1)%set('integrate_rates')
    if (.not. associated(self%atom%gamma_ij)) then
      m(1:3) = self%n
      m(4) = self%atom%nrad
      call os%allocate(self%atom%gamma_ij, m, zero=.true., label='atom%gamma_ij')
      call os%allocate(self%atom%gamma_ji, m, zero=.true., label='atom%gamma_ji')
      self%atom%stiff=.true.
    endif
    if (self%i_mu>0) then
      do i=1, self%nv
        k = self%idx_remap(i)
        self%atom%gamma_ij(:,:,:,k) = self%atom%gamma_ij(:,:,:,k) + self%wmu * gamma_ij(:,:,:,i)
        self%atom%gamma_ji(:,:,:,k) = self%atom%gamma_ji(:,:,:,k) + self%wmu * gamma_ji(:,:,:,i)
      end do
    else
      do i=1, self%nv
        k = self%idx_remap(i)
        self%atom%gamma_ij(:,:,:,k) = self%atom%gamma_ij(:,:,:,k) + gamma_ij(:,:,:,i)
        self%atom%gamma_ji(:,:,:,k) = self%atom%gamma_ji(:,:,:,k) + gamma_ji(:,:,:,i)
      end do
    endif
    if (.not. self%static_mem) then
      self%atom%stiff=.false.
    endif
    call self%m3d_parent%mem_lock(1)%unset('integrate_rates')
  END SUBROUTINE integrate_rates
  
  !===============================================================================
  SUBROUTINE rate_patch_nbor_relations(self, another, needed, needs_me, download)
    class(rate_patch_t) :: self
    class(task_t), pointer :: another         ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                             !    does it need me, should I download?
    
    call trace%begin ('rad_patch_t%nbor_relations')
    call self%ray_patch_t%nbor_relations(another, needed, needs_me, download)

    if (needed .or. needs_me) then
      select type (another)
      class is (rad_patch_t)
        if (self%i_mu == another%i_mu) then
          needs_me=.false.
          self%rad => another
        else
          needed = .false.; needs_me = .false.
        endif
      class is (lines_patch_t)
        if (self%i_mu > 0 .or. comoving_frame) then
          needs_me=.false.
          self%lines => another
        else
          needed = .false.; needs_me = .false.
        endif
      class default
        if (another%kind == 'atom_patch') then
          ! is needed and needs me
        else
          needed = .false.; needs_me = .false.
        endif
      end select
    endif
    call trace%end()
  END SUBROUTINE rate_patch_nbor_relations

  !===============================================================================
  SUBROUTINE print_post_update_info (self)
    class(rate_patch_t) :: self
    integer :: p(3), l
    !-----------------------------------------------------------------------------
    p = self%pcheck
    l = self%idx_map(self%linecheck)
    if (self%update_skipped) then
      call self%print_to_stdout(text='update skipped')
    else
      call self%print_to_stdout(arr=self%mem(p(1),p(2),p(3),l,1,:))
    endif
  END SUBROUTINE print_post_update_info
  
  !===============================================================================
  FUNCTION check_in_bounds (self) result (cib)
    class(rate_patch_t) :: self
    logical :: cib
    !-----------------------------------------------------------------------------
    cib = .false.
    if (self%ray_patch_t%check_in_bounds()) then
      if (self%linecheck <= self%atom%nline) then
        cib = self%idx_map(self%linecheck) > 0
      endif
    endif
  END FUNCTION check_in_bounds
  
  !===============================================================================
  !> Fast method to generate nbor lists in cases with very long task lists
  !===============================================================================
  SUBROUTINE init_nbors_hash (self, nbors)
    class(rate_patch_t):: self
    class(nbor_t), pointer:: nbors
    !.............................................................................
    class(task_t), pointer:: task
    class(nbor_t), pointer:: nbor
    integer, save:: itimer=0
    integer :: key(4), i
    !-----------------------------------------------------------------------------
    call trace%begin('rate_patch_t%init_nbors_hash', itimer=itimer)
    call self%m3d_patch_t%init_nbors_hash(nbors)

    key(1:3) = self%key(1:3)
    key(4) = hash_string('lines_patch')
    call self%check_nbor(key, nbors)
  
    if (comoving_frame) then
      key(4) = hash_string('mean_rad')
    else
      key(4) = hash_string('angle_rad')
    endif
    call self%check_nbor(key, nbors)
    call trace%end (itimer)
  END SUBROUTINE init_nbors_hash
  
  END MODULE
