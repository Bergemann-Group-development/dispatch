MODULE mean_rad_patch_mod
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
    USE m3d_patch_mod, only:maxiter, m3d_patch_t, comoving_frame, ilambd
    USE ray_patch_mod, only:ray_patch_t
    USE angle_rad_patch_mod, only:rad_patch_t
    implicit none
    private
  
    !-----------------------------------------------------------------------------
    type, public, extends(rad_patch_t) :: mean_rad_patch_t
    contains
      procedure :: pre_init
      procedure :: init
      procedure :: update
      procedure :: nbor_relations
      procedure :: init_nbors_hash
    end type mean_rad_patch_t
      
  contains
  
  !===============================================================================
  SUBROUTINE pre_init (self)
    class(mean_rad_patch_t), target:: self
    integer :: nw, fs, fe, i
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin('mean_rad_patch_t%pre_init', itimer=itimer)
    self%kind = 'mean_rad'
    call self%ray_patch_t%pre_init()
    call self%set_nv_and_nw (self%fe-self%fs+1, 1)

    if (.not. associated(self%mem_lock)) then
      allocate (self%mem_lock(self%nt))
      do i=1,self%nt
        call self%mem_lock(i)%init ('meml')
      end do
    endif

    self%ng = 0
    call self%post_pre_init()
    call trace%end(itimer)
  END SUBROUTINE pre_init

  !===============================================================================
  SUBROUTINE init (self)
    class(mean_rad_patch_t) :: self
    ! ----------------------------------------------------------------------------
    call self%ray_patch_t%init()
    call self%allocate_mem(self)
  END SUBROUTINE init
  
  !===============================================================================
  SUBROUTINE update (self)
    class(mean_rad_patch_t) :: self
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: nbtask
    real(4), dimension(:,:,:,:), pointer :: Inu, Jnu
    real(4), dimension(:,:,:), pointer :: ss, eta_back, xx
    real(8) :: dt, ex
    real(4) :: Psi, xx0, xxup
    integer, save:: itimer=0
    integer :: i,j,k,iw, it, f
    real(4) :: wmu
    !-----------------------------------------------------------------------------
    if (self%is_finished()) return
    call trace%begin('mean_rad_patch_t%update', itimer=itimer)
    if (comoving_frame) then
      call self%check_access(self%opac)
      Jnu => self%mem(:,:,:,:,self%it,1)

      do f=1, self%nv
        xx => self%opac%mem(:,:,:,f,self%opac%it,self%opac%i_chi)
        self%Psi(:,:,:,f) = self%Psi(:,:,:,f) / xx
      end do

      do f=1, self%nv
        xx => self%opac%mem(:,:,:,f,self%opac%it,self%opac%i_chi)
        eta_back => self%back_cont%mem(:,:,:,f,self%back_cont%it,self%back_cont%i_eta)
        self%Ieff(:,:,:,f) = self%Ieff(:,:,:,f) + self%Psi(:,:,:,f) * eta_back
      end do
    endif
    call trace%end(itimer)
  END SUBROUTINE update
  
  !===============================================================================
  SUBROUTINE nbor_relations(self, another, needed, needs_me, download)
    class(mean_rad_patch_t) :: self
    class(task_t), pointer :: another        ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                             ! does it need me, should I download?
    class(task_t), pointer :: other
    logical :: default
    ! ----------------------------------------------------------------------------
    call trace%begin ('mean_rad_patch_t%nbor_relations')
    call self%ray_patch_t%nbor_relations(another, needed, needs_me, download)
  
    if (needed .or. needs_me) then
      select type(another)
      class is (ray_patch_t)
        if (another%kind == 'angle_rad') then
            needs_me = .false.
        elseif (another%kind == 'atom_rate') then
            if (self%i_mu /= another%i_mu) needs_me = .false.
            needed = .false.
        elseif (another%kind == 'net_opac') then
          needs_me = .false.
        elseif (another%kind == 'back_cont') then
          !if (.not. comoving_frame) needed = .false.
          !if (ilambd == 0) needs_me = .false.
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
    class(mean_rad_patch_t):: self
    class(nbor_t), pointer:: nbors
    !.............................................................................
    class(task_t), pointer:: task
    class(nbor_t), pointer:: nbor
    integer, save:: itimer=0
    integer :: key(4), i
    !-----------------------------------------------------------------------------
    call trace%begin('mean_rad_patch_t%init_nbors_hash', itimer=itimer)
    call self%ray_patch_t%init_nbors_hash(nbors)
  
    key(1:3) = self%key(1:3)
    key(4) = hash_string('atom_rate')
    call self%check_nbor(key, nbors)

    key(1:3) = self%key(1:3)
    key(4) = hash_string('back_cont')
    call self%check_nbor(key, nbors)

    key(1:3) = self%key(1:3)
    key(4) = hash_string('net_opac')
    call self%check_nbor(key, nbors)

    do i=1, self%n_mu
      key(3) = i
      key(4) = hash_string('angle_rad')
      call self%check_nbor(key, nbors)
    end do
    call trace%end (itimer)
  END SUBROUTINE init_nbors_hash
  
  END MODULE mean_rad_patch_mod
  