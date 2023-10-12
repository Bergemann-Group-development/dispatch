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
    USE m3d_halo_mod, only:lock_mem, unlock_mem
    USE m3d_patch_mod, only:maxiter, m3d_patch_t, comoving_frame, ilambd, preconditioning, save_Qrad
    USE ray_patch_mod, only:ray_patch_t
    USE angle_rad_patch_mod, only:rad_patch_t
    implicit none
    private
  
    !-----------------------------------------------------------------------------
    type, public, extends(rad_patch_t) :: mean_rad_patch_t
      real(4), dimension(:), contiguous, pointer :: old_Jnu => null()
    contains
      procedure :: pre_init
      procedure :: init
      procedure :: execute
      procedure :: nbor_relations
      procedure :: init_nbors_hash
      procedure :: allocate_mem
      procedure :: deallocate_mem
    end type mean_rad_patch_t

    real(4), parameter :: sp_convlim = 0
      
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

    nw = 3
    call self%set_nv_and_nw (self%fe-self%fs+1, nw)

    self%ng = 0
    call self%post_pre_init()
    call trace%end(itimer)
  END SUBROUTINE pre_init

  !===============================================================================
  SUBROUTINE init (self)
    class(mean_rad_patch_t) :: self
    ! ----------------------------------------------------------------------------
    call self%ray_patch_t%init()
    if (self%patch%boundaries%is_set(bits%zu)) call self%boundaries%set(bits%zu)
    if (self%patch%boundaries%is_set(bits%zl)) call self%boundaries%set(bits%zl)
    call self%allocate_mem(self)
    if (sp_convlim /= 0 .and. self%patch%boundaries%is_set(bits%zu)) then
      call os%allocate(self%old_Jnu, self%nv, label=trim(self%kind)//'%old_Jnu')
      self%old_Jnu = huge(1.0)
    endif
  END SUBROUTINE init

  !===============================================================================
  SUBROUTINE allocate_mem (self, allocator, m_in, label)
    class(mean_rad_patch_t):: self
    class(m3d_patch_t), intent(in), target :: allocator
    integer, intent(in), optional :: m_in(6)
    character(len=*), optional :: label
    !-----------------------------------------------------------------------------
    call self%ray_patch_t%allocate_mem(allocator, m_in)
    self%Psi  => self%mem(:,:,:,:,self%it,2)
    self%Ieff => self%mem(:,:,:,:,self%it,3)
    if (save_Qrad) then
      call os%allocate(self%Qrad, self%n, label=trim(self%kind)//'%Qrad')
      self%Qrad = 0.d0
    endif
  END SUBROUTINE allocate_mem

  !===============================================================================
  SUBROUTINE deallocate_mem (self, deallocator, print)
    class(mean_rad_patch_t):: self
    class(m3d_patch_t), intent(in), target :: deallocator
    logical, optional :: print
    !-----------------------------------------------------------------------------
    call self%ray_patch_t%deallocate_mem(self)
    if (associated(self%Qrad)) then
      call os%deallocate(self%Qrad, label=trim(self%kind)//'%Qrad')
    endif
  END SUBROUTINE deallocate_mem
  
  !===============================================================================
  SUBROUTINE execute (self)
    class(mean_rad_patch_t) :: self
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: nbtask
    real(4), dimension(:,:,:,:), pointer :: Inu
    real(4), dimension(:,:,:), pointer :: ss, eta_back, xx
    real(8) :: dt, ex
    real(4) :: Psi, xx0, xxup
    integer :: i,j,k,iw, it, f
    real(4) :: wmu, new_Jnu, sp_conv, max_conv
    logical :: sp_converged
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin (trim(self%kind)//'%update', itimer=itimer)
    Inu => self%authority_mem(:,:,:,:,self%it,1)
    if (comoving_frame) then
      call lock_mem( self%opac, 'mean_rad_update')
      call self%check_access(self%opac)

      do f=1, self%nv
        xx => self%opac%mem(:,:,:,f,self%it,self%i_chi)
        self%Psi(:,:,:,f) = self%Psi(:,:,:,f) / xx
      end do
      call unlock_mem( self%opac, 'mean_rad_update', read=.true.)

      if (preconditioning) then
        call lock_mem( self%back_cont, 'mean_rad_update', read=.true.)
        call self%check_access(self%back_cont)
      endif

      do f=1, self%nv
        if (preconditioning) then
          eta_back => self%back_cont%mem(:,:,:,f,self%it,self%i_eta)
          self%Ieff(:,:,:,f) = self%Ieff(:,:,:,f) + self%Psi(:,:,:,f) * eta_back
        else
          self%Ieff(:,:,:,f) = Inu(:,:,:,f)
        endif
      end do

      if (preconditioning) call unlock_mem( self%back_cont, 'mean_rad_update', read=.true.)
    endif
    if (associated(self%old_Jnu)) then
      max_conv = 0.
      do f=1, self%nv
        new_Jnu = sum(Inu(:,:,self%n(3),f))
        sp_conv = abs(new_Jnu / self%old_Jnu(f) - 1.)
        max_conv = max(max_conv, sp_conv)
        self%old_Jnu(f) = new_Jnu
      end do
      sp_converged = max_conv < sp_convlim
      if (sp_converged) then
        print *, 'Spectrum has converged, i_nu = ', self%i_nu
      endif
      print *, 'Spectrum max_conv, i_nu = ', max_conv, self%i_nu
    endif
    call trace%end(itimer)
  END SUBROUTINE execute
  
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
  