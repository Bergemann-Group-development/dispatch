MODULE atom_rate_patch_mod
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
    USE handyfuncs_mod
    USE m3d_halo_mod, only:lock_mem, unlock_mem
    USE m3d_patch_mod, only: comoving_frame, preconditioning
    USE ray_patch_mod, only: ray_patch_t
    USE angle_rad_patch_mod, only: rad_patch_t
    USE trans_mod, only: transition_type
    USE atom_mod, only: atom_type, atomic_spectrum_type, unique_int
    USE atom_line_patch_mod, only: atom_line_t
    USE mean_rad_patch_mod, only:mean_rad_patch_t
    USE mconst
    USE m3d_debug_mod
    USE math_mod
    implicit none

    !-----------------------------------------------------------------------------

    type, public, extends(ray_patch_t) :: atom_rate_t
      class(rad_patch_t), pointer :: rad => null()
      class(atom_line_t), pointer :: lines => null()
      integer, dimension(:), pointer, contiguous :: idx_map, idx_remap
      logical :: update_skipped=.false.
      contains
        procedure :: pre_init => atom_rate_pre_init
        procedure :: execute
        procedure :: nbor_relations => atom_rate_nbor_relations
        procedure :: integrate_rates
        procedure :: print_post_update_info
        procedure :: check_in_bounds
        procedure :: init_nbors_hash
        procedure :: full_preconditioning
    end type atom_rate_t
  
    !-----------------------------------------------------------------------------
  
  contains
    
  !===============================================================================
  SUBROUTINE atom_rate_pre_init (self)
    class(atom_rate_t), target:: self
    class(transition_type), pointer :: trans
    integer :: i, k, nrad, fs, fe
    logical :: skip
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin('atom_rate_t%atom_rate_pre_init', itimer=itimer)
    self%kind = 'atom_rate'
    call self%ray_patch_t%pre_init()

    call os%allocate(self%idx_map, self%atom%nrad, zero=.true., label='atom_rate:idx_map')

    nrad = 0
    if (self%sp%maxal > 0) then
      do i=1, size(self%atom%line)
        trans => self%atom%line(i)
        if (.not. trans%active) cycle
        call trans%check_range(self%fs,self%fe,fs,fe,skip)
        if (skip) cycle
        nrad = nrad+1
        self%idx_map(trans%ntrans) = nrad
      end do
    endif

    if (self%sp%maxac > 0) then
      do i=1, size(self%atom%cont)
        trans => self%atom%cont(i)
        if (.not. trans%active) cycle
        call trans%check_range(self%fs,self%fe,fs,fe,skip)
        if (skip) cycle
        nrad = nrad+1
        self%idx_map(trans%ntrans) = nrad
      end do
    endif

    call os%allocate(self%idx_remap, nrad, label='atom_rate:idx_remap')
    do i=1, self%atom%nrad
      k = self%idx_map(i)
      if (k > 0) self%idx_remap(k) = i
    end do

    call self%set_nv_and_nw(nrad, 2)

    call self%post_pre_init()
    call trace%end(itimer)
  END SUBROUTINE atom_rate_pre_init

  !===============================================================================
  SUBROUTINE execute (self)
    class(atom_rate_t) :: self
    real(4), dimension(:,:,:,:), pointer :: gamma_ij, gamma_ji
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin (trim(self%kind)//'%update', itimer=itimer)
    call lock_mem(self%rad, 'atom_rate_update', read=.true.)
    call self%check_access(self%rad)
    if (self%sp%maxal > 0) then
      call self%check_access(self%lines)
    endif

    gamma_ij => self%mem(:,:,:,:, self%it, self%i_Gij)
    gamma_ji => self%mem(:,:,:,:, self%it, self%i_Gji)

    if (.not. self%static_mem) then
      call self%full_preconditioning(gamma_ij, gamma_ji)
      self%update_skipped = .false.
    else
      self%update_skipped = .true.
    endif
    
    call unlock_mem(self%rad, 'atom_rate_update', read=.true.)

    call self%integrate_rates(gamma_ij, gamma_ji)
    self%free_self = .not. self%static_mem
    call trace%end(itimer)
  END SUBROUTINE execute
  
  !====================================================================================
  SUBROUTINE full_preconditioning(self, gamma_ij, gamma_ji)
    class(atom_rate_t), target, intent(in) :: self
    real(4), dimension(:,:,:,:), intent(inout), contiguous :: gamma_ij, gamma_ji
    class(atom_type), pointer :: atom
    class(atomic_spectrum_type), pointer :: sp
    real(4), dimension(:,:,:), pointer, contiguous :: temp
    real(4), dimension(:,:,:,:), pointer, contiguous :: Psi, Ieff
    real(4), dimension(:,:,:,:), pointer :: Inu
    real(4), dimension(:,:,:), allocatable :: blzfac, hk_over_temp
    real(4), dimension(:,:,:), pointer, contiguous :: ni, nj
    real(4), dimension(:,:,:,:), allocatable :: SumChiDn, SumChiUp, phi, SumU
    real(4), dimension(:,:,:,:), allocatable :: vij, vji
    real(4), dimension(:), allocatable :: wnu
    integer, dimension(:), allocatable :: lower, upper, trans, ii, jj
    logical, dimension(:), allocatable :: cross
    real(4) :: chi, twohnu3_c2, alpha, Bij, Bji, ufac, alpha2, dn, up, h4pi
    integer :: i, j, k, kr, kr2, kc, kl, fs, fe, n_active, ff, f, p(3)
    integer :: x,y,z, nx, ny, nz, nupper, nlower, i2
    logical :: cross_coupling, skip_admin
    integer, save :: itimer1, itimer2, itimer3, itimer4, itimer5, itimer6
    !------------------------------------------------------------------------------------
    atom => self%atom
    sp => self%sp

    if (self%print_update_info) then
      memcheck = 0
    endif

    ! Zero the rates
    gamma_ij = 0.0d0
    gamma_ji = 0.0d0

    fs = self%fs - self%rad%fs + 1
    fe = self%fe - self%rad%fs + 1
    Psi  => self%rad%Psi(:,:,:,fs:fe)
    Ieff => self%rad%Ieff(:,:,:,fs:fe)
    Inu  => self%rad%authority_mem(:,:,:,fs:fe,self%rad%it,self%rad%i_inu)
    temp => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_temp, self%it,1)

    fs = self%fs
    fe = self%fe
    nx = self%n(1)
    ny = self%n(2)
    nz = self%n(3)

    !n_active = maxval(sp%nac(fs:fe) + sp%nal(fs:fe))
    n_active = sp%nmax

    allocate(wnu(n_active))
    allocate(lower(n_active), upper(n_active), trans(n_active), cross(n_active))
    allocate(ii(n_active), jj(n_active))
    allocate(Vij(nx,ny,nz,n_active))
    allocate(Vji(nx,ny,nz,n_active))
    
    nlower = sp%nlower
    nupper = sp%nupper

    allocate(phi, mold=Psi(:,:,:,1:1))
    allocate(blzfac, hk_over_temp, mold=temp)

    if (preconditioning) then
      allocate(SumChiUp(nx,ny,nz,nlower))
      allocate(SumChiDn(nx,ny,nz,nupper))
      allocate(SumU    (nx,ny,nz,nupper))
    endif

    hk_over_temp = h_k / temp
    h4pi = self%wmu*4*PI / hh

    do ff=fs,fe
        !call trace%begin('administration', itimer=itimer1)
        f = ff - fs + 1

        n_active = sp%nac(ff) + sp%nal(ff)
        trans(1:sp%nac(ff)) = sp%ac(1:sp%nac(ff),ff)
        trans(sp%nac(ff)+1:n_active) = sp%al(1:sp%nal(ff), ff)

        do kr=1, n_active
          k = trans(kr)
          if (kr <= sp%nac(ff)) then
            wnu(kr) = atom%cont(k)%wnu(ff-atom%cont(k)%ired+1)
          else
            wnu(kr) = atom%line(k)%wnu(ff-atom%line(k)%ired+1)
          endif
        end do

        twohnu3_c2 = twoh_c2 * sp%nu(ff)**3
        wnu(1:n_active) = h4pi / sp%nu(ff) * wnu(1:n_active)

        skip_admin = .false.
        if (ff > fs) then
          if (all(sp%al(:,ff) == sp%al(:,ff-1))) then
            if (all(sp%ac(:,ff) == sp%ac(:,ff-1))) then
              skip_admin = .true.
            endif
          endif
        endif

        if (.not. skip_admin) then
          do kr=1, n_active
            k = trans(kr)
            if (kr <= sp%nac(ff)) then
              lower(kr) = atom%cont(k)%i
              upper(kr) = atom%cont(k)%j
            else
              lower(kr) = atom%line(k)%i
              upper(kr) = atom%line(k)%j
            endif
          end do

          if (preconditioning) then 
            cross(1:n_active) = .false.
            cross_coupling = .false.

            do kr=1, n_active
              do kr2=1, n_active
                if (lower(kr) == upper(kr2)) then
                  cross_coupling = .true.
                  cross(kr) = .true.
                endif
              end do
            end do
          endif

          call unique_int(lower(1:n_active), nlower, idx=ii)
          call unique_int(upper(1:n_active), nupper, idx=jj)
        endif

        !call trace%end(itimer1)
        !call trace%begin('nullify', itimer=itimer4)

        if (preconditioning) then
          do i=1,nlower
            SumChiUp(:,:,:,i) = 0
          end do

          do j=1,nupper
            SumU(:,:,:,j) = 0
            if (cross_coupling) SumChiDn(:,:,:,j) = 0
          end do
        endif

        if (sp%nac(ff) > 0) then
          blzfac = exp(- hk_over_temp * real(sp%nu(ff)))
        endif

        !call trace%end(itimer4)
        !call trace%begin('vij', itimer=itimer2)

        do kr=1, n_active
          i = ii(kr)
          j = jj(kr)
          ni => atom%data%n(:,:,:,lower(i))
          nj => atom%data%n(:,:,:,upper(j))

          ufac = twohnu3_c2 / wnu(kr)

          if (kr <= sp%nac(ff)) then
            kc = sp%ac(kr,ff)
            alpha = atom%cont(kc)%alpha(ff-atom%cont(kc)%ired+1) * wnu(kr)
            do z=1,nz
            do y=1,ny
            do x=1,nx
              vij(x,y,z,kr) = alpha
              vji(x,y,z,kr) = alpha * atom%data%ni_over_nj_lte(x,y,z,kc) * blzfac(x,y,z)
              if (preconditioning) then
                chi = vij(x,y,z,kr) * ni(x,y,z) - vji(x,y,z,kr) * nj(x,y,z)
                SumU(x,y,z,j) = SumU(x,y,z,j) + vji(x,y,z,kr) * ufac
                SumChiUp(x,y,z,i) = SumChiUp(x,y,z,i) + chi
                if (cross(kr)) SumChiDn(x,y,z,j) = SumChiDn(x,y,z,j) + chi
              endif
            end do
            end do
            end do
          else
            kl = sp%al(kr - sp%nac(ff),ff)
            kl = self%lines%idx_map(kl)
            call self%lines%line(kl)%interpolate_phi(sp%nu(ff:ff), phi)
            Bij = self%lines%line(kl)%Bij * wnu(kr)
            Bji = self%lines%line(kl)%Bji * wnu(kr)
            do z=1,nz
            do y=1,ny
            do x=1,nx
              vij(x,y,z,kr) = phi(x,y,z,1) * Bij
              vji(x,y,z,kr) = phi(x,y,z,1) * Bji
              if (preconditioning) then
                chi = vij(x,y,z,kr) * ni(x,y,z) - vji(x,y,z,kr) * nj(x,y,z)
                SumU(x,y,z,j) = SumU(x,y,z,j) + vji(x,y,z,kr) * ufac
                SumChiUp(x,y,z,i) = SumChiUp(x,y,z,i) + chi
                if (cross(kr)) SumChiDn(x,y,z,j) = SumChiDn(x,y,z,j) + chi
              endif
            end do
            end do
            end do
          endif
        end do

        !call trace%end(itimer2)
        !call trace%begin('Rij', itimer=itimer3)

        if (preconditioning) then
          do j=1,nupper
            SumU(:,:,:,j) = SumU(:,:,:,j) * Psi(:,:,:,f)
          end do
        endif

        !if (.not. comoving_frame) then
        !  Ieff(:,:,:,f) = Inu(:,:,:,f) - sum(SumU(:,:,:,1:nupper), 4)
        !endif

        do kr=1, n_active
          i = ii(kr)
          j = jj(kr)
          if (cross(kr)) i2 = findloc(upper(1:nupper), lower(i), 1)
          k = trans(kr)
          if (kr <= sp%nac(ff)) k = k + atom%nline
          k = self%idx_map(k)

          if (self%print_update_info) then
            if (f == self%fcheck) then
              if (k==self%linecheck) then
                p = self%pcheck
                memcheck(1) = gamma_ij(p(1),p(2),p(3),k)
                memcheck(2) = gamma_ji(p(1),p(2),p(3),k)
              endif
            endif
          endif

          if (preconditioning) then
            do z=1,nz
            do y=1,ny
            do x=1,nx
              gamma_ji(x,y,z,k) = gamma_ji(x,y,z,k) + vji(x,y,z,kr)*(twohnu3_c2 + Ieff(x,y,z,f)) - SumChiUp(x,y,z,i)*SumU(x,y,z,j)
              if (cross(kr)) then
                gamma_ij(x,y,z,k) = gamma_ij(x,y,z,k) + vij(x,y,z,kr)*Ieff(x,y,z,f) + SumChiDn(x,y,z,j)*SumU(x,y,z,i2)
              else
                gamma_ij(x,y,z,k) = gamma_ij(x,y,z,k) + vij(x,y,z,kr)*Ieff(x,y,z,f)
              endif
            end do
            end do
            end do
          else
            do z=1,nz
            do y=1,ny
            do x=1,nx
              gamma_ji(x,y,z,k) = gamma_ji(x,y,z,k) + vji(x,y,z,kr)*(twohnu3_c2 + Ieff(x,y,z,f))
              gamma_ij(x,y,z,k) = gamma_ij(x,y,z,k) + vij(x,y,z,kr)*Ieff(x,y,z,f)
            end do
            end do
            end do
          endif

          if (self%print_update_info) then
            if (f == self%fcheck) then
              if (k==self%linecheck) then
                p = self%pcheck
                memcheck(1) = gamma_ij(p(1),p(2),p(3),k) - memcheck(1)
                memcheck(2) = gamma_ji(p(1),p(2),p(3),k) - memcheck(2)
              endif
            endif
          endif
        end do

        !call trace%end(itimer3)
    end do
  END SUBROUTINE full_preconditioning
  !------------------------------------------------------------------------------------

  !===============================================================================
  SUBROUTINE integrate_rates (self, gamma_ij, gamma_ji)
    class(atom_rate_t) :: self
    class(atom_type), pointer :: atom
    integer :: m(4), i, k
    real(4), dimension(:,:,:,:), pointer :: gamma_ij, gamma_ji
    real(4) :: wmu
    !-----------------------------------------------------------------------------
    call lock_mem(self%m3d_parent, 'integrate_rates')
    atom => self%atom
    if (.not. associated(atom%data%gamma_ij)) then
      m = [self%n, atom%nrad]
      call os%allocate(atom%data%gamma_ij, m, zero=.true., label='atom%gamma_ij')
      call os%allocate(atom%data%gamma_ji, m, zero=.true., label='atom%gamma_ji')
    endif
    do i=1, self%nv
      k = self%idx_remap(i)
      atom%data%gamma_ij(:,:,:,k) = atom%data%gamma_ij(:,:,:,k) + gamma_ij(:,:,:,i)
      atom%data%gamma_ji(:,:,:,k) = atom%data%gamma_ji(:,:,:,k) + gamma_ji(:,:,:,i)
    end do
    call unlock_mem(self%m3d_parent, 'integrate_rates')
  END SUBROUTINE integrate_rates
  
  !===============================================================================
  SUBROUTINE atom_rate_nbor_relations(self, another, needed, needs_me, download)
    class(atom_rate_t) :: self
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
      class is (atom_line_t)
        if (self%i_mu /= 0 .or. comoving_frame) then
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
  END SUBROUTINE atom_rate_nbor_relations

  !===============================================================================
  SUBROUTINE print_post_update_info (self)
    class(atom_rate_t) :: self
    integer :: p(3)
    !-----------------------------------------------------------------------------
    p = self%pcheck
    if (self%update_skipped) then
      call self%print_to_stdout(text='update skipped')
    elseif (any(memcheck > 0)) then
      call self%print_to_stdout(arr=memcheck(1:2))
    else
      call self%print_to_stdout()
    endif
  END SUBROUTINE print_post_update_info
  
  !===============================================================================
  FUNCTION check_in_bounds (self) result (cib)
    class(atom_rate_t) :: self
    logical :: cib
    !-----------------------------------------------------------------------------
    cib = .false.
    if (self%ray_patch_t%check_in_bounds()) then
      cib = self%idx_map(self%linecheck) > 0
      self%linecheck = self%idx_map(self%linecheck)
    endif
  END FUNCTION check_in_bounds
  
  !===============================================================================
  !> Fast method to generate nbor lists in cases with very long task lists
  !===============================================================================
  SUBROUTINE init_nbors_hash (self, nbors)
    class(atom_rate_t):: self
    class(nbor_t), pointer:: nbors
    !.............................................................................
    class(task_t), pointer:: task
    class(nbor_t), pointer:: nbor
    integer, save:: itimer=0
    integer :: key(4), i
    !-----------------------------------------------------------------------------
    call trace%begin('atom_rate_t%init_nbors_hash', itimer=itimer)
    call self%m3d_patch_t%init_nbors_hash(nbors)

    key(1:3) = self%key(1:3)
    key(4) = hash_string('atom_line')
    call self%check_nbor(key, nbors)
  
    if (self%i_mu == 0) then
      key(4) = hash_string('mean_rad')
    else
      key(4) = hash_string('angle_rad')
    endif
    call self%check_nbor(key, nbors)
    call trace%end (itimer)
  END SUBROUTINE init_nbors_hash
  
  END MODULE
