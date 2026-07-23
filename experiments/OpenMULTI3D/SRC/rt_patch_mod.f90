MODULE rt_patch_mod
  USE os_mod
  USE halo_mod
  USE gpatch_mod
  USE link_mod
  USE task_mod
  USE io_mod
  USE timer_mod
  USE scaling_mod
  USE trace_mod
  USE bits_mod
  USE os_mod
  USE const_mod, only: CC, h_k, twoh_c2, half, one_third, one_sixth, one_twelfth, HZ_TO_AA, pi
  USE patch_mod, only: patch_t
  USE mmath_mod, only: cc_kernel
  USE interpol_mod, only: pchip_mono
  USE handyfuncs_mod
  USE atom_mod, only: atom_type, atom_data_t, atom_trans_t
  USE ray_patch_mod, only: ray_patch_t
  USE m3d_halo_mod, only: lock_mem, unlock_mem
  USE m3d_patch_mod, only: use_vel, m3d_patch_t, &
    dlam, &
    premature_finish, m3d_patch_pre_update=>pre_update, save_resolved, &
    column_by_column, rotate_atmos, save_cntrbf, depression_cntrbf, &
    decouple_continuum, initial_Jnu,  &
    Jnu_to_file, preconditioning, make_eos, amr, save_Inu, save_chi, save_Snu, &
    fixed_Jnu, n_nu, save_Jnu, &
    convlim, converged, flux_convergence, flux_conv_signal, &
    global_flux_err, global_flux_pos, global_flux_lambda, flux_update_counter, &
    maser_chi_floor, maser_dtau_max, maser_Inu_max, suppress_masers
  USE rt_params_mod, only: rt_params_t, rt_iter, rt_out, ray_offset, pivot_z
  USE ray_mod, only: one_dimensional, ray_type, rstar, rorigin
  USE trans_mod, only: line_type, cont_type
  USE bbtrans_mod, only: voigt_function2, voigt_profile, tolerance
  USE linelist_mod, only: ll_add_opacity => add_opacity, ll_add_rates => add_rates, &
    ll_add_opacity_equi => add_opacity_equi, voigt_peak
  USE element_mod, only: element_type
  USE molecule_mod, only: molec_type
  USE ion_mod, only: particle_type
  USE bftrans_mod, only: H_occupation_probability
  USE hprofl_mod, only: hprofile
  USE m3d_debug_mod
  USE shared_mod
  USE composition_mod, only: composition_type, scattering_absorption, ref_lam, &
    departure_pointer_t, read_departures, &
    line_scattering
  USE osmet_mod, only: compute_epsilon_stage
  USE spectrum_mod, only: spectrum_type
  USE optab_mod, only: sample_optab_type
  USE omp_lock_mod, only: lock_t
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer :: n_negative_chi = 0                 ! diagnostic: cells with chi<=0 in eta_to_Snu

  !-----------------------------------------------------------------------------
  type, public, extends(ray_patch_t) :: rt_patch_t
    class(rt_patch_t), pointer :: upstream=>null()
    class(rt_patch_t), pointer :: dnstream=>null()
    class(rt_patch_t), pointer :: up_twin=>null()
    class(rt_patch_t), pointer :: dn_twin=>null()
    logical :: Jnu_freed=.false.
    class(atom_trans_t), dimension(:), pointer :: trans=>null()
    type(spectrum_type) :: sub_sp
    real(4), dimension(:,:,:,:,:), allocatable :: Jnu
    real(4), dimension(:,:,:,:), allocatable :: Psi, PsiJ
    ! points into the composition patch's shared table for this i_nu, or, on the
    ! make_eos path, at a private table (save_opac_table overwrites it in place)
    type(sample_optab_type), pointer :: opac => null()
    type(lock_t), pointer :: opac_lock => null()
    ! Separate sampled LTE line-list table (sample_linelist): non-departure
    ! species only, interpolated into the foreground work arrays with a per-cell
    ! Doppler shift; shares one instance per i_nu across i_mu twins, like opac.
    type(sample_optab_type), pointer :: ll_opac => null()
    type(lock_t), pointer :: ll_opac_lock => null()
    real(4), dimension(:), allocatable :: twohnu3_c2
    real(8), dimension(:), pointer :: nu
    integer :: mus, mue, n_set, n_mu, natom
    logical :: bound_bound, bound_free
    logical, dimension(:), allocatable :: get_rates
    integer :: nx, ny, nz
    real(4), dimension(:,:,:,:,:), allocatable :: work, back
    ! line_scattering: per-(f,x,y,z) accumulator of the (1-eps)-weighted LTE
    ! line opacity built during foreground, and (direct path only) the per-stage
    ! scattering fraction 1-eps on the patch frequency grid (index 1=neutral,
    ! 2=ionized), refreshed per nu-chunk
    real(4), dimension(:,:,:,:),   allocatable :: chi_scat
    real(4), dimension(:,:,:,:,:), allocatable :: one_m_eps
    real(4), dimension(:,:,:,:,:,:), allocatable :: cntrbf_data
    real(4), dimension(:,:,:,:,:), allocatable :: cntrbf_tmp
    ! SoA line kernel path for model-atom bb transitions (see linelist_mod)
    logical, dimension(:), allocatable :: use_ll        ! per atom: SoA kernel path active
    real(4), dimension(:,:,:,:), allocatable :: logref  ! log(tol*cont/voigt_peak) per set
    real(4), dimension(:,:,:,:), allocatable :: b4      ! real(4) copy of adata%b
    class(atom_data_t), pointer :: b4_ptr => null()     ! adata the b4 copy belongs to
    real(8), dimension(:), allocatable :: llam          ! log(lambda[AA]) of the patch grid
    real(8), dimension(:), allocatable :: flux_prev       ! previous-iteration top emergent flux, per frequency (fs:fe)
    real(8) :: flux_global=0.d0, flux_maxerr=0.d0          ! globally-agreed / local worst relative flux change this iteration
  contains
    procedure :: pre_init => rt_patch_pre_init
    procedure :: init_patch_trans
    procedure :: pre_update => rt_patch_pre_update
    procedure :: execute
    procedure :: check_flux_convergence
    procedure :: background
    procedure :: foreground
    procedure :: aniso_integrations
    procedure :: iso_integrations
    procedure :: save_debug_arrays
    procedure :: radiative_transfer
    procedure :: lte_background_opacities
    procedure :: lte_background_emission
    procedure :: generate_opacity_table
    procedure :: generate_ll_opacity_table
    procedure :: save_opac_table
    procedure :: add_bound_free_opacities
    procedure :: bound_bound_damping
    procedure :: add_bound_bound_opacities
    procedure :: bound_bound_rates
    procedure :: bound_free_rates
    procedure :: fill_b4
    procedure :: reduce_ll_diag
    procedure :: convert_to_work_units
    procedure :: reset_array
    procedure :: back_to_work
    procedure :: eta_to_Snu
    procedure :: Snu_to_eta
    procedure :: set_upstream_boundary
    procedure :: integrate_intensity
    procedure :: integrate_rates
    procedure :: set_downstream_boundary
    procedure :: save_intensity
    procedure :: save_contribution_function
    procedure :: rotate_array
    procedure :: short_char
    procedure :: set_fixed_Jnu
    procedure :: long_char_linear
    procedure :: long_char_linear_spherical
    procedure :: long_char_quadratic
    procedure :: long_char_hermite
    procedure :: long_char_hermite2
    procedure :: short_char_hermite
    procedure :: get_vproj
    procedure :: print_post_update_info
    procedure :: nbor_relations => rt_patch_nbor_relations
    procedure :: init_nbors_hash
    procedure :: finish
    procedure :: allocate_mem
    procedure :: copy_to
  end type rt_patch_t

  !-----------------------------------------------------------------------------
  integer, save:: verbose=0                           ! verbosity
  integer, public :: n_cntrbf=0
  integer, public, dimension(:), allocatable :: cntrbf_lines
  ! NB: work, back, cntrbf_data, cntrbf_tmp and the patch dimensions nx,ny,nz are
  ! now components of rt_patch_t (owned per patch); the loop indices x,y,z,f are
  ! plain locals in each routine. They used to live here as module threadprivate
  ! variables, which on macOS/GCC routed every access through emulated TLS
  ! (pthread_getspecific). One thread per patch makes instance/local data both
  ! race-free and TLS-free.

  !-----------------------------------------------------------------------------
  !> Formal solvers, implemented in the rt_solvers_mod submodule (rt_solvers_mod.f90)
  interface
    module subroutine short_char (self, vec)
      class(rt_patch_t) :: self
      real(8), intent(in) :: vec(3)
    end subroutine short_char

    module subroutine long_char_linear_spherical (self, vec)
      class(rt_patch_t) :: self
      real(8) :: vec(3)
    end subroutine long_char_linear_spherical

    module subroutine long_char_linear (self, vec)
      class(rt_patch_t) :: self
      real(8) :: vec(3)
    end subroutine long_char_linear

    pure module subroutine parab_moments (dt, ex, rdt, U0, U1, q)
      real(8), intent(in)  :: dt, ex, rdt
      real(8), intent(out) :: U0, U1, q
    end subroutine parab_moments

    module subroutine long_char_quadratic (self, vec)
      class(rt_patch_t), target :: self
      real(8) :: vec(3)
    end subroutine long_char_quadratic

    module subroutine long_char_hermite (self, comp, vec, comp_up)
      class(rt_patch_t), target :: self
      class(composition_type) :: comp
      class(composition_type), pointer :: comp_up
      real(8), intent(in) :: vec(3)
    end subroutine long_char_hermite

    module subroutine long_char_hermite2 (self, comp, vec, comp_up)
      class(rt_patch_t), target :: self
      class(composition_type) :: comp
      class(composition_type), pointer :: comp_up
      real(8), intent(in) :: vec(3)
    end subroutine long_char_hermite2

    module subroutine short_char_hermite (self, comp, vec)
      class(rt_patch_t) :: self
      class(composition_type) :: comp
      real(8), intent(in) :: vec(3)
    end subroutine short_char_hermite
  end interface

contains

  !===============================================================================
  SUBROUTINE rt_patch_pre_init (self)
    class(rt_patch_t), target:: self
    integer, save :: itimer=0
    logical :: is_nlte
    ! ----------------------------------------------------------------------------
    call trace%begin('rt_patch_t%_pre_init', itimer=itimer)
    self%kind = 'rt_patch'
    call self%m3d_patch_t%pre_init()

    self%dir = merge(1, -1, abs(self%i_mu)==1)
    self%n_nu = n_nu
    self%n_mu = merge(size(rt_out%wmu), size(rt_iter%wmu)/2, self%do_final)
    self%mus = merge(1, self%n_mu+1, self%dir>0)
    self%mue = merge(self%n_mu, size(rt_iter%wmu), self%dir>0)
    self%natom = merge(size(self%atoms), 0, associated(self%atoms))
    is_nlte = rt_iter%maxiter > 0 .or. rt_out%maxiter > 0 .or. read_departures

    if (self%do_out_sp) then
      ! Set 1 is NLTE (or LTE)
      ! Set 2 is Continuum
      ! Set 3 is LTE if set 1 is NLTE, otherwise it is not used
      self%n_set = merge(3, 2, is_nlte)
    else
      ! For the ALI we do not need to compute the continuum or the LTE set
      self%n_set = 1
    endif

    call self%set_freq_indices()
    call self%init_patch_trans()
    allocate(self%get_rates(self%natom))
    self%get_rates = self%do_full_sp .and. self%i_nu>0

    self%ng(3) = 1
    call self%post_pre_init()

    call trace%end(itimer)
  END SUBROUTINE rt_patch_pre_init

  !===============================================================================
  SUBROUTINE init_patch_trans (self)
    class(rt_patch_t) :: self
    class(line_type), pointer :: line
    class(cont_type), pointer :: cont
    integer :: i, j, k, nline, ncont, fs, fe
    logical :: skip
    !-----------------------------------------------------------------------------
    if (.not. self%has_atom) then
      self%bound_bound = .false.
      self%bound_free = .false.
      return
    endif

    self%bound_bound = any([(associated(self%atoms(i)%trans%line), i=1, self%natom)])
    self%bound_free = any([(associated(self%atoms(i)%trans%cont), i=1, self%natom)])
    allocate(self%trans(self%natom))
    allocate(self%use_ll(self%natom))
    ! SoA kernel path for atom lines: disabled for now (no speed gain over the old
    ! formalism for dense model atoms); atom%build_ll still runs globally so the
    ! transitions are exported in linelist_out.txt. Re-enable by restoring
    !   associated(self%atoms(i)%species)
    self%use_ll = .false.

    do k=1, self%natom
      associate(atom=>self%atoms(k), trans=>self%trans(k))
        nline = 0
        ncont = 0
        if (self%bound_bound) then
          do i=1, atom%nline
            line => atom%trans%line(i)
            if (.not. line%active) cycle
            call line%check_range(self%fs, self%fe, fs, fe, skip)
            if (skip) cycle
            nline = nline + 1
          end do
        endif
        if (self%bound_free) then
          do i=1, atom%ncont
            cont => atom%trans%cont(i)
            if (.not. cont%active) cycle
            call cont%check_range(self%fs, self%fe, fs, fe, skip)
            if (skip) cycle
            ncont = ncont + 1
          end do
        endif
        call trans%allocate(nline, ncont)
        if (self%bound_bound) then
          j = 1
          do i=1, atom%nline
            line => atom%trans%line(i)
            if (.not. line%active) cycle
            call line%check_range(self%fs, self%fe, fs, fe, skip)
            if (skip) cycle
            trans%line(j) = line
            trans%line(j)%iblue = fs
            trans%line(j)%ired = fe
            j = j + 1
          end do
        endif
        if (self%bound_free) then
          j = 1
          do i=1, atom%ncont
            cont => atom%trans%cont(i)
            if (.not. cont%active) cycle
            call cont%check_range(self%fs, self%fe, fs, fe, skip)
            if (skip) cycle
            trans%cont(j) = cont
            trans%cont(j)%iblue = fs
            trans%cont(j)%ired = fe
            j = j + 1
          end do
        endif
        ! windowed SoA copy of the (non-hydrogen) patch line subset for the
        ! SoA opacity/rate kernels
        if (self%use_ll(k) .and. self%bound_bound) then
          call atom%build_ll(trans, self%fs)
        endif
      end associate
    end do
  END SUBROUTINE init_patch_trans

  !===============================================================================
  SUBROUTINE rt_patch_pre_update (self)
    class(rt_patch_t), target :: self
    logical :: conv
    integer :: i
    !-----------------------------------------------------------------------------
    call m3d_patch_pre_update(self)
    if (self%has_atom .and. self%do_full_sp) then
      conv = all([(self%atoms(i)%data%converged, i=1, self%natom)])
      if (premature_finish .and. conv) then
        self%static_mem = .true.
      endif
    endif
  END SUBROUTINE rt_patch_pre_update

  !===============================================================================
  SUBROUTINE execute (self)
    class(rt_patch_t) :: self
    class(atom_type), pointer :: atom
    class(atom_trans_t), pointer :: trans
    class(atom_data_t), pointer :: adata
    class(composition_type), pointer :: comp
    real(4), dimension(:,:,:), allocatable :: vproj
    real(4) :: wmu
    real(8) :: vec(3)
    logical :: skip, get_cntrbf, end_boundary, do_rotate_atmos, top_of_domain
    integer :: i, j, i_mu, zs, ze, fs, fe, ng, i_set, nrad
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%execute', itimer=itimer)
    fs = self%fs
    fe = self%fe
    self%nx = self%patch%n(1)
    self%ny = self%patch%n(2)
    self%nz = self%patch%n(3)
    ng = self%ng(3)
    zs = 1 - ng
    ze = self%nz + ng

    self%nu(fs:fe) => self%sp%nu(fs:fe)
    call self%sub_sp%reset(nu=self%nu(fs:fe))
    do_rotate_atmos = rotate_atmos .and. self%do_final
    get_cntrbf = self%do_final .and. save_cntrbf
    self%get_rates = self%get_rates .and. [(.not. self%atoms(i)%data%converged, i=1, self%natom)]
    end_boundary = self%patch%boundaries%is_set(bits%zl) .and. (self%dir<0)
    end_boundary = self%patch%boundaries%is_set(bits%zu) .and. (self%dir>0) .or. end_boundary
    if (.not. end_boundary .and. associated(self%dnstream)) then
      end_boundary = self%dnstream%is_finished()
    endif
    top_of_domain = self%patch%boundaries%is_set(bits%zu) .and. (self%dir > 0)

    call allocate_arrays()

    if (do_rotate_atmos) then ! Final output spectrum from ray-aligned atmosphere
      do i_mu=self%mus, self%mue
        wmu = rt_out%wmu(i_mu)
        vec = rt_out%vec(:,i_mu)
        comp => self%m3d_parent%rotated_composition(i_mu)
        call self%get_vproj(vproj, vec, do_rotate_atmos, i_mu)
        do i_set=1, self%n_set
          call self%background(comp, i_mu, i_set, do_rotate_atmos)
          call self%back_to_work()
          if (i_set /= 2) call self%foreground(comp, i_mu, i_set, do_rotate_atmos, wmu, vproj)
          call self%eta_to_Snu()
          call self%set_upstream_boundary(i_mu, i_set, ng)
          call self%radiative_transfer(comp, vec, do_rotate_atmos, i_mu)
          call self%set_downstream_boundary(i_mu, i_set)
          if (get_cntrbf) call self%save_contribution_function(i_set, i_mu)
          call self%save_debug_arrays(i_set, i_mu)
        end do
      end do
    elseif (self%do_final) then ! Final output spectrum from non-ray-aligned atmosphere
      do i_set=1, self%n_set
        comp => self%m3d_parent%composition
        call self%background(comp, 0, i_set, do_rotate_atmos)
        do i_mu=self%mus, self%mue
          wmu = rt_out%wmu(i_mu)
          vec = rt_out%vec(:,i_mu)
          call self%get_vproj(vproj, vec, do_rotate_atmos, i_mu)
          call self%back_to_work()
          if (i_set /= 2) call self%foreground(comp, i_mu, i_set, do_rotate_atmos, wmu, vproj)
          call self%eta_to_Snu()
          call self%set_upstream_boundary(i_mu, i_set, ng)
          call self%radiative_transfer(comp, vec, do_rotate_atmos, i_mu)
          call self%set_downstream_boundary(i_mu, i_set)
          call self%aniso_integrations(comp, i_mu, i_set, do_rotate_atmos, wmu, vproj)
          if (get_cntrbf) call self%save_contribution_function(i_set, i_mu)
          call self%save_debug_arrays(i_set, i_mu)
        end do
      end do
    else ! ALI iteration for NLTE populations or scattering iteration for final spectrum
      do i_set=1, self%n_set
        comp => self%m3d_parent%composition
        call self%background(comp, 0, i_set, do_rotate_atmos)
        do i_mu=self%mus, self%mue
          wmu = rt_iter%wmu(i_mu)
          vec = rt_iter%vec(:,i_mu)
          call self%get_vproj(vproj, vec, do_rotate_atmos, i_mu)
          if (self%do_full_sp .and. fixed_Jnu > 0) then
            call self%set_fixed_Jnu()
          else
            call self%back_to_work()
            call self%foreground(comp, i_mu, i_set, do_rotate_atmos, wmu, vproj)
            call self%eta_to_Snu()
            call self%set_upstream_boundary(i_mu, i_set, ng)
            call self%radiative_transfer(comp, vec, do_rotate_atmos, i_mu)
          endif
          call self%set_downstream_boundary(i_mu, i_set)
          call self%aniso_integrations(comp, i_mu, i_set, do_rotate_atmos, wmu, vproj)
        end do
        call self%iso_integrations()
      end do
    endif

    if (allocated(self%use_ll) .and. self%do_final) then
      if (any(self%use_ll)) call self%reduce_ll_diag()
    endif

    if (end_boundary) then
      if (self%do_final .and. self%n_mu>0) then
        call self%save_intensity()
        call self%print_linelist(comp, self%n_nu)
      elseif (self%do_full_sp .and. top_of_domain) then
        call self%check_flux_convergence()
      endif
      call self%deallocate_mem()
    endif

    call deallocate_arrays()

    if (self%has_atom .and. self%iter > rt_out%maxiter-1) then
      ! per-atom population convergence (convmax) is an early-exit shortcut for
      ! the population-governed case only -- it must not terminate the task
      ! ahead of the flux metric when flux_convergence is what's authoritative,
      ! or a run can wind down (reach do_final) with flux nowhere near convlim
      ! and no "Converged!" ever printed, simply because a single NG-accelerated
      ! level snapped its local population error below the (very tight) convmax.
      if (all([(self%atoms(i)%data%converged, i=1, self%natom)]) .and. premature_finish &
        .and. .not. flux_convergence) then
        call self%finish()
      endif
    endif

    call trace%end(itimer)

  contains

    !-----------------------------------------------------------------------------
    subroutine allocate_arrays()
      allocate(self%back(fs:fe,self%nx,self%ny,self%nz,0:3), source=0.)
      allocate(self%work(fs:fe,self%nx,self%ny,zs:ze,3), source=0.)
      allocate(self%twohnu3_c2(fs:fe))
      allocate(vproj(self%nx,self%ny,self%nz))
      self%b4_ptr => null()  ! departure coefficients may have changed since last execute
      self%twohnu3_c2 = twoh_c2 * self%nu**3
      call mem_register(self%back,      label='self%back')
      call mem_register(self%work,      label='self%work')
      if (.not. allocated(self%Jnu)) then
        allocate(self%Jnu(fs:fe,self%nx,self%ny,self%nz,self%n_set), source=0.)
        call mem_register(self%Jnu, label='Jnu')
      endif
      do i=1, self%natom
        trans => self%trans(i)
        if (self%bound_bound) then
          ! with the SoA kernel path, adamp is only read by the old-formalism loops,
          ! which then handle just the hydrogen lines
          if (.not. self%use_ll(i) .or. any(trans%line%hydrogen)) then
            call os%allocate(trans%adamp, [self%nx,self%ny,self%nz,size(trans%line)], label='adamp', zero=.true.)
          endif
        endif
        if (self%get_rates(i)) then
          nrad = size(trans%line) + size(trans%cont)
          call os%allocate(trans%gamma_ij, [self%nx,self%ny,self%nz,nrad], label='gamma_ij', zero=.true.)
          call os%allocate(trans%gamma_ji, [self%nx,self%ny,self%nz,nrad], label='gamma_ji', zero=.true.)
        endif
      end do
      if (allocated(self%use_ll)) then
        if (self%bound_bound .and. any(self%use_ll)) then
          allocate(self%logref(fs:fe,self%nx,self%ny,self%nz))
          call mem_register(self%logref, label='self%logref')
          allocate(self%llam(fs:fe))
          self%llam = log(HZ_TO_AA / self%nu(fs:fe))
        endif
      endif
      if (line_scattering) then
        allocate(self%chi_scat(fs:fe,self%nx,self%ny,self%nz), source=0.)
        call mem_register(self%chi_scat, label='self%chi_scat')
        ! one_m_eps is filled per foreground for the detailed line path; the
        ! sampled path bakes epsilon into the table instead
        allocate(self%one_m_eps(fs:fe,self%nx,self%ny,self%nz,2), source=0.)
        call mem_register(self%one_m_eps, label='self%one_m_eps')
      endif
      if (any(self%get_rates)) then
        allocate(self%Psi(fs:fe,self%nx,self%ny,self%nz), source=0.)
        call mem_register(self%Psi, label='Psi')
        if (self%bound_free) then
          allocate(self%PsiJ(fs:fe,self%nx,self%ny,self%nz), source=0.)
          call mem_register(self%PsiJ, label='PsiJ')
        endif
      endif
      if (get_cntrbf) then
        allocate(self%cntrbf_data(fs:fe,self%nx,self%ny,self%nz,2,self%mus:self%mue))
        call mem_register(self%cntrbf_data, 'self%cntrbf_data')
        allocate(self%cntrbf_tmp(fs:fe,self%nx,self%ny,self%nz,2))
        call mem_register(self%cntrbf_tmp, 'self%cntrbf_tmp')
      endif
    end subroutine allocate_arrays

    !-----------------------------------------------------------------------------
    subroutine deallocate_arrays()
      call mem_unregister(self%work     , label='self%work')
      call mem_unregister(self%back     , label='self%back')
      deallocate(self%work, self%back)
      deallocate(self%twohnu3_c2)
      do i=1, self%natom
        trans => self%trans(i)
        if (self%bound_bound .and. allocated(trans%adamp)) then
          call os%deallocate(trans%adamp, label='adamp')
        endif
        if (self%get_rates(i)) then
          call os%deallocate(trans%gamma_ij, label='gamma_ij')
          call os%deallocate(trans%gamma_ji, label='gamma_ji')
        endif
      end do
      if (allocated(self%logref)) then
        call mem_unregister(self%logref, label='self%logref')
        deallocate(self%logref)
        deallocate(self%llam)
      endif
      if (allocated(self%chi_scat)) then
        call mem_unregister(self%chi_scat, label='self%chi_scat')
        call mem_unregister(self%one_m_eps, label='self%one_m_eps')
        deallocate(self%chi_scat, self%one_m_eps)
      endif
      if (allocated(self%b4)) deallocate(self%b4)
      if (allocated(self%Jnu) .and. Jnu_to_file) then
        call mem_unregister(self%Jnu, label='Jnu')
        deallocate(self%Jnu)
      endif
      if (any(self%get_rates)) then
        call mem_unregister(self%Psi, label='Psi')
        deallocate(self%Psi)
        if (self%bound_free) then
          call mem_unregister(self%PsiJ, label='PsiJ')
          deallocate(self%PsiJ)
        endif
      endif
      if (get_cntrbf) then
        call mem_unregister(self%cntrbf_data, label='self%cntrbf_data')
        deallocate(self%cntrbf_data)
        call mem_unregister(self%cntrbf_tmp, label='self%cntrbf_tmp')
        deallocate(self%cntrbf_tmp)
      endif
    end subroutine deallocate_arrays
  END SUBROUTINE execute

  !===============================================================================
  !> The rt_params instance governing this patch: rt_out for any ray working on
  !> the output spectrum, rt_iter for the internal (iteration) spectrum.
  !===============================================================================
  FUNCTION pass_params (self) RESULT (rtp)
    class(rt_patch_t), intent(in) :: self
    type(rt_params_t), pointer :: rtp
    !-----------------------------------------------------------------------------
    if (self%do_out_sp) then
      rtp => rt_out
    else
      rtp => rt_iter
    endif
  END FUNCTION pass_params

  !===============================================================================
  SUBROUTINE background(self, comp, i_mu, i_set, do_rotate_atmos)
    class(rt_patch_t) :: self
    class(composition_type), pointer :: comp
    class(atom_data_t), pointer :: adata
    integer, intent(in) :: i_mu, i_set
    logical, intent(in) :: do_rotate_atmos
    type(rt_params_t), pointer :: rtp
    real(8) :: vec(3)
    logical :: include_absmet, do_sample_opac, do_make_eos
    integer :: i
    !-----------------------------------------------------------------------------
    rtp => pass_params(self)
    ! with a decoupled continuum, only the bound-free half of the internal
    ! spectrum (i_nu > n_nu) can carry the absmet background
    include_absmet = rtp%include_absmet .and. (self%do_out_sp .or. &
                     self%i_nu > self%n_nu .or. .not. decouple_continuum)
    do_sample_opac = rtp%sample_back_opac
    do_make_eos = make_eos .and. self%do_out_sp

    if (do_rotate_atmos .and. i_mu /= 0) then
      vec = rt_out%vec(:,i_mu)
    endif

    if (i_set==1) then
      call self%lte_background_opacities(comp, include_absmet, do_sample_opac, do_make_eos)
    endif
    call self%lte_background_emission(comp, i_set, vec, include_absmet, do_rotate_atmos)
    call self%convert_to_work_units()
    call self%reset_array(self%Jnu(:,:,:,:,i_set:i_set))
    call self%reset_array(self%work)
    
    do i=1, self%natom
      associate(atom=>self%atoms(i), trans=>self%trans(i))
        adata => atom%data
        if (do_rotate_atmos .and. i_mu /= 0) then
          adata => atom%rotated_data(i_mu)
        endif

        if (self%bound_free) then
          call self%add_bound_free_opacities(atom, trans, adata)
        endif
        ! adamp is only allocated where the old-formalism loops still need it
        ! (atoms with hydrogen lines)
        if (self%bound_bound .and. i_set==1 .and. allocated(trans%adamp)) then
          call self%bound_bound_damping(comp, atom, trans, adata)
        endif
      end associate
    end do

    ! log-space tolerance reference for the SoA line kernels; the continuum (back
    ! work units, incl. bound-free) is fixed per i_set, so compute it once here and
    ! reuse it for all mu directions
    if (allocated(self%logref) .and. i_set /= 2) then
      call fill_logref()
    endif

  contains
    subroutine fill_logref()
      real(4) :: log_tol, b
      integer :: x, y, z, f
      integer, save :: itimer=0
      call trace%begin('rt_patch%fill_logref', itimer=itimer)
      log_tol = log(tolerance / voigt_peak)
      do z=1, self%nz
        do y=1, self%ny
          do x=1, self%nx
            do f=self%fs, self%fe
              b = self%back(f,x,y,z,1)
              if (b > 0.) then
                self%logref(f,x,y,z) = log_tol + log(b)
              else
                self%logref(f,x,y,z) = -huge(1.)
              endif
            end do
          end do
        end do
      end do
      call trace%end(itimer)
    end subroutine fill_logref
  END SUBROUTINE background

  !===============================================================================
  SUBROUTINE foreground(self, comp, i_mu, i_set, do_rotate_atmos, wmu, vproj)
    class(rt_patch_t) :: self
    class(composition_type), pointer :: comp
    class(atom_data_t), pointer :: adata
    class(departure_pointer_t), pointer :: departures(:)
    integer, intent(in) :: i_mu, i_set
    logical, intent(in) :: do_rotate_atmos
    real(4), intent(in) :: wmu
    real(4), intent(in), dimension(:,:,:) :: vproj
    type(rt_params_t), pointer :: rtp
    logical :: use_linelist, sample_ll, do_scat
    integer :: i, x, y, z, f
    real(4) :: beta, dJ, dEta
    !-----------------------------------------------------------------------------
    ! detailed line list vs sampled LTE line table, per pass (rt_iter or rt_out)
    rtp => pass_params(self)
    use_linelist = comp%ll_initialised .and. rtp%include_linelist .and. .not. rtp%sample_linelist
    sample_ll = rtp%sample_linelist .and. associated(self%ll_opac)
    if (sample_ll) sample_ll = self%ll_opac%initialised

    ! line_scattering: accumulate the (1-eps)-weighted LTE line opacity, then
    ! redirect that fraction of the line source from Bnu to Jnu (see below)
    do_scat = line_scattering .and. (use_linelist .or. sample_ll)
    if (do_scat) then
      self%chi_scat = 0.
      if (use_linelist) then
        ! per-stage scattering fraction on the detailed-path frequency grid
        call compute_epsilon_stage(1, self%sub_sp%nu, comp%temp, comp%ne, &
          self%one_m_eps(:,:,:,:,1))
        call compute_epsilon_stage(2, self%sub_sp%nu, comp%temp, comp%ne, &
          self%one_m_eps(:,:,:,:,2))
      endif
    endif

    do i=1, self%natom
      associate(atom=>self%atoms(i), trans=>self%trans(i))
        adata => atom%data
        if (do_rotate_atmos .and. i_mu /= 0) then
          adata => atom%rotated_data(i_mu)
        endif
        if (self%bound_bound) then
          call self%add_bound_bound_opacities(comp, atom, trans, adata, wmu, vproj, &
                                              LTE=i_set==3, use_ll=self%use_ll(i))
        endif
      end associate
    end do

    ! NLTE departure coefficients for this (i_set, i_mu), shared by both the full
    ! detailed line list and the sampled-mode detailed departure-species pass
    departures => null()
    if (associated(self%m3d_parent%departures) .and. i_set==1) then
      i = merge(i_mu+1, 1, do_rotate_atmos)
      departures => self%m3d_parent%departures(:,i)
    endif

    if (use_linelist) then
      ! full detailed line list (non-sampled mode): every species, with departures
      associate(xhat=>self%work(:,:,:,1:self%nz,1), ehat=>self%work(:,:,:,1:self%nz,2))
        if (do_scat) then
          call comp%add_linelist_opacities2(self%sub_sp, xhat, ehat, dlam, vproj, &
            departures=departures, one_m_eps=self%one_m_eps, chi_scat=self%chi_scat)
        else
          call comp%add_linelist_opacities2(self%sub_sp, xhat, ehat, dlam, vproj, departures=departures)
        endif
      end associate
    elseif (sample_ll) then
      ! (a) non-departure (LTE) species from the sampled table, Doppler-shifted per
      ! cell, added to both work slots (xhat == ehat in LTE); the table also
      ! carries the (1-eps)-weighted scattering component when line_scattering
      if (do_scat) then
        call self%ll_opac%interpolate_ll(comp%temp, comp%ne, comp%chi500, &
          self%work(:,:,:,1:self%nz,1:2), self%sp%lam, vproj, self%fs, self%fe, &
          chi_scat=self%chi_scat)
      else
        call self%ll_opac%interpolate_ll(comp%temp, comp%ne, comp%chi500, &
          self%work(:,:,:,1:self%nz,1:2), self%sp%lam, vproj, self%fs, self%fe)
      endif
      ! (b) departure-carrying species in detail, with departures + velocity
      if (associated(departures)) then
        associate(xhat=>self%work(:,:,:,1:self%nz,1), ehat=>self%work(:,:,:,1:self%nz,2))
          call comp%add_linelist_opacities2(self%sub_sp, xhat, ehat, dlam, vproj, &
            departures=departures, only_departures=.true.)
        end associate
      endif
    endif

    ! line_scattering: the LTE lines were added as pure thermal absorbers (source
    ! Bnu).  Redirect the scattering fraction (1-eps) of their opacity from Bnu
    ! to the mean intensity: Snu_line = eps*Bnu + (1-eps)*Jnu.  With chi_scat the
    ! (1-eps)-weighted raw line opacity, the emission correction Delta_eta =
    ! chi_scat*(1-beta)*(Jnu-Bnu) is added to eta without altering chi, which in
    ! working units means adding Delta_eta to xhat and Delta_eta/beta to ehat
    ! (eta_to_Snu reconstructs eta=ehat*beta, chi=xhat-eta).  Jnu already holds
    ! Jnu-Bnu here when scattering_absorption; this is a Lambda-iteration term
    ! (no effect at maxiter=0 / Planck Jnu), exactly like continuum scattering.
    if (do_scat) then
      do z=1, self%nz
        do y=1, self%ny
          do x=1, self%nx
            do f=self%fs, self%fe
              beta = self%back(f,x,y,z,3)
              dJ   = self%Jnu(f,x,y,z,i_set)
              if (.not. scattering_absorption) dJ = dJ - beta/(1 - beta)
              dEta = self%chi_scat(f,x,y,z) * (1 - beta) * dJ
              self%work(f,x,y,z,1) = self%work(f,x,y,z,1) + dEta
              self%work(f,x,y,z,2) = self%work(f,x,y,z,2) + dEta/beta
            end do
          end do
        end do
      end do
    endif
  END SUBROUTINE foreground

  !===============================================================================
  SUBROUTINE radiative_transfer(self, comp, vec, do_rotate_atmos, i_mu)
    class(rt_patch_t) :: self
    class(composition_type) :: comp
    real(8), intent(in) :: vec(3)
    logical, intent(in) :: do_rotate_atmos
    integer, intent(in) :: i_mu
    class(composition_type), pointer :: comp_up
    !-----------------------------------------------------------------------------
    ! composition of the real upstream neighbour patch, one plane before self's
    ! ghost boundary -- used by long_char_hermite in place of a same-patch proxy
    ! when a genuine neighbour is available (null at a physical domain boundary)
    comp_up => null()
    if (associated(self%upstream)) then
      if (do_rotate_atmos) then
        comp_up => self%upstream%m3d_parent%rotated_composition(i_mu)
      else
        comp_up => self%upstream%m3d_parent%composition
      endif
    endif

    if (self%do_final) then
      if (.not. (one_dimensional .or. column_by_column .or. do_rotate_atmos)) then
        call self%rotate_array(vec, self%work(:,:,:,1:self%nz,1:2))
      endif

      select case(rt_out%solver)
        case('linear')
        !if (one_dimensional .and. rstar > 2) then
        !  call self%long_char_linear_spherical(vec)
        !else
        call self%long_char_linear(vec)
        !endif
        case('quadratic')
        call self%long_char_quadratic(vec)
        case('hermite')
        ! the chain rule needs the composition gradients along the ray; in 3D
        ! that means ray-aligned (rotated) compositions -- the rotate_array
        ! shortcut above rotates only chi/Snu, not the thermodynamic state
        if (.not. (one_dimensional .or. column_by_column .or. do_rotate_atmos)) &
          stop "rt_output_params solver='hermite' in 3D requires rotate_atmos=.true."
        call self%long_char_hermite(comp, vec, comp_up)
        case('hermite2')
        if (.not. (one_dimensional .or. column_by_column .or. do_rotate_atmos)) &
          stop "rt_output_params solver='hermite2' in 3D requires rotate_atmos=.true."
        call self%long_char_hermite2(comp, vec, comp_up)
        case default
        stop 'Unknown rt_output_params solver, use "linear", "quadratic", "hermite" or "hermite2"'
      end select
    else
      if (one_dimensional .or. column_by_column) then
        select case(rt_iter%solver)
          case('linear')
            call self%long_char_linear(vec)
          case('quadratic')
            call self%long_char_quadratic(vec)
          case('hermite')
            call self%long_char_hermite(comp, vec, comp_up)
          case('hermite2')
            call self%long_char_hermite2(comp, vec, comp_up)
          case default
            stop 'Unknown rt_iteration_params solver, use "linear", "quadratic", "hermite" or "hermite2"'
        end select
      else
        select case(rt_iter%solver)
          case('linear')
            call self%short_char(vec)
          case('hermite', 'quadratic')
            call self%short_char_hermite(comp, vec)
          case default
            stop 'Unknown rt_iteration_params solver, use "linear", "quadratic" or "hermite"'
        end select
      endif
    endif
  END SUBROUTINE radiative_transfer

  !===============================================================================
  SUBROUTINE aniso_integrations(self, comp, i_mu, i_set, do_rotate_atmos, wmu, vproj)
    class(rt_patch_t) :: self
    class(composition_type), pointer :: comp
    class(atom_data_t), pointer :: adata
    integer, intent(in) :: i_mu, i_set
    logical, intent(in) :: do_rotate_atmos
    real(4), intent(in) :: wmu
    real(4), dimension(:,:,:), intent(in) :: vproj
    logical :: get_Jnu, get_PsiJ
    integer :: i
    !-----------------------------------------------------------------------------
    get_Jnu = ((rt_iter%maxiter > 0) .or. (rt_out%maxiter > 0)) .and. .not. self%do_final
    get_PsiJ = any(self%get_rates) .and. self%bound_free

    if (get_Jnu) then
      call self%integrate_intensity(wmu, i_set, save_PsiJ=get_PsiJ)
    endif

    do i=1, self%natom
      associate (atom=>self%atoms(i), trans=>self%trans(i))
        adata => atom%data
        if (do_rotate_atmos) then
          adata => atom%rotated_data(i_mu)
        endif
        if (self%get_rates(i) .and. self%bound_bound) then
          call self%bound_bound_rates(comp, atom, trans, adata, wmu, vproj, self%use_ll(i))
        endif
      end associate
    end do
  END SUBROUTINE aniso_integrations

  !===============================================================================
  SUBROUTINE iso_integrations(self)
    class(rt_patch_t) :: self
    character(8) :: Jnu_label
    integer :: i
    !-----------------------------------------------------------------------------
    do i=1, self%natom
      associate(atom=>self%atoms(i), trans=>self%trans(i))
        if (self%get_rates(i)) then
          if (self%bound_free) call self%bound_free_rates(atom, trans, atom%data)
          call self%integrate_rates(atom, trans, atom%data)
        endif
      end associate
    end do

    ! twin sub-patches only exist for do_full_sp when rt_iter%maxiter>0, and for
    ! do_final when rt_out%maxiter>0 (see composition_patch_mod.f90 create_sub_patches
    ! calls) -- saving Jnu outside those conditions leaves it with no reader.
    if (Jnu_to_file .and. ((self%do_full_sp .and. rt_iter%maxiter > 0) .or. &
                           (self%do_final .and. rt_out%maxiter > 0))) then
      Jnu_label = merge('Jnu_out', 'Jnu    ', self%do_out_sp)
      call self%save5(self%Jnu, trim(Jnu_label), key=self%key(1:3), add2meta=save_Jnu)
    endif
  END SUBROUTINE iso_integrations

  !===============================================================================
  SUBROUTINE save_debug_arrays(self, i_set, i_mu)
    class(rt_patch_t) :: self
    integer, intent(in) :: i_set, i_mu
    character(6) :: sav_label
    !-----------------------------------------------------------------------------
    if (save_chi) then
      if (i_set==1) then
        write(sav_label, '(A, I2.2)') 'chi_', i_mu
        call self%save4(self%work(:,:,:,1:self%nz,1), sav_label, key=self%key(1:3))
      endif
    endif
    if (save_Snu) then
      if (i_set==1) then
        write(sav_label, '(A, I2.2)') 'Snu_', i_mu
        call self%save4(self%work(:,:,:,1:self%nz,2), sav_label, key=self%key(1:3))
      endif
    endif
    if (save_Inu) then
      if (i_set==1) then
        write(sav_label, '(A, I2.2)') 'Inu_', i_mu
        call self%save4(self%work(:,:,:,1:self%nz,3), sav_label, key=self%key(1:3))
      elseif (i_set==2) then
        write(sav_label, '(A, I2.2)') 'Cnt_', i_mu
        call self%save4(self%work(:,:,:,1:self%nz,3), sav_label, key=self%key(1:3))
      endif
    endif
  END SUBROUTINE save_debug_arrays


  !===============================================================================
  SUBROUTINE lte_background_opacities (self, comp, include_absmet, do_sample_opac, do_make_eos)
    class(rt_patch_t) :: self
    class(composition_type) :: comp
    real(4) :: h_kT, hv_kT, bnu
    type(spectrum_type) :: tab_sp
    type(rt_params_t), pointer :: rtp
    logical :: include_absmet, do_sample_opac, do_make_eos
    integer :: fs, fe
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%lte_background_opacities', itimer=itimer)
    rtp => pass_params(self)
    fs = self%fs
    fe = self%fe

    do z=1, self%nz
      do y=1, self%ny
        do x=1, self%nx
          ! computing boltzmann factor
          h_kT = h_k / real(comp%temp(x,y,z))
          do f=fs, fe
            hv_kT = self%nu(f) * h_kT
            ! Clamp the LTE Boltzmann factor to the smallest positive real(4)
            ! so it never underflows to exactly 0. NLTE emissivity (eta) can
            ! stay finite where beta would underflow, and downstream code forms
            ! eta/beta; a hard zero here produces Inf/NaN (see convert_to_work_units).
            self%back(f,x,y,z,3) = max(exp(-hv_kT), tiny(1.0))
          end do
        end do
      end do
    end do

    ! A table's content depends on i_nu alone, so the rt_patches that differ only
    ! in i_mu share the one the composition patch holds for this i_nu: the first
    ! of them to get here fills it while the others wait on its lock.  The
    ! make_eos path keeps a private table, since save_opac_table overwrites it.
    if (do_make_eos .or. do_sample_opac) then
      if (.not. associated(self%opac)) then
        if (do_make_eos .or. .not. associated(self%m3d_parent%opac_tables)) then
          allocate(self%opac)
        else
          self%opac      => self%m3d_parent%opac_tables(self%i_nu)
          self%opac_lock => self%m3d_parent%opac_locks (self%i_nu)
        endif
      endif

      if (.not. self%opac%initialised) then
        if (associated(self%opac_lock)) call self%opac_lock%set('opac_table')
        if (.not. self%opac%initialised) then     ! a twin may have filled it meanwhile
          if (do_make_eos) then
            call self%generate_opacity_table(self%sub_sp, include_absmet, comp%ll_initialised, &
              build_deriv=.false.)
          else
            if (include_absmet) then
              tab_sp = self%sub_sp
            else
              call comp%get_table_sp(self%sub_sp, tab_sp)
            endif
            ! The background table is continuum-only. With sample_linelist the line
            ! opacity lives in the separate ll_opac table (built below), not baked
            ! into the background, so it can be interpolated with a per-cell Doppler
            ! shift into the foreground rather than the velocity-blind self%back.
            call self%generate_opacity_table(tab_sp, include_absmet, include_llist=.false., &
              build_deriv=.true.)
          endif
        endif
        if (associated(self%opac_lock)) call self%opac_lock%unset('opac_table')
      endif
    endif

    ! Separate LTE line-list opacity table (non-departure species only), keyed by
    ! i_nu like opac and shared across i_mu twins. Built velocity-free here; the
    ! Doppler shift is applied when it is interpolated into the foreground.
    if (do_sample_opac .and. rtp%sample_linelist .and. &
        associated(self%m3d_parent%ll_opac_tables)) then
      if (.not. associated(self%ll_opac)) then
        self%ll_opac      => self%m3d_parent%ll_opac_tables(self%i_nu)
        self%ll_opac_lock => self%m3d_parent%ll_opac_locks (self%i_nu)
      endif
      if (.not. self%ll_opac%initialised) then
        if (associated(self%ll_opac_lock)) call self%ll_opac_lock%set('ll_opac_table')
        if (.not. self%ll_opac%initialised) then   ! a twin may have filled it meanwhile
          call self%generate_ll_opacity_table(self%sub_sp, include_absmet)
        endif
        if (associated(self%ll_opac_lock)) call self%ll_opac_lock%unset('ll_opac_table')
      endif
    endif

    ! getting background LTE absorption
    if (do_make_eos .and. self%iter==0) then
      call self%save_opac_table()
    elseif (do_sample_opac) then
      call self%opac%interpolate(comp%temp, comp%ne, comp%chi500, &
        self%back(:,:,:,:,0), self%sp%lam, self%fs, self%fe)
    else
      call comp%absorption(self%sub_sp, self%back(:,:,:,:,0), &
        beta=self%back(:,:,:,:,3), include_absmet=include_absmet, sorted=self%sp%sorted)
    endif
    call trace%end(itimer)
  END SUBROUTINE lte_background_opacities

  !===============================================================================
  SUBROUTINE lte_background_emission (self, comp, i_set, vec, include_absmet, do_rotate_atmos)
    class(rt_patch_t), target :: self
    class(composition_type) :: comp
    integer, intent(in) :: i_set
    real(8), intent(in) :: vec(3)
    logical, intent(in) :: include_absmet, do_rotate_atmos
    real(4), dimension(:,:,:,:,:), allocatable, target :: Jup_mem, Jdn_mem
    real(4), dimension(:,:,:,:,:), pointer, contiguous :: Jup, Jdn
    real(4) :: Bnu
    character(8) :: Jnu_label
    logical :: Jnu_to_Bnu
    integer :: s, fs, fe
    integer, save :: itimer
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%lte_background_emission', itimer=itimer)
    if (.not. associated(self%up_twin)) self%up_twin => self

    fs = self%fs
    fe = self%fe

    select case (uppercase(initial_Jnu))
    case ('B', 'BNU', 'PLANCK')
      Jnu_to_Bnu = .true.
    case ('0', 'ZERO')
      Jnu_to_Bnu = .false.
    case default
      print *, 'WARNING: Unknown initial_Jnu, use "Planck" or "zero". Setting to "Planck".'
      Jnu_to_Bnu = .true.
    end select

    if (self%iter > 0 .or. (self%do_final .and. rt_out%maxiter > 0)) then
      Jnu_to_Bnu = .false.
      if (Jnu_to_file) then
        Jnu_label = merge('Jnu_out', 'Jnu    ', self%do_out_sp)
        allocate(Jup_mem(fs:fe,self%nx,self%ny,self%nz,self%up_twin%n_set))
        allocate(Jdn_mem(fs:fe,self%nx,self%ny,self%nz,self%dn_twin%n_set))
        call self%up_twin%load5(Jnu_label, key=self%up_twin%key(1:3), mem=Jup_mem)
        call self%dn_twin%load5(Jnu_label, key=self%dn_twin%key(1:3), mem=Jdn_mem)
        Jup => Jup_mem
        Jdn => Jdn_mem
      else
        Jup => self%up_twin%Jnu
        Jdn => self%dn_twin%Jnu
      endif
      do z=1, self%nz
        do y=1, self%ny
          do x=1, self%nx
            do f=fs, fe
              self%Jnu(f,x,y,z,i_set) = 0.5*(Jup(f,x,y,z,i_set) + Jdn(f,x,y,z,i_set))
            end do
          end do
        end do
      end do
      !-------------------------------------------------------------------------
      ! Jnu enters eta via comp%scattering, so a corrupt Jnu -- from the twins'
      ! save/load, whichever way it is passed -- shows up as an infinite source
      ! function.  Separate that from a corrupt boundary (see
      ! set_upstream_boundary) and from a corrupt opacity.
      !-------------------------------------------------------------------------
      if (check_NaN) then
        do z=1, self%nz
          do y=1, self%ny
            do x=1, self%nx
              do f=fs, fe
                if (.not. ieee_is_finite(self%Jnu(f,x,y,z,i_set))) then
                  write(stdout,'(A,I8,A,2(I8),A,4(I4),A,2(ES12.4))') &
                    'RACE: bad Jnu in', self%id, ' twins', self%up_twin%id, self%dn_twin%id, &
                    ' at f,x,y,z=', f, x, y, z, ' Jup,Jdn=', Jup(f,x,y,z,i_set), Jdn(f,x,y,z,i_set)
                  flush(stdout)
                  call shared%abort ('rt_patch_t%lte_background_emission: non-finite Jnu')
                endif
              end do
            end do
          end do
        end do
      endif
      if (do_rotate_atmos .and. .not. (one_dimensional .or. column_by_column)) then
        call self%rotate_array(vec, self%Jnu(:,:,:,:,i_set:i_set))
      endif
    elseif (.not. Jnu_to_Bnu) then
      call self%reset_array(self%Jnu(:,:,:,:,i_set:i_set))
    endif

    do z=1, self%nz
      do y=1, self%ny
        do x=1, self%nx
          do f=fs, fe
            associate(Jnu=>self%Jnu(f,x,y,z,i_set), chi=>self%back(f,x,y,z,0), eta=>self%back(f,x,y,z,2), beta=>self%back(f,x,y,z,3))
              Bnu = beta / (1 - beta)
              if (Jnu_to_Bnu) then
                Jnu = Bnu
              endif
              if (scattering_absorption) then
                ! because we already included scattering opacities in chi
                ! rather than doing eta = true*Bnu + scat*Jnu
                ! we can also do eta = (true + scat)*Bnu + scat*(Jnu - Bnu)
                Jnu = Jnu - Bnu
              endif
              eta = chi * Bnu
            end associate
          end do
        end do
      end do
    end do

    call comp%scattering(self%sub_sp, self%back(:,:,:,:,2), self%Jnu(:,:,:,:,i_set), include_absmet=include_absmet)
    call trace%end(itimer)
  END SUBROUTINE

  !===============================================================================
  SUBROUTINE generate_opacity_table (self, sp, include_absmet, include_llist, build_deriv)
    class(rt_patch_t), target :: self
    class(composition_type), pointer :: comp
    type(spectrum_type) :: sp
    real(4), dimension(:,:,:,:), pointer :: all_chi
    real(4), dimension(:,:,:,:), allocatable :: all_eta, all_beta
    logical, intent(in) :: include_absmet, include_llist, build_deriv
    real(4) :: over_chi500
    integer :: i,j,k, ntemp, nne, nnu
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%generate_opacity_table', itimer=itimer)
    comp => self%m3d_parent%tab_comp
    nnu = sp%nnu
    ntemp = size(comp%temp,1)
    nne   = size(comp%temp,2)
    allocate(self%opac%val(nnu, ntemp, nne), source=0.)
    allocate(self%opac%lam, source=sp%lam)
    self%opac%nnu = nnu
    self%opac%ntemp = ntemp
    self%opac%nne = nne
    call mem_register(self%opac%val, label='opac_table')

    all_chi(1:nnu, 1:ntemp, 1:nne, 1:1) => self%opac%val
    call comp%absorption(sp, all_chi, include_absmet=include_absmet, sorted=self%sp%sorted)

    if (include_llist) then
      ! absorption returns the actual chi, while add_linelist_opacities2
      ! adds the working versions of chi and eta
      allocate(all_eta, all_beta, mold=all_chi)
      do k=1, nne
        do j=1, ntemp
          associate(chi=>all_chi(:,j,k,1), eta=>all_eta(:,j,k,1), beta=>all_beta(:,j,k,1), &
                    xhat=>all_chi(:,j,k,1), ehat=>all_eta(:,j,k,1))
            beta = exp(-h_k * sp%nu / comp%temp(j,k,1))
            ehat = chi / (1 - beta)
            xhat = chi + ehat*beta
          end associate
        end do
      end do
      call comp%add_linelist_opacities2(sp, all_chi, all_eta, dlam)
      do k=1, nne
        do j=1, ntemp
          associate(chi=>all_chi(:,j,k,1), eta=>all_eta(:,j,k,1), beta=>all_beta(:,j,k,1), &
                    xhat=>all_chi(:,j,k,1), ehat=>all_eta(:,j,k,1))
            eta = ehat * beta
            chi = xhat - eta
          end associate
        end do
      end do
    endif

    do k=1, nne
      do j=1, ntemp
        over_chi500 = 1 / comp%chi500(j,k,1)
        self%opac%val(:,j,k) = self%opac%val(:,j,k) * over_chi500
      end do
    end do

    ! On the sampling path, capture the PCA (lnT, lnNe) axis transform built in
    ! init_tabulated_comp and precompute the monotone node derivatives along every
    ! dimension. The EOS path (build_deriv=.false.) uses a temp-rho grid and only
    ! needs val/lam.
    if (build_deriv) then
      self%opac%pca = comp%pca
      call self%opac%build_derivatives()
    endif
    self%opac%initialised = .true.
    call trace%end(itimer)
  END SUBROUTINE generate_opacity_table

  !===============================================================================
  ! Build the separate sampled LTE line-list table (sample_linelist). tab_comp is
  ! masked to the non-departure species (see mask_departure_species), so the line
  ! opacity added here is pure LTE. In LTE the two working forms coincide
  ! (xhat_line = ehat_line = chi_line/(1-beta)); we store that single quantity,
  ! normalised by chi500, and at interpolation time add it to both work slots.
  ! The pure line delta is isolated by subtracting the continuum working xhat that
  ! add_linelist needs in place for its opacity-threshold reference.
  !===============================================================================
  SUBROUTINE generate_ll_opacity_table (self, sp, include_absmet)
    class(rt_patch_t), target :: self
    class(composition_type), pointer :: comp
    type(spectrum_type) :: sp, psp
    logical, intent(in) :: include_absmet
    real(4), dimension(:,:,:,:), pointer :: all_chi
    real(4), dimension(:,:,:,:), allocatable :: all_eta, all_beta, xhat_cont
    real(4), dimension(:,:,:,:), allocatable :: chi_scat
    real(4), dimension(:,:,:,:,:), allocatable :: one_m_eps
    real(8), dimension(:), allocatable :: plam
    real(8) :: vmax, dblue, dred
    real(4) :: over_chi500
    integer :: i, j, k, ntemp, nne, nnu, nnu0, npb, npr
    integer, parameter :: pad_guard = 2
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%generate_ll_opacity_table', itimer=itimer)
    comp => self%m3d_parent%tab_comp
    ntemp = size(comp%temp,1)
    nne   = size(comp%temp,2)

    ! ---- Doppler padding of the table wavelength grid ------------------------
    ! At interpolation each cell queries lam*(1+vproj), vproj = vec.v (|vec|=1) in
    ! units of c, so |vproj| <= vmax = max cell speed / c. The chunk output grid
    ! spans [lam(fs), lam(fe)]; without padding a shifted query at the chunk edge
    ! falls off the table (searchsorted clamps) and any line that Doppler-shifts
    ! into the window there is missed. Extend the grid geometrically at its edge
    ! log-spacing by enough nodes to cover +/-vmax (plus a guard for line wings),
    ! then build the table on the padded grid. vmax comes from the un-rotated
    ! composition; rotation preserves the speed, so it bounds every direction.
    vmax = 0
    if (associated(self%m3d_parent%composition)) then
      associate(bc => self%m3d_parent%composition)
        if (allocated(bc%vx)) vmax = sqrt(real(maxval(bc%vx**2 + bc%vy**2 + bc%vz**2), 8))
      end associate
    endif
    nnu0 = sp%nnu
    if (vmax > 0 .and. nnu0 >= 2) then
      dblue = log(sp%lam(2))    - log(sp%lam(1))
      dred  = log(sp%lam(nnu0)) - log(sp%lam(nnu0-1))
      npb = ceiling(vmax/abs(dblue)) + pad_guard
      npr = ceiling(vmax/abs(dred))  + pad_guard
      allocate(plam(nnu0 + npb + npr))
      do i=1, npb
        plam(i) = sp%lam(1) * exp(-(npb-i+1)*dblue)     ! blue extension (ascending lam)
      end do
      plam(npb+1:npb+nnu0) = sp%lam
      do i=1, npr
        plam(npb+nnu0+i) = sp%lam(nnu0) * exp(i*dred)    ! red extension
      end do
      call psp%reset(lam=plam)
    else
      call psp%reset(lam=sp%lam)
    endif
    psp%sorted = sp%sorted
    nnu = psp%nnu

    allocate(self%ll_opac%val(nnu, ntemp, nne), source=0.)
    allocate(self%ll_opac%lam, source=psp%lam)
    self%ll_opac%nnu = nnu
    self%ll_opac%ntemp = ntemp
    self%ll_opac%nne = nne
    call mem_register(self%ll_opac%val, label='ll_opac_table')

    ! continuum reference (needed for add_linelist's opacity threshold), converted
    ! into the working xhat/ehat units the line routines accumulate onto
    all_chi(1:nnu, 1:ntemp, 1:nne, 1:1) => self%ll_opac%val
    call comp%absorption(psp, all_chi, include_absmet=include_absmet, sorted=self%sp%sorted)
    allocate(all_eta, all_beta, xhat_cont, mold=all_chi)
    do k=1, nne
      do j=1, ntemp
        associate(chi=>all_chi(:,j,k,1), eta=>all_eta(:,j,k,1), beta=>all_beta(:,j,k,1), &
                  xhat=>all_chi(:,j,k,1), ehat=>all_eta(:,j,k,1))
          beta = exp(-h_k * psp%nu / comp%temp(j,k,1))
          ehat = chi / (1 - beta)
          xhat = chi + ehat*beta
        end associate
      end do
    end do
    xhat_cont = all_chi     ! snapshot of the continuum working xhat

    ! LTE line list (no departures/vproj); tab_comp is masked to non-departure
    ! species, so only LTE species contribute here
    if (line_scattering) then
      ! scattering split: precompute the per-stage scattering fraction 1-eps at
      ! the node (T,ne) on the (padded) table grid, accumulate the
      ! (1-eps)-weighted line opacity into chi_scat, and store it as the table's
      ! second component (normalised the same way as the total)
      allocate(one_m_eps(nnu, ntemp, nne, 1, 2))
      allocate(chi_scat(nnu, ntemp, nne, 1), source=0.)
      call compute_epsilon_stage(1, psp%nu, comp%temp, comp%ne, one_m_eps(:,:,:,:,1))
      call compute_epsilon_stage(2, psp%nu, comp%temp, comp%ne, one_m_eps(:,:,:,:,2))
      call comp%add_linelist_opacities2(psp, all_chi, all_eta, dlam, &
        one_m_eps=one_m_eps, chi_scat=chi_scat)
      allocate(self%ll_opac%scat(nnu, ntemp, nne), source=0.)
      call mem_register(self%ll_opac%scat, label='ll_opac_table')
    else
      call comp%add_linelist_opacities2(psp, all_chi, all_eta, dlam)
    endif

    ! isolate the pure line xhat delta and normalise by chi500 (undone per cell in
    ! sample_interpolate_ll). In LTE the ehat delta equals the xhat delta.
    do k=1, nne
      do j=1, ntemp
        over_chi500 = 1 / comp%chi500(j,k,1)
        self%ll_opac%val(:,j,k) = (all_chi(:,j,k,1) - xhat_cont(:,j,k,1)) * over_chi500
        if (allocated(chi_scat)) &
          self%ll_opac%scat(:,j,k) = chi_scat(:,j,k,1) * over_chi500
      end do
    end do

    self%ll_opac%pca = comp%pca
    call self%ll_opac%build_derivatives()
    self%ll_opac%initialised = .true.
    call trace%end(itimer)
  END SUBROUTINE generate_ll_opacity_table

  !===============================================================================
  SUBROUTINE save_opac_table (self)
    class(rt_patch_t), target :: self
    class(composition_type), pointer :: comp
    real(4), dimension(:,:,:), allocatable :: eos
    real(4), dimension(:,:,:,:), pointer :: view
    real(4) :: over_rho
    integer :: i,j,k, ntemp, nrho, nnu
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%save_opac_table', itimer=itimer)
    comp => self%m3d_parent%tab_comp
    ntemp = size(self%opac%val, 2)
    nrho = size(self%opac%val, 3)
    nnu = size(self%opac%lam)

    ! save absorption part first (which has been normalised by chi500)
    do k=1, nrho
      do j=1, ntemp
        over_rho = comp%chi500(j,k,1) / comp%rho(j,k,1)
        self%opac%val(:,j,k) = self%opac%val(:,j,k) * over_rho
      end do
    end do
    call self%save3(self%opac%val, 'opa')

    ! then save scattering part
    self%opac%val = 0.0
    view(1:nnu, 1:ntemp, 1:nrho, 1:1) => self%opac%val
    call comp%electron_scattering(view)
    call comp%rayleigh_scattering(self%opac%lam, view)
    do k=1, nrho
      do j=1, ntemp
        over_rho = 1 / comp%rho(j,k,1)
        self%opac%val(:,j,k) = self%opac%val(:,j,k) * over_rho
      end do
    end do
    call self%save3(self%opac%val, 'scat')
    call trace%end(itimer)
  END SUBROUTINE save_opac_table

  !===============================================================================
  SUBROUTINE reset_array (self, array)
    class(rt_patch_t), target :: self
    real(4), dimension(:,:,:,:,:), intent(inout) :: array
    integer :: i
    integer :: x, y, z
    !-----------------------------------------------------------------------------
    do i=lbound(array, 5), ubound(array,5)
      do z=lbound(array, 4), ubound(array,4)
        do y=1, self%ny
          do x=1, self%nx
            array(:,x,y,z,i) = 0.
          end do
        end do
      end do
    end do
  END SUBROUTINE reset_array

  !===============================================================================
  SUBROUTINE add_bound_free_opacities (self, atom, trans, adata)
    class(rt_patch_t), target :: self
    class(atom_type) :: atom
    class(atom_trans_t) :: trans
    class(atom_data_t), target :: adata
    class(cont_type), pointer :: cont
    real(4) :: ni, nj, tmp, alpha
    integer :: i, fs, fe
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%add_bound_free_opacities', itimer=itimer)
    
    ! Summing contributions from all bound-free transitions
    do i=1, size(trans%cont)
      cont => trans%cont(i)
      if (.not. cont%active) cycle
      fs = cont%iblue
      fe = cont%ired

      if (allocated(adata%b)) then
        do z=1, self%nz
          do y=1, self%ny
            do x=1, self%nx
              ni  = adata%nstar(x,y,z,cont%i) * adata%b(x,y,z,cont%i)
              tmp = adata%nstar(x,y,z,cont%i) * adata%b(x,y,z,cont%j)
              do f=fs, fe
                alpha = cont%alpha(f)
                self%back(f,x,y,z,1) = self%back(f,x,y,z,1) + alpha * ni  ! chi
                self%back(f,x,y,z,2) = self%back(f,x,y,z,2) + alpha * tmp ! eta
              end do
            end do
          end do
        end do
      else ! b = 1
        do z=1, self%nz
          do y=1, self%ny
            do x=1, self%nx
              ni  = adata%nstar(x,y,z,cont%i)
              tmp = adata%nstar(x,y,z,cont%i)
              do f=fs, fe
                alpha = cont%alpha(f)
                self%back(f,x,y,z,1) = self%back(f,x,y,z,1) + alpha * ni  ! chi
                self%back(f,x,y,z,2) = self%back(f,x,y,z,2) + alpha * tmp ! eta
              end do
            end do
          end do
        end do
      endif
    end do
    call trace%end(itimer)
  END SUBROUTINE add_bound_free_opacities

  !===============================================================================
  SUBROUTINE bound_bound_damping (self, comp, atom, trans, adata)
    class(rt_patch_t), target :: self
    class(composition_type) :: comp
    class(atom_type) :: atom
    class(atom_trans_t) :: trans
    class(atom_data_t) :: adata
    class(line_type), pointer :: line
    real(4), dimension(:,:,:), allocatable :: tothi, vbar2, ne, temp, nh6
    real(4), dimension(:,:,:), allocatable :: stark_pe, unsold
    integer :: i
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin (trim(self%kind)//'%bound_bound_damping', itimer=itimer)
    allocate(temp, ne, tothi, vbar2, stark_pe, unsold, nh6, mold=adata%inv_df)
    ne = real(comp%ne)
    temp = real(comp%temp)
    tothi = real(comp%tothi) ! total neutral hydrogen density
    nh6 = real(comp%elements(1)%ions(1)%levels(6)%n) ! neutral hydrogen density in level 6
    stark_pe = ne * temp**one_sixth
    vbar2 = atom%compute_vbar2(temp)
    unsold = vbar2**0.3 * tothi

    ! computing damping coefficients for all bb-transitions
    do i=1, size(trans%line)
      line => trans%line(i)
      if (.not. line%active) cycle
      if (.not. line%hydrogen) then
        ! damping for hydrogen lines is computed in hprofl_mod
        call line%damping(trans%adamp(:,:,:,i), ne, tothi, vbar2, stark_pe, unsold)
        trans%adamp(:,:,:,i) = trans%adamp(:,:,:,i) * adata%inv_df
      endif
    end do
    call trace%end(itimer)
  END SUBROUTINE bound_bound_damping

  !===============================================================================
  SUBROUTINE add_bound_bound_opacities (self, comp, atom, trans, adata, wmu, vproj, LTE, use_ll)
    class(rt_patch_t), target :: self
    class(composition_type) :: comp
    class(atom_type) :: atom
    class(atom_trans_t) :: trans
    class(atom_data_t), target :: adata
    class(line_type), pointer :: line
    real(4), contiguous, dimension(:,:,:) :: vproj
    real(4), dimension(:,:,:), allocatable :: wphi, temp, ne, fo, vbar2h, tothi, t3nhe
    real(4), dimension(:,:,:), allocatable :: e_kT4, ne4, tothi4, sp4, idf3, vb3
    real(4), dimension(:,:,:,:), allocatable :: lnbase, dummy_dep
    real(4) :: eps_dummy(0,0,0,0), scat_dummy(0,0,0,0)
    real(4), dimension(:), allocatable :: vphi
    real(4) :: wmu, idf, vpr, a, idf_4pi, ni, nj, Bji_nj, Bij_ni, widf, wp, amin, amax
    real(8) :: q
    real(4) :: lower_off, upper_off, ekt, wj, chi_fac, eta_fac, bi, bj, Rij
    integer :: i, ih, fs, fe
    logical :: LTE, do_LTE, get_wphi, use_ll
    integer, save :: itimer=0, itimer1=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin (trim(self%kind)//'%add_bound_bound_opacities', itimer=itimer)
    ! Copying continuum opacities to angle dependent opacities

    allocate(vphi(self%fs:self%fe))
    allocate(wphi, mold=vproj)
    get_wphi = self%iter == 0 .and. .not. self%do_out_sp
    do_LTE = LTE .or. .not. allocated(adata%b)

    if (any(trans%line%hydrogen)) then
      allocate(temp, ne, fo, vbar2h, tothi, t3nhe, mold=adata%inv_df)
      temp = comp%temp
      ne = comp%ne
      fo = comp%fo
      vbar2h = comp%vbar2h
      tothi = comp%tothi
      t3nhe = comp%t3nhe
    endif

    do i=1, size(trans%line)
      line => trans%line(i)
      ! non-hydrogen lines are handled by the SoA kernel below when use_ll
      if (use_ll .and. .not. line%hydrogen) cycle
      fs = line%iblue
      fe = line%ired

      lower_off = line%lower_ev - atom%ev(line%i)
      upper_off = line%upper_ev - atom%ev(line%j)

      ! Adding bound-bound opacities
      do z=1, self%nz
        do y=1, self%ny
          do x=1, self%nx
            idf = adata%inv_df(x,y,z)
            ekt = comp%e_kT(x,y,z)
            ni = adata%nstar(x,y,z,line%i) ! these are LTE populations
            nj = adata%nstar(x,y,z,line%j) ! these are LTE populations

            if (do_LTE) then
              bi = 1.0
              bj = 1.0
            else
              bi = adata%b(x,y,z,line%i)
              bj = adata%b(x,y,z,line%j)
            endif

            if (lower_off /= 0) bi = bi * exp(-lower_off*ekt)
            if (upper_off /= 0) bj = bj * exp(-upper_off*ekt)

            associate(phi => vphi(fs:fe))
              if (line%hydrogen) then
                do f=fs, fe
                  q = line%q(f) - vproj(x,y,z)
                  vphi(f) = hprofile(line, q, temp(x,y,z), ne(x,y,z), fo(x,y,z), &
                    vbar2h(x,y,z), tothi(x,y,z), t3nhe(x,y,z), idf)
                end do
                if (H_occupation_probability) then
                  wj = comp%elements(1)%ions(1)%levels(line%j)%w(x,y,z)
                  bi = bi * wj
                  bj = bj * wj
                endif
              elseif (line%io > 1) then
                phi = 0
                do ih=1, line%io
                  a = trans%adamp(x,y,z,i) + idf*line%da(ih)
                  vpr = vproj(x,y,z) + line%dq(ih)
                  phi = phi + voigt_profile(a, line%q(fs:fe), vpr, idf) * line%wt(ih)
                end do
              else
                vpr = vproj(x,y,z)
                a = trans%adamp(x,y,z,i)
                phi = voigt_profile(a, line%q(fs:fe), vpr, idf)
              endif

              Rij = ni * line%Bij
              eta_fac = Rij * bj
              chi_fac = Rij * bi

              associate(xhat=>self%work(fs:fe,x,y,z,1), ehat=>self%work(fs:fe,x,y,z,2))
                xhat = xhat + phi * chi_fac
                ehat = ehat + phi * eta_fac
              end associate

              if (get_wphi) then
                wphi(x,y,z) = sum(phi * line%wq(fs:fe))
              endif
            end associate
          end do
        end do
      end do

      if (get_wphi) then
        call lock_mem(self%m3d_parent, 'wphi')
        associate(sum_wphi => adata%wphi(:,:,:,line%ntrans))
          sum_wphi = sum_wphi + wmu * wphi
        end associate
        call unlock_mem(self%m3d_parent, 'wphi')
      endif
    end do

    ! SoA kernel for the non-hydrogen lines: exp(arg) populations from the
    ! composition ion densities, damping recomputed per cell, LTE wing cutoff
    if (use_ll .and. trans%ll%nline > 0) then
      call get_ll_env(comp, atom, trans, adata, e_kT4, ne4, tothi4, sp4, idf3, vb3, lnbase)
      associate(xhat=>self%work(:,:,:,1:self%nz,1), ehat=>self%work(:,:,:,1:self%nz,2), &
                cont=>self%back(:,:,:,:,1))
        if (do_LTE) then
          allocate(dummy_dep(1,1,1,1), source=0.)
          if (self%sub_sp%equidistant) then
            call ll_add_opacity_equi(trans%ll, self%llam, xhat, ehat, cont, self%logref, &
              vproj, idf3, vb3, e_kT4, lnbase, ne4, tothi4, sp4, dummy_dep, .false., &
              eps_dummy, scat_dummy, .false.)
          else
            call ll_add_opacity(trans%ll, self%llam, xhat, ehat, cont, self%logref, &
              vproj, idf3, vb3, e_kT4, lnbase, ne4, tothi4, sp4, dummy_dep, .false., &
              eps_dummy, scat_dummy, .false.)
          endif
        else
          call self%fill_b4(adata, atom%nlevel)
          if (self%sub_sp%equidistant) then
            call ll_add_opacity_equi(trans%ll, self%llam, xhat, ehat, cont, self%logref, &
              vproj, idf3, vb3, e_kT4, lnbase, ne4, tothi4, sp4, &
              self%b4(:,:,:,1:atom%nlevel), .true., eps_dummy, scat_dummy, .false.)
          else
            call ll_add_opacity(trans%ll, self%llam, xhat, ehat, cont, self%logref, &
              vproj, idf3, vb3, e_kT4, lnbase, ne4, tothi4, sp4, &
              self%b4(:,:,:,1:atom%nlevel), .true., eps_dummy, scat_dummy, .false.)
          endif
        endif
      end associate
    endif
    call trace%end(itimer)
  END SUBROUTINE add_bound_bound_opacities

  !===============================================================================
  ! Gather the per-cell inputs of the SoA line kernels (linelist_mod) for one atom:
  ! thermodynamic fields, Doppler arrays (single isotope, index 0) and the per-ion-
  ! stage population base
  !   lnbase(:,:,:,s) = log(n_s) - log_part_s + dxi*(s-1)*e_kT + log(inv_df)
  ! from the composition particle (element ion stage s, or the molecule itself) --
  ! the same formalism as the bulk linelist path in add_linelist_opacities2.
  !-------------------------------------------------------------------------------
  SUBROUTINE get_ll_env (comp, atom, trans, adata, e_kT, ne, tothi, stark_pe, &
                         idf3, vb3, lnbase)
    class(composition_type), target :: comp
    class(atom_type) :: atom
    class(atom_trans_t) :: trans
    class(atom_data_t) :: adata
    real(4), dimension(:,:,:),   allocatable, intent(out) :: e_kT, ne, tothi, stark_pe
    real(4), dimension(:,:,:),   allocatable, intent(out) :: idf3, vb3
    real(4), dimension(:,:,:,:), allocatable, intent(out) :: lnbase
    class(element_type), pointer :: element
    class(molec_type), pointer :: mol
    integer :: m(3), np, s
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%get_ll_env', itimer=itimer)
    m = shape(adata%inv_df)
    allocate(e_kT(m(1),m(2),m(3)))
    allocate(ne, tothi, stark_pe, mold=e_kT)
    e_kT     = comp%e_kT
    ne       = real(comp%ne)
    tothi    = real(comp%tothi)
    stark_pe = real(comp%stark_pe)

    allocate(idf3(m(1),m(2),m(3)), vb3(m(1),m(2),m(3)))
    idf3 = adata%inv_df
    vb3  = atom%compute_vbar2(real(comp%temp, 4))

    np = maxval(trans%ll%ipart(1:trans%ll%nline))
    allocate(lnbase(m(1),m(2),m(3),np))
    lnbase = -huge(1.)   ! stages without an EOS particle contribute no opacity

    call comp%element_by_id(atom%id, element, verbose=.false.)
    if (associated(element)) then
      do s=1, np
        if (s < lbound(element%ions,1) .or. s > ubound(element%ions,1)) cycle
        associate(ion => element%ions(s))
          if (.not. allocated(ion%n)) cycle
          lnbase(:,:,:,s) = log(ion%n) - ion%part &
                          + real(comp%dxi) * (s-1) * e_kT + log(idf3)
        end associate
      end do
    else
      call comp%molecule_by_id(atom%id, mol, verbose=.false.)
      if (associated(mol)) then
        if (allocated(mol%n)) then
          lnbase(:,:,:,1) = log(mol%n) - mol%part + log(idf3)
        endif
      endif
    endif
    call trace%end(itimer)
  END SUBROUTINE get_ll_env

  !===============================================================================
  ! real(4) copy of the departure coefficients for the SoA kernels; cached per
  ! adata within one execute (b is constant during an execute)
  !-------------------------------------------------------------------------------
  SUBROUTINE fill_b4 (self, adata, nlevel)
    class(rt_patch_t) :: self
    class(atom_data_t), target :: adata
    integer, intent(in) :: nlevel
    integer :: l
    !-----------------------------------------------------------------------------
    if (associated(self%b4_ptr, adata)) then
      if (size(self%b4,4) >= nlevel) return
    endif
    if (allocated(self%b4)) then
      if (size(self%b4,4) < nlevel) deallocate(self%b4)
    endif
    if (.not. allocated(self%b4)) then
      allocate(self%b4(self%nx,self%ny,self%nz,nlevel))
    endif
    do l=1, nlevel
      self%b4(:,:,:,l) = adata%b(:,:,:,l)
    end do
    self%b4_ptr => adata
  END SUBROUTINE fill_b4

  !===============================================================================
  ! Push the used/strength diagnostics of the per-patch SoA lists back to the
  ! shared source lists (self%atoms(k)%trans%ll), which print_linelist exports.
  !-------------------------------------------------------------------------------
  SUBROUTINE reduce_ll_diag (self)
    class(rt_patch_t) :: self
    logical, allocatable :: gused(:)
    real(4), allocatable :: gstr(:)
    integer :: k, l, gl
    !-----------------------------------------------------------------------------
    do k=1, self%natom
      if (.not. self%use_ll(k)) cycle
      associate(pll => self%trans(k)%ll, sll => self%atoms(k)%trans%ll, &
                lines => self%trans(k)%line)
        if (pll%nline == 0) cycle
        if (.not. any(pll%used(1:pll%nline))) cycle
        ! collect per (global) transition, then mark all its source entries
        ! (multi-component entries share the same transition index)
        allocate(gused(size(self%atoms(k)%trans%line)), source=.false.)
        allocate(gstr (size(self%atoms(k)%trans%line)), source=0.)
        do l=1, pll%nline
          if (.not. pll%used(l)) cycle
          gl = lines(pll%itrans(l))%ntrans
          gused(gl) = .true.
          gstr(gl)  = max(gstr(gl), pll%strength(l))
        end do
        call lock_mem(self%m3d_parent, 'll_used')
        do l=1, sll%nline
          gl = sll%itrans(l)
          if (gused(gl)) then
            sll%used(l)     = .true.
            sll%strength(l) = max(sll%strength(l), gstr(gl))
          endif
        end do
        call unlock_mem(self%m3d_parent, 'll_used')
        deallocate(gused, gstr)
      end associate
    end do
  END SUBROUTINE reduce_ll_diag

  !===============================================================================
  SUBROUTINE bound_bound_rates (self, comp, atom, trans, adata, wmu, vproj, use_ll)
    class(rt_patch_t), target :: self
    class(composition_type) :: comp
    class(atom_type) :: atom
    class(atom_trans_t) :: trans
    class(atom_data_t), target :: adata
    class(line_type), pointer :: line
    real(4), contiguous, dimension(:,:,:) :: vproj
    real(4), dimension(:,:,:), allocatable :: temp, ne, fo, vbar2h, tothi, t3nhe
    real(4), dimension(:,:,:), allocatable :: e_kT4, ne4, tothi4, sp4, idf3, vb3
    real(4), dimension(:,:,:,:), allocatable :: lnbase, wphi_buf
    real(4), dimension(:), allocatable :: wmunu
    real(4) :: wmu, idf, vpr, a, idf_4pi, ni, nj, chi_fac, phi, vel, vgt
    real(4) :: vvv, a2, sum_I, sum_P, sum_A
    real(4) :: g_ij, g_ji, q
    real(4) :: lower_off, upper_off, ekt
    integer :: i, l, ih, fs, fe
    logical :: use_ll, get_wphi
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin (trim(self%kind)//'%bound_bound_rates', itimer=itimer)
    allocate(wmunu(self%fs:self%fe))
    ! Computing bound-bound radiative rates

    if (any(trans%line%hydrogen)) then
      allocate(temp, ne, fo, vbar2h, tothi, t3nhe, mold=adata%inv_df)
      temp = comp%temp
      ne = comp%ne
      fo = comp%fo
      vbar2h = comp%vbar2h
      tothi = comp%tothi
      t3nhe = comp%t3nhe
    endif

    do i=1, size(trans%line)
      line => trans%line(i)
      ! non-hydrogen lines are handled by the SoA rate kernel below when use_ll
      if (use_ll .and. .not. line%hydrogen) cycle
      fs = line%iblue
      fe = line%ired
      wmunu(fs:fe) = wmu * line%wnu(fs:fe) * line%Aij
      lower_off = line%lower_ev - atom%ev(line%i)
      upper_off = line%upper_ev - atom%ev(line%j)

      do z=1,self%nz
        do y=1,self%ny
          do x=1,self%nx
            ni = adata%b(x,y,z,line%i) * adata%nstar(x,y,z,line%i)
            nj = adata%b(x,y,z,line%j) * adata%nstar(x,y,z,line%j)
            idf = adata%inv_df(x,y,z)
            vpr = vproj(x,y,z)
            a = trans%adamp(x,y,z,i)
            ekt = comp%e_kT(x,y,z)
            if (lower_off /= 0) ni = ni * exp(-lower_off*ekt)
            if (upper_off /= 0) nj = nj * exp(-upper_off*ekt)

            sum_A = 0.0
            sum_I = 0.0
            sum_P = 0.0
            if (line%hydrogen) then
              do f=fs, fe
                associate(Inu=>self%work(f,x,y,z,3), Psi=>self%Psi(f,x,y,z))
                  q = line%q(f) - vpr
                  phi = hprofile(line, dble(q), temp(x,y,z), ne(x,y,z), fo(x,y,z), &
                    vbar2h(x,y,z), tothi(x,y,z), t3nhe(x,y,z), idf)
                  vvv = phi * wmunu(f)
                  sum_A = sum_A + vvv
                  sum_I = sum_I + vvv * Inu
                  sum_P = sum_P + vvv * phi * Psi
                end associate
              end do
            elseif (line%io > 1) then
              do f=fs, fe
                associate(Inu=>self%work(f,x,y,z,3), Psi=>self%Psi(f,x,y,z))
                  phi = 0
                  do ih=1, line%io
                    q = line%q(f) - line%dq(ih)
                    a2 = a + idf*line%da(ih)
                    vgt = voigt_profile(a2, q, vpr, idf)
                    phi = phi + vgt*line%wt(ih)
                  end do
                  vvv = phi * wmunu(f)
                  sum_A = sum_A + vvv
                  sum_I = sum_I + vvv * Inu
                  sum_P = sum_P + vvv * phi * Psi
                end associate
              end do
            else
              do f=fs, fe
                associate(Inu=>self%work(f,x,y,z,3), Psi=>self%Psi(f,x,y,z))
                  phi = voigt_profile(a, line%q(f), vpr, idf)
                  vvv = phi * wmunu(f)
                  sum_A = sum_A + vvv
                  sum_I = sum_I + vvv * Inu
                  sum_P = sum_P + vvv * phi * Psi
                end associate
              end do
            endif

            g_ij = sum_I
            g_ji = sum_I + sum_A

            if (preconditioning) then
              g_ij = g_ij - sum_P * line%Bji * nj
              g_ji = g_ji - sum_P * line%Bij * ni
            endif

            trans%gamma_ij(x,y,z,i) = trans%gamma_ij(x,y,z,i) + g_ij
            trans%gamma_ji(x,y,z,i) = trans%gamma_ji(x,y,z,i) + g_ji
          end do
        end do
      end do
    end do

    ! SoA rate kernel for the non-hydrogen lines; also accumulates the profile
    ! normalization wphi (moved here from the opacity routine) on the first iteration
    if (use_ll .and. trans%ll%nline > 0) then
      call get_ll_env(comp, atom, trans, adata, e_kT4, ne4, tothi4, sp4, idf3, vb3, lnbase)
      call self%fill_b4(adata, atom%nlevel)
      get_wphi = self%iter == 0 .and. .not. self%do_out_sp .and. fixed_Jnu <= 0
      associate(Inu=>self%work(:,:,:,1:self%nz,3), cont=>self%back(:,:,:,:,1))
        if (get_wphi) then
          allocate(wphi_buf(self%nx,self%ny,self%nz,size(trans%line)), source=0.)
          call ll_add_rates(trans%ll, self%llam, Inu, self%Psi, cont, self%logref, &
            vproj, idf3, vb3, e_kT4, lnbase, ne4, tothi4, sp4, &
            self%b4(:,:,:,1:atom%nlevel), allocated(adata%b), wmu, preconditioning, &
            trans%gamma_ij, trans%gamma_ji, wphi_buf, .true.)
          call lock_mem(self%m3d_parent, 'wphi')
          do i=1, size(trans%line)
            if (trans%line(i)%hydrogen) cycle
            associate(sum_wphi => adata%wphi(:,:,:,trans%line(i)%ntrans))
              sum_wphi = sum_wphi + wphi_buf(:,:,:,i)
            end associate
          end do
          call unlock_mem(self%m3d_parent, 'wphi')
        else
          call ll_add_rates(trans%ll, self%llam, Inu, self%Psi, cont, self%logref, &
            vproj, idf3, vb3, e_kT4, lnbase, ne4, tothi4, sp4, &
            self%b4(:,:,:,1:atom%nlevel), allocated(adata%b), wmu, preconditioning, &
            trans%gamma_ij, trans%gamma_ji, get_wphi=.false.)
        endif
      end associate
    endif
    call trace%end(itimer)
  END SUBROUTINE bound_bound_rates

  !===============================================================================
  SUBROUTINE bound_free_rates (self, atom, trans, adata)
    class(rt_patch_t), target :: self
    class(atom_type) :: atom
    class(atom_trans_t) :: trans
    class(atom_data_t), target :: adata
    class(cont_type), pointer :: cont
    real(4), dimension(:), allocatable :: w_alpha, wmunu
    real(4) :: ni, nj, tmp, eta_fac, q, phi
    real(4) :: g_ij, g_ji, vij, vji, sum_abs, sum_ems, sum_pre
    integer :: i, j, fs, fe
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin (trim(self%kind)//'%bound_free_rates', itimer=itimer)
    allocate(wmunu(self%fs:self%fe))
    allocate(w_alpha(self%fs:self%fe))
    ! Computing bound-free radiative rates based on half integrated mean intensity
    do j=1, size(trans%cont)
      cont => trans%cont(j)
      if (.not. cont%active) cycle
      i = j + size(trans%line)
      fs = cont%iblue
      fe = cont%ired
      wmunu(fs:fe) = 0.5 * cont%wnu(fs:fe)
      w_alpha(fs:fe) = self%twohnu3_c2(fs:fe) * cont%alpha(fs:fe) * wmunu(fs:fe)

      do z=1,self%nz
        do y=1,self%ny
          do x=1,self%nx
            sum_abs = 0.0
            sum_ems = 0.0
            sum_pre = 0.0

            do f=fs, fe
              associate(Jnu=>self%Jnu(f,x,y,z,1), beta=>self%back(f,x,y,z,3), Psi=>self%PsiJ(f,x,y,z))
                vij = w_alpha(f)
                vji = vij * beta
                phi = Psi * cont%alpha(f)
                sum_abs = sum_abs + vij * Jnu
                sum_ems = sum_ems + vji * (Jnu + 1.0)
                sum_pre = sum_pre + vji * phi
              end associate
            end do

            g_ij = sum_abs
            g_ji = sum_ems

            if (preconditioning) then
              g_ij = g_ij - sum_pre * adata%nstar(x,y,z,cont%i) * adata%b(x,y,z,cont%j)
              g_ji = g_ji - sum_pre * adata%nstar(x,y,z,cont%i) * adata%b(x,y,z,cont%i)
            endif

            trans%gamma_ij(x,y,z,i) = trans%gamma_ij(x,y,z,i) + g_ij
            trans%gamma_ji(x,y,z,i) = trans%gamma_ji(x,y,z,i) + g_ji
          end do
        end do
      end do
    end do
    call trace%end(itimer)
  END SUBROUTINE bound_free_rates

  !===============================================================================
  SUBROUTINE set_upstream_boundary (self, i_mu, i_set, ng)
    class(rt_patch_t) :: self
    integer :: i_mu, i_set, ng, i
    real(4) :: xxold(3), xxnew(ng)
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%set_upstream_boundary', itimer=itimer)
    z = merge(lbound(self%work,4), ubound(self%work,4), self%dir > 0)

    if (self%boundaries%is_set(bits%zu)) then
      self%work(:,:,:,z,:) = 0.
    elseif (self%boundaries%is_set(bits%zl)) then
      xxold = [(i*self%ds(3), i=1, 3)]
      xxnew = [(i*self%ds(3), i=1-ng, 0)]
      do i=1, 2
        do y=1,  self%ny
          do x=1,  self%nx
            do f=self%fs, self%fe
              call pchip_mono(xxold, self%work(f,x,y,z+1:z+3,i), &
                xxnew, self%work(f,x,y,z:z,i), sorted=.true.)
              self%work(f,x,y,z,i) = max(self%work(f,x,y,z,i), 0.)
            end do
          end do
        end do
      end do
      do y=1,  self%ny
        do x=1,  self%nx
          do f=self%fs, self%fe
            self%work(f,x,y,z,3) = self%work(f,x,y,z,2)
          end do
        end do
      end do
    else
      do y=1,  self%ny
        do x=1,  self%nx
          do f=self%fs, self%fe
            self%work(f,x,y,z,1) = self%mem(f,x,y,i_mu,i_set,1)
            self%work(f,x,y,z,2) = self%mem(f,x,y,i_mu,i_set,2)
            self%work(f,x,y,z,3) = self%mem(f,x,y,i_mu,i_set,3)
          end do
        end do
      end do
      !-------------------------------------------------------------------------
      ! The incoming boundary comes straight from the upstream patch, via the
      ! buffer handed over in copy_to().  If it is already non-finite here, the
      ! corruption did not happen in this patch, and the upstream id says where
      ! to look next.  If it is clean, but mem is non-finite in post_update, the
      ! garbage was produced locally (opacities, Jnu, or the transfer solve).
      !-------------------------------------------------------------------------
      if (check_NaN) then
        do y=1,  self%ny
          do x=1,  self%nx
            do f=self%fs, self%fe
              do i=1, 3
                if (.not. ieee_is_finite(self%work(f,x,y,z,i))) then
                  write(stdout,'(A,I8,A,I8,A,3(I4),A,I3,A,2(I4),A,ES12.4)') &
                    'RACE: bad upstream boundary in', self%id, ' from', self%upstream%id, &
                    ' at f,x,y=', f, x, y, ' w=', i, ' i_mu,i_set=', i_mu, i_set, &
                    ' value=', self%work(f,x,y,z,i)
                  flush(stdout)
                  call shared%abort ('rt_patch_t%set_upstream_boundary: non-finite boundary')
                endif
              end do
            end do
          end do
        end do
      endif
    endif
    call trace%end(itimer)
  END SUBROUTINE set_upstream_boundary

  !===============================================================================
  SUBROUTINE integrate_intensity (self, wmu, i_set, save_PsiJ)
    class(rt_patch_t) :: self
    real(4) :: wmu, wmu2
    logical :: save_PsiJ
    integer :: i_set
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%integrate_intensity', itimer=itimer)
    wmu2 = 2*wmu
    do z=1, self%nz
      do y=1, self%ny
        do x=1, self%nx
          do f=self%fs, self%fe
            self%Jnu(f,x,y,z,i_set) = self%Jnu(f,x,y,z,i_set) + wmu2 * self%work(f,x,y,z,3)
          end do
        end do
      end do
    end do
    if (save_PsiJ) then
      do z=1, self%nz
        do y=1, self%ny
          do x=1, self%nx
            do f=self%fs, self%fe
              self%PsiJ(f,x,y,z) = self%PsiJ(f,x,y,z) + wmu2 * self%Psi(f,x,y,z)
            end do
          end do
        end do
      end do
    endif
    call trace%end(itimer)
  END SUBROUTINE integrate_intensity

  !===============================================================================
  SUBROUTINE integrate_rates (self, atom, trans, adata)
    class(rt_patch_t) :: self
    class(atom_type) :: atom
    class(atom_trans_t) :: trans
    class(atom_data_t) :: adata
    integer :: m(4), i, j, k
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch_t%integrate_rates', itimer=itimer)
    call lock_mem(self%m3d_parent, 'integrate_rates')
    if (.not. allocated(adata%gamma_ij)) then
      m = [self%patch%n, atom%nrad]
      call os%allocate(adata%gamma_ij, m, label='atom%gamma_ij', zero=.true.)
      call os%allocate(adata%gamma_ji, m, label='atom%gamma_ji', zero=.true.)
    endif
    do i=1, size(trans%line)
      k = trans%line(i)%ntrans
      adata%gamma_ij(:,:,:,k) = adata%gamma_ij(:,:,:,k) + trans%gamma_ij(:,:,:,i)
      adata%gamma_ji(:,:,:,k) = adata%gamma_ji(:,:,:,k) + trans%gamma_ji(:,:,:,i)
    end do
    do i=1, size(trans%cont)
      j = i + size(trans%line)
      k = trans%cont(i)%ntrans
      adata%gamma_ij(:,:,:,k) = adata%gamma_ij(:,:,:,k) + trans%gamma_ij(:,:,:,j)
      adata%gamma_ji(:,:,:,k) = adata%gamma_ji(:,:,:,k) + trans%gamma_ji(:,:,:,j)
    end do
    call unlock_mem(self%m3d_parent, 'integrate_rates')
    call trace%end(itimer)
  END SUBROUTINE integrate_rates

  !===============================================================================
  SUBROUTINE set_downstream_boundary (self, i_mu, i_set)
    class(rt_patch_t) :: self
    integer :: i_mu, i_set
    real(4) :: nan_check(3)
    integer :: p(3)
    integer, save :: itimer=0
    integer :: x, y, z
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%set_downstream_boundary', itimer=itimer)

    if (self%print_update_info) then
      p = self%pcheck
      if (check_NaN) then
        nan_check = self%work(self%fcheck, p(1), p(2), p(3),1:3)
        if (any(nan_check/=nan_check)) then
          print *, 'NaN in rt_patch', nan_check
        endif
      endif
      if (i_mu==self%mus .and. i_set==1) then
        memcheck(1:3) = self%work(self%fcheck, p(1), p(2), p(3),:)
        !memcheck(2:3) = memcheck(2:3) * self%twohnu3_c2(self%fcheck)
        associate(beta=>self%back(self%fcheck, p(1), p(2), p(3), 3))
          memcheck(2:3) = memcheck(2:3) * (1 - beta) / beta
        end associate
      endif
    endif

    if (.not. associated(self%mem)) then
      call self%allocate_mem()
    endif

    z = merge(self%patch%n(3), 1, self%dir > 0)
    do y=1, self%ny
      do x=1, self%nx
        self%mem(:,x,y,i_mu,i_set,1) = self%work(:,x,y,z,1)
        self%mem(:,x,y,i_mu,i_set,2) = self%work(:,x,y,z,2)
        self%mem(:,x,y,i_mu,i_set,3) = self%work(:,x,y,z,3)
      end do
    end do
    call trace%end(itimer)
  END SUBROUTINE set_downstream_boundary

  !===============================================================================
  SUBROUTINE convert_to_work_units(self)
    class(rt_patch_t) :: self
    real(4), dimension(:,:,:,:), pointer :: chi, eta, beta
    integer :: fs, fe
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    fs = self%fs
    fe = self%fe
    do z=1, self%nz
      do y=1, self%ny
        do x=1, self%nx
          do f=fs, fe
            associate(chi=>self%back(f,x,y,z,0), eta=>self%back(f,x,y,z,2), beta=>self%back(f,x,y,z,3), &
                      xhat=>self%back(f,x,y,z,1), ehat=>self%back(f,x,y,z,2))
              xhat = chi + eta
              ehat = eta / beta            ! beta is clamped nonzero (see lte_background_opacities)
            end associate
          end do
        end do
      end do
    end do
  END SUBROUTINE convert_to_work_units

  !===============================================================================
  SUBROUTINE back_to_work (self)
    class(rt_patch_t) :: self
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%back_to_work', itimer=itimer)
    ! copy self%back opacities into angle dependent opacities
    do z=1, self%nz
      do y=1, self%ny
        do x=1, self%nx
          do f=self%fs, self%fe
            associate(xhat=>self%work(f,x,y,z,1), ehat=>self%work(f,x,y,z,2))
              xhat = self%back(f,x,y,z,1)
              ehat = self%back(f,x,y,z,2)
            end associate
          end do
        end do
      end do
    end do
    call trace%end(itimer)
  END SUBROUTINE back_to_work

  !===============================================================================
  SUBROUTINE eta_to_Snu (self)
    class(rt_patch_t) :: self
    integer, save :: itimer=0
    integer :: x, y, z, f, i
    real(4) :: xhat, ehat, eta, chi, floor
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%eta_to_Snu', itimer=itimer)
    ! Converting eta into source function Snu=eta/chi
    do z=1, self%nz
      do y=1, self%ny
        do x=1, self%nx
          do f=self%fs, self%fe
            xhat = self%work(f,x,y,z,1)
            ehat = self%work(f,x,y,z,2)
            eta  = ehat * self%back(f,x,y,z,3)          ! ehat * beta
            chi  = xhat - eta
            ! Each line adds phi*Rij*(bi - bj*beta) to chi, which goes negative for
            ! an inverted population.  chi may then pass through zero, and eta/chi
            ! diverges.  Keep the sign (the inversion is real) but hold |chi| away
            ! from zero, so that Snu stays bounded by 1/maser_chi_floor.
            floor = max(real(maser_chi_floor) * abs(xhat), tiny(1.0))
            if (suppress_masers) then
              chi = max(chi, floor)
            elseif (abs(chi) < floor) then
              chi = sign(floor, chi)
              if (check_NaN) then
                !$omp atomic capture
                n_negative_chi = n_negative_chi + 1
                i = n_negative_chi
                !$omp end atomic
                if (i <= 20) then
                  write(stdout,'(A,A,4(I5),4(1x,ES13.5))') self%m3d_patch_info(), &
                    '  CHI floored f,x,y,z=', f, x, y, z, chi, xhat, eta, &
                    self%back(f,x,y,z,3)
                  flush(stdout)
                endif
              endif
            endif
            self%work(f,x,y,z,1) = chi
            self%work(f,x,y,z,2) = eta / chi
          end do
        end do
      end do
    end do
    call trace%end(itimer)
  END SUBROUTINE eta_to_Snu

  !===============================================================================
  SUBROUTINE Snu_to_eta (self)
    class(rt_patch_t) :: self
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%Snu_to_eta', itimer=itimer)
    do z=1, self%nz
      do y=1, self%ny
        do x=1, self%nx
          do f=self%fs, self%fe
            associate(chi=>self%work(f,x,y,z,1), eta=>self%work(f,x,y,z,2), beta=>self%back(f,x,y,z,3), &
                      xhat=>self%work(f,x,y,z,1), ehat=>self%work(f,x,y,z,2))
              eta = self%work(f,x,y,z,2) * chi
              ehat = eta / beta            ! beta is clamped nonzero (see lte_background_opacities)
              xhat = chi + ehat
            end associate
          end do
        end do
      end do
    end do
    call trace%end(itimer)
  END SUBROUTINE Snu_to_eta

  !===============================================================================
  SUBROUTINE rotate_array (self, vec, array)
    class(rt_patch_t) :: self
    real(8), intent(in) :: vec(3)
    real(4), intent(inout), dimension(:,:,:,:,:) :: array
    real(4), dimension(2) :: delta, deltak
    real(4), dimension(:,:,:), allocatable :: shifted
    real(4) :: bottom, bot2pivot, height, f1, f2, a, b, c, d
    integer :: up1, up2, x1, x2, y1, y2, i
    integer :: x, y, z
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin('net_opac_t%rotate_array', itimer=itimer)
    allocate(shifted, mold=array(:,:,:,1,1))
    delta = vec(:2) / (vec(3) * self%ds(:2)) ! shift in units of grid cells

    ! bot2pivot is the distance between the bottom of the patch and the pivot plane
    bottom = self%position(3) - 0.5*self%size(3)
    bot2pivot = pivot_z - bottom

    do z=1, self%nz
      height = z*self%ds(3) ! height of plane as measured from bottom
      deltak = (bot2pivot - height)*delta + ray_offset
      up1 = floor(deltak(1)) !integer part of shift
      up2 = floor(deltak(2)) !integer part of shift
      f1 = deltak(1) - up1 !fractional part
      f2 = deltak(2) - up2 !fractional part

      ! bilinear interpolation to avoid overshooting
      a = (1.0-f1)*(1.0-f2)
      b =      f1 *(1.0-f2)
      c = (1.0-f1)*     f2
      d =      f1 *     f2

      do i=1, size(array, 5)
        associate(plane => array(:,:,:,z,i))
          do y=1, self%ny
            y1 = modulo(y+up2-1, self%ny)+1
            y2 = modulo(y+up2,   self%ny)+1

            do x=1, self%nx
              x1 = modulo(x+up1-1, self%nx)+1
              x2 = modulo(x+up1,   self%nx)+1

              shifted(:,x,y) = a*plane(:,x1,y1)
              shifted(:,x,y) = b*plane(:,x2,y1) + shifted(:,x,y)
              shifted(:,x,y) = c*plane(:,x1,y2) + shifted(:,x,y)
              shifted(:,x,y) = d*plane(:,x2,y2) + shifted(:,x,y)
            end do
          end do
          plane = shifted
        end associate
      end do
    end do
    call trace%end(itimer)
  END SUBROUTINE rotate_array

  !===============================================================================
  SUBROUTINE set_fixed_Jnu (self)
    class(rt_patch_t) :: self
    real(4) :: Bnu
    integer :: nl, i, j
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%set_fixed_Jnu', itimer=itimer)
    do z=1, self%nz
      do y=1, self%ny
        do x=1, self%nx
          do f=self%fs, self%fe
            associate(beta=>self%back(f,x,y,z,3), Inu=>self%work(f,x,y,z,3))
              Bnu = beta / (1 - beta)
              Inu = fixed_Jnu * Bnu
            end associate
          end do
        end do
      end do
    end do
    do i=1, size(self%atoms, 1)
      associate(adata => self%atoms(i)%data)
        nl = size(adata%wphi, 4)
        if (.not. allocated(adata%wphi)) nl = 0
        do j=1, nl
          do z=1, self%nz
            do y=1, self%ny
              do x=1, self%nx
                adata%wphi(x,y,z,j) = 1
              end do
            end do
          end do
        end do
      end associate
    end do
    call trace%end(itimer)
  END SUBROUTINE set_fixed_Jnu

  !===============================================================================
  SUBROUTINE get_vproj (self, vproj, vec, rotate, i_mu)
    class(rt_patch_t) :: self
    real(4), dimension(:,:,:) :: vproj
    real(8), dimension(:) :: vec
    class(composition_type), pointer :: comp
    integer :: i_mu
    logical :: rotate
    integer, save :: itimer
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch_t%get_vproj', itimer=itimer)
    ! the velocities are already in units of c and already masked by use_vel, done
    ! once when the composition was loaded
    if (rotate) then
      comp => self%m3d_parent%rotated_composition(i_mu)
    else
      comp => self%m3d_parent%composition
    endif

    vproj = vec(1)*comp%vx + vec(2)*comp%vy + vec(3)*comp%vz
    call trace%end(itimer)
  END SUBROUTINE get_vproj

  !===============================================================================
  SUBROUTINE save_contribution_function (self, i_set, i_mu)
    class(rt_patch_t) :: self
    integer :: i, j, i_set, i_mu, ff, im
    real(4) :: net_chi, net_src, lin_chi, lin_eta
    real(4) :: cnt_inu, cnt_chi, cnt_eta
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch_t%save_contribution_function', itimer=itimer)
    do z=1,self%nz
      do y=1,self%ny
        do x=1,self%nx
          do f=self%fs,self%fe
            associate(cf => self%cntrbf_tmp(f,x,y,z,:), cf_out => self%cntrbf_data(f,x,y,z,:,i_mu))
              if (i_set == 1) then ! NLTE set
                net_chi = self%work(f,x,y,z,1)
                net_src = self%work(f,x,y,z,2)
                cf(1) = net_chi
                cf(2) = merge(net_src*net_chi, net_src, depression_cntrbf)
              elseif (i_set == 2) then ! continuum set
                cnt_chi = self%work(f,x,y,z,1)
                cnt_eta = self%work(f,x,y,z,2)*cnt_chi
                cnt_inu = self%work(f,x,y,z,3)

                lin_chi = cf(1) - cnt_chi
                lin_eta = cf(2) - cnt_eta

                cf_out(1) = merge(cnt_inu*lin_chi - lin_eta, cf(1)*cf(2), depression_cntrbf)
                cf_out(1) = cf_out(1) * self%twohnu3_c2(f)
                cf_out(2) = cf(1)
              endif
            end associate
          end do
        end do
      end do
    end do
    if (i_set==2 .and. i_mu==self%mue) then
      call self%save6(self%cntrbf_data, 'cntrbf')
    endif
    call trace%end(itimer)
  END SUBROUTINE save_contribution_function

  !===============================================================================
  SUBROUTINE print_post_update_info (self)
    class(rt_patch_t) :: self
    !-----------------------------------------------------------------------------
    if (self%verbose > 1) then
      call self%print_to_stdout(arr=memcheck(1:3))
    endif
  END SUBROUTINE print_post_update_info

  !===============================================================================
  SUBROUTINE rt_patch_nbor_relations(self, another, needed, needs_me, download)
    class(rt_patch_t) :: self
    class(task_t), pointer :: another        ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed,
    ! does it need me, should I download?
    class(task_t), pointer :: other
    logical :: default
    ! ----------------------------------------------------------------------------
    call trace%begin ('rt_patch_t%nbor_relations')
    call self%ray_patch_t%nbor_relations(another, needed, needs_me, download)

    if (needed .or. needs_me) then
      select type(another)
       class is (rt_patch_t)
        ! do nothing
        if (self%ipos(3) /= another%ipos(3)) then
          if (needed) then
            self%upstream => another
          endif
          if (needs_me) then
            self%dnstream => another
          endif
        else
          if (self%do_final) then
            needs_me = .false.
            if (another%i_mu==1) then
              self%up_twin => another
            elseif (another%i_mu==2) then
              self%dn_twin => another
            endif
          elseif(another%do_final) then
            needed = .false.
          else
            ! we only use dn_twin here as it doesn't really matter
            ! which one goes up and which one goes down. We assign
            ! the other direction with self in lte_background_emission
            self%dn_twin => another
          endif
        endif
       class is (m3d_patch_t)
        if (another%id == self%m3d_parent%id) then
          needs_me = .false.
        else
          select case(another%kind)
           case ('nlte_patch')
            needs_me = .true.; needed = .true.
          end select
        endif
       class default
        needed = .false.; needs_me = .false.
      end select
    endif

    call trace%end()
  END SUBROUTINE rt_patch_nbor_relations

  !===============================================================================
  !> Fast method to generate nbor lists in cases with very long task lists
  !===============================================================================
  SUBROUTINE init_nbors_hash (self, nbors)
    class(rt_patch_t):: self
    class(nbor_t), pointer:: nbors
    !.............................................................................
    class(task_t), pointer:: task
    class(nbor_t), pointer:: nbor
    integer, save:: itimer=0
    integer :: key(4)
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch_t%init_nbors_hash', itimer=itimer)
    call self%ray_patch_t%init_nbors_hash(nbors)

    key = self%key
    if (rt_out%maxiter > 0) then
      if (self%do_out_sp) then
        if (self%do_final) then
          key(3) = 1
          call self%check_nbor(key, nbors)
          key(3) = 2
          call self%check_nbor(key, nbors)
        else
          key(3) = -1
          call self%check_nbor(key, nbors)
        endif
      endif
    endif
    if (.not. self%do_final) then
      key(3) = merge(1, 2, self%i_mu==2)
      call self%check_nbor(key, nbors)
    endif

    key(2) = 0
    key(3) = 0
    key(4) = hash_string('nlte_patch')
    call self%check_nbor(key, nbors)

    key(4) = hash_string('composition')
    call self%check_nbor(key, nbors)
    call trace%end (itimer)
  END SUBROUTINE init_nbors_hash

  !===============================================================================
  SUBROUTINE finish(self)
    class(rt_patch_t) :: self
    character(8) :: Jnu_label
    logical :: cleanup
    ! ----------------------------------------------------------------------------
    call self%m3d_patch_t%finish()
    !-----------------------------------------------------------------------------
    ! Free the upstream boundary memory, but only if the upstream task will not
    ! run again.  For i_mu=2 the upstream is the patch above, which converges
    ! later than this one and is then still iterating -- taking its mem away
    ! races with the set_downstream_boundary() writing it.  Such an upstream
    ! frees its own mem on its next execute anyway, where end_boundary has now
    ! become true, because this task, its dnstream, has finished.
    !-----------------------------------------------------------------------------
    if (associated(self%upstream)) then
      if (self%upstream%is_finished()) then
        call lock_mem (self%upstream, 'finish')
        if (associated(self%upstream%mem)) then
          call self%upstream%deallocate_mem()
        endif
        call unlock_mem (self%upstream, 'finish')
      endif
    endif
    !-----------------------------------------------------------------------------
    ! Both twins hold a Jnu, and both are freed here, by whichever of the two
    ! finishes last.  They can reach this point at the same time, though, so the
    ! ownership has to be settled between them before either one frees anything.
    !-----------------------------------------------------------------------------
    if ((self%do_full_sp .and. rt_iter%maxiter > 0) .or. (self%do_final .and. rt_out%maxiter > 0)) then
      cleanup = .false.
      !$omp critical (twin_Jnu_cleanup)
      if (self%up_twin%is_finished() .and. self%dn_twin%is_finished()) then
        cleanup = .not. self%up_twin%Jnu_freed .and. .not. self%dn_twin%Jnu_freed
        if (cleanup) then
          self%up_twin%Jnu_freed = .true.
          self%dn_twin%Jnu_freed = .true.
        endif
      endif
      !$omp end critical (twin_Jnu_cleanup)
      if (cleanup) then
        if (Jnu_to_file) then
          if (.not. save_Jnu) then
            Jnu_label = merge('Jnu_out', 'Jnu    ', self%do_out_sp)
            call self%up_twin%remove(Jnu_label, key=self%up_twin%key(1:3))
            call self%dn_twin%remove(Jnu_label, key=self%dn_twin%key(1:3))
          endif
        else
          if (allocated(self%up_twin%Jnu)) then
            call mem_unregister(self%up_twin%Jnu, label='Jnu')
            deallocate(self%up_twin%Jnu)
          endif
          if (allocated(self%dn_twin%Jnu)) then
            call mem_unregister(self%dn_twin%Jnu, label='Jnu')
            deallocate(self%dn_twin%Jnu)
          endif
        endif
      endif
    endif
  END SUBROUTINE finish

  !===============================================================================
  SUBROUTINE allocate_mem(self, m_in)
    class(rt_patch_t) :: self
    class(line_type), pointer :: line
    real(4), dimension(:,:,:,:,:,:), pointer, contiguous :: tmp
    integer, intent(in), optional :: m_in(6)
    integer :: i, fs, fe, m(6)
    logical :: ok
    ! ----------------------------------------------------------------------------
    ok = self%is_set(bits%no_download) .or. self%iter > 0
    if (ok) then
      m = [self%nf, self%n(1), self%n(2), self%n_set, self%n_mu, 3]
      call self%patch_t%allocate_mem(m)
      fs = self%fs
      fe = self%fe
      self%nx = self%n(1)
      self%ny = self%n(2)
      tmp => self%mem
      self%mem(fs:fe,1:self%nx,1:self%ny,self%mus:self%mue,1:self%n_set,1:3) => tmp
    endif
  END SUBROUTINE allocate_mem

  !===============================================================================
  SUBROUTINE save_intensity (self)
    class(rt_patch_t) :: self
    real(4), dimension(:,:,:,:,:), allocatable :: ie, i3, flux
    real(8), parameter :: conversion = 2 * HH / (CC**3 * CM_TO_AA)
    integer :: fs, fe
    integer, save :: itimer=0
    integer :: x, y, f, mu, i_set
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch_t%save_intensity', itimer=itimer)
    fs = self%fs
    fe = self%fe
    allocate(i3(fs:fe,self%nx,self%ny, self%mus:self%mue, self%n_set), source=0.)
    allocate(ie(fs:fe, 1, 1, self%mus:self%mue, self%n_set), source=0.)

    do f=fs, fe
      ! this yields erg/s/cm^2/AA/sr
      i3(f,:,:,:,:) = conversion * self%nu(f)**5 * self%mem(f,:,:,:,:,3)
    end do

    do y=1, self%ny
      do x=1, self%nx
        ie(:,1,1,:,:) = ie(:,1,1,:,:) + i3(:,x,y,:,:)
      end do
    end do
    ie = ie / real(self%nx*self%ny)

    call self%save5(ie, 'ie')

    if (save_resolved) then
      call self%save5(i3, 'i3')
    endif

    ! disk-integrated flux: F(lam) = 2*pi * sum_mu rt_out%wmu(mu) * ie(lam,mu)
    ! (rt_out%wmu already carries the cosine/mu factor, see ray_mod.f90's
    ! nphi_copies), matching tumult/flux.py's integrate_flux exactly so the
    ! Python post-processing step is no longer needed for this quantity.
    allocate(flux(fs:fe, 1, 1, 1, self%n_set), source=0.)
    do i_set=1, self%n_set
      do mu=self%mus, self%mue
        flux(:,1,1,1,i_set) = flux(:,1,1,1,i_set) + real(2*pi*rt_out%wmu(mu), 4) * ie(:,1,1,mu,i_set)
      end do
    end do
    call self%save5(flux, 'flux')

    call trace%end(itimer)
  END SUBROUTINE save_intensity

  !===============================================================================
  !> Track the relative change of the top-of-domain emergent FLUX (mu-integrated,
  !> using the same weighting convention as the astrophysical flux, see
  !> integrate_intensity/save_intensity) between iterations for this wavelength
  !> chunk, and fold it into a global (all-chunk) worst-case value via the same
  !> RMA counter idiom used for the population-based convergence check in
  !> nlte_patch_t%execute. Runs every ALI iteration (not just the final pass) so
  !> both convergence metrics stay observable regardless of which one
  !> (flux_convergence) actually governs the stop decision. The convergence
  !> report itself (worst error, wavelength, i_nu) is printed by nlte_patch_t, by
  !> reading the shared RMA counters set here.
  !>
  !> rt_iter%wmu carries no cosine (mu) weighting -- it is the isotropic weight used
  !> for the mean-intensity (Jnu) angle sum, normalized so that
  !> sum(rt_iter%wmu(1:size)) = 1 over the full sphere (see integrate_intensity,
  !> which undoes the top/bottom-hemisphere-mirror halving via wmu2=2*wmu). Flux
  !> requires the additional mu = rt_iter%vec(3,i) factor, giving:
  !>   F(lam) = 2*pi * sum_i [2*rt_iter%wmu(i)] * rt_iter%vec(3,i) * I(lam,i)
  !>          = 4*pi * sum_i rt_iter%wmu(i) * rt_iter%vec(3,i) * I(lam,i)
  !> summed over the top-hemisphere ray range self%mus:self%mue. This mirrors the
  !> rt_out%wmu flux weights used by save_intensity's do_final pass (rt_out%wmu already
  !> bakes the mu factor in, see ray_mod.f90's nphi_copies), but is only inferred
  !> by analogy for the (coarser) short/ALI quadrature -- there is no other code
  !> path in the repo computing flux from rt_iter%wmu to check it against.
  !===============================================================================
  SUBROUTINE check_flux_convergence (self)
    class(rt_patch_t) :: self
    real(8), dimension(:), allocatable :: flux_now, e
    real(8) :: local_err, floor, lambda_max
    integer :: fs, fe, count, expected, round_pos, i, imax
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch_t%check_flux_convergence', itimer=itimer)
    fs = self%fs
    fe = self%fe

    allocate(flux_now(fs:fe), source=0.d0)
    do i=self%mus, self%mue
      flux_now = flux_now + (4*pi) * rt_iter%wmu(i) * rt_iter%vec(3,i) * &
        sum(sum(dble(self%mem(fs:fe,:,:,i,1,3)), dim=3), dim=2) / dble(self%nx * self%ny)
    end do

    if (.not. allocated(self%flux_prev)) then
      allocate(self%flux_prev(fs:fe), source=flux_now)
    else
      allocate(e(fs:fe))
      floor = 1d-10 * maxval(abs(flux_now))
      where (abs(self%flux_prev) > floor)
        e = flux_now / self%flux_prev - 1.d0
      elsewhere
        e = 0.d0
      endwhere
      imax = fs - 1 + maxloc(abs(e), dim=1)
      local_err = abs(e(imax))
      lambda_max = HZ_TO_AA / self%nu(imax)
      self%flux_prev = flux_now
      self%flux_maxerr = local_err

      ! expected = number of top-of-domain wavelength-chunk reporters per round
      ! (bb + bf chunks when decouple_continuum). round_pos is this report's
      ! 1-based position within its round, from the monotonic report counter:
      ! round_pos==1 starts a fresh round (reset the running max to just this
      ! report), and the decision is only ever evaluated at round_pos==expected,
      ! i.e. once EVERY reporter has folded its value into the round's running
      ! max. This replaces an earlier "same reporter seen again" reset idiom
      ! that could fire (and evaluate convergence) after only one of several
      ! reporters had contributed to the current round -- harmless with many
      ! reporters (the population check uses this same idiom, keyed by
      ! z-column, where it is very unlikely to misfire) but a near-certain
      ! false convergence with as few as 2 reporters (e.g. n_nu=1 with
      ! decouple_continuum=T).
      expected = merge(2*n_nu, n_nu, decouple_continuum)

      !$omp critical (flux_maxerr_reset)
      count = flux_update_counter%update(+1)
      round_pos = mod(count-1, expected) + 1
      if (round_pos == 1) then
        call global_flux_pos%reset(self%i_nu)
        call global_flux_err%reset(local_err)
        call global_flux_lambda%reset(lambda_max)
        self%flux_global = local_err
      else
        self%flux_global = global_flux_err%update(0.d0)
        if (abs(local_err) > abs(self%flux_global)) then
          call global_flux_pos%reset(self%i_nu)
          call global_flux_err%reset(local_err)
          call global_flux_lambda%reset(lambda_max)
          self%flux_global = local_err
        endif
      endif
      !$omp end critical (flux_maxerr_reset)

      if (round_pos == expected) then
        if ((abs(self%flux_global) > 0) .and. abs(self%flux_global) < convlim &
          .and. abs(local_err) < convlim) then
          if (flux_convergence) then
            ! Do NOT stop here: the flux was computed from populations that
            ! nlte_patch solved against the PREVIOUS radiation field earlier this
            ! iteration. Raise the cross-rank signal so nlte_patch performs one
            ! more population solve (against the now-converged field) before the
            ! run winds down. Guarded so it is raised at most once.
            if (self%iter < self%maxiter .and. flux_conv_signal%update(+0) == 0) then
              i = flux_conv_signal%update(+1)
            endif
          endif
        endif
      endif
    endif
    call trace%end(itimer)
  END SUBROUTINE check_flux_convergence

  !===============================================================================
  SUBROUTINE copy_to (self, tgt)
    class(rt_patch_t):: self
    class(patch_t):: tgt
    integer :: i, i_set, i_mu
    integer :: x, y
    !-----------------------------------------------------------------------------
    ! Invariant: the source of a download always holds a boundary buffer.  It is
    ! handed over (not copied) below, and the source re-creates one in
    ! set_downstream_boundary().  A source without mem here means someone else
    ! took it -- the download only holds a *read* lock on the source, while
    ! copy_to modifies it, so a concurrent download by the source's own dnstream,
    ! or by the source itself, can steal the buffer.  Report loudly rather than
    ! hand out a null pointer, which set_upstream_boundary() would then read.
    !-----------------------------------------------------------------------------
    if (.not. associated(self%mem)) then
      write(stdout,'(A,2(A,I8),2(A,L2))') 'RACE: copy_to with unassociated source mem', &
        ' src=', self%id, ' tgt=', tgt%id, ' src_static=', self%static_mem, &
        ' tgt_static=', tgt%static_mem
      flush(stdout)
      call shared%abort ('rt_patch_t%copy_to: source mem gone')
    endif
    if (.not. self%static_mem) then
      tgt%mem => self%mem
      nullify(self%mem)
    else
      if (.not. associated(tgt%mem)) then
        call tgt%allocate_mem()
      endif
      do i=1, 3
        do i_set=1, self%n_set
          do i_mu=self%mus, self%mue
            do y=1, self%n(2)
              do x=1, self%n(1)
                tgt%mem(:,x,y,i_mu,i_set,i) = self%mem(:,x,y,i_mu,i_set,i)
              end do
            end do
          end do
        end do
      end do
      if (tgt%static_mem) then
        call self%deallocate_mem()
      endif
    endif
  END SUBROUTINE copy_to

  !===============================================================================
  SUBROUTINE mem_register(arr, label)
    real(4), dimension(..) :: arr
    character(len=*), optional :: label
    !-----------------------------------------------------------------------------
    call os%mem_register(storage_size(arr), product(shape(arr)), label)
  END SUBROUTINE mem_register

  !===============================================================================
  SUBROUTINE mem_unregister(arr, label)
    real(4), dimension(..) :: arr
    character(len=*), optional :: label
    !-----------------------------------------------------------------------------
    call os%mem_register(storage_size(arr), -product(shape(arr)), label)
  END SUBROUTINE mem_unregister

END MODULE rt_patch_mod
