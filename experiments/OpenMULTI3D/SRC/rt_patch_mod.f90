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
  USE const_mod, only: CC, h_k, twoh_c2, half, one_third, one_sixth
  USE patch_mod, only: patch_t
  USE mmath_mod, only: cc_kernel
  USE op_math_mod, only: pchip_mono
  USE handyfuncs_mod
  USE atom_mod, only: atom_type, atom_data_t, atom_trans_t
  USE ray_patch_mod, only: ray_patch_t
  USE m3d_halo_mod, only: lock_mem, unlock_mem
  USE m3d_patch_mod, only: short_vec=>vec, short_wmu=>wmu, use_vel, m3d_patch_t, &
    absmet_in_output, long_wmu, long_vec, ray_offset, pivot_z, dlam, &
    premature_finish, m3d_patch_pre_update=>pre_update, save_resolved, &
    column_by_column, rotate_atmos, rank_lvlerr, save_cntrbf, depression_cntrbf, &
    decouple_continuum, linelist_in_output, ilambd, initial_Jnu,  &
    Jnu_to_file, preconditioning, make_eos, amr, save_Inu, save_chi, save_Snu, &
    long_solver, fixed_Jnu, n_nu, save_Jnu, maxiter
  USE ray_mod, only: one_dimensional, ray_type, rstar, rorigin
  USE trans_mod, only: line_type, cont_type
  USE bbtrans_mod, only: voigt_function2, voigt_profile
  USE bftrans_mod, only: H_occupation_probability
  USE hprofl_mod, only: hprofile
  USE m3d_debug_mod
  USE shared_mod
  USE composition_mod, only: composition_type, scattering_absorption, ref_lam, &
    departure_pointer_t, read_departures, sample_back_opac, &
    sample_linelist, sample_back_opac_in_output
  implicit none
  private

  !-----------------------------------------------------------------------------
  type, public, extends(ray_patch_t) :: rt_patch_t
    class(rt_patch_t), pointer :: upstream=>null()
    class(rt_patch_t), pointer :: dnstream=>null()
    class(rt_patch_t), pointer :: up_twin=>null()
    class(rt_patch_t), pointer :: dn_twin=>null()
    class(atom_trans_t), dimension(:), pointer :: trans=>null()
    real(4), dimension(:,:,:,:,:), allocatable :: Jnu
    real(4), dimension(:,:,:,:), allocatable :: Psi, PsiJ
    real(4), dimension(:,:,:), allocatable :: opac_table
    real(8), dimension(:), allocatable :: opac_table_nu
    real(4), dimension(:), allocatable :: twohnu3_c2
    real(8), dimension(:), pointer :: nu
    integer :: mus, mue, n_set, n_mu, natom
    logical :: bound_bound, bound_free
    logical, dimension(:), allocatable :: get_rates
  contains
    procedure :: pre_init => rt_patch_pre_init
    procedure :: init_patch_trans
    procedure :: pre_update => rt_patch_pre_update
    procedure :: execute
    procedure :: background
    procedure :: foreground
    procedure :: aniso_integrations
    procedure :: iso_integrations
    procedure :: save_debug_arrays
    procedure :: radiative_transfer
    procedure :: lte_background_opacities
    procedure :: lte_background_emission
    procedure :: generate_opacity_table
    procedure :: save_opac_table
    procedure :: interpolate_opacity_table
    procedure :: add_bound_free_opacities
    procedure :: bound_bound_damping
    procedure :: add_bound_bound_opacities
    procedure :: bound_bound_rates
    procedure :: bound_free_rates
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
    procedure :: long_char_quadratic2
    procedure :: long_char_quadratic3
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
  integer :: f, x, y, z, nx, ny, nz, nf
  real(4), dimension(:,:,:,:,:), allocatable, target :: work, back
  real(4), dimension(:,:,:,:,:,:), allocatable :: cntrbf_data
  real(4), dimension(:,:,:,:,:), allocatable :: cntrbf_tmp
  !$omp threadprivate(work, back, x, y, z, nx, ny, nz, f, nf)
  !$omp threadprivate(cntrbf_data, cntrbf_tmp)

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
    self%n_mu = merge(size(long_wmu), size(short_wmu)/2, self%do_final)
    self%mus = merge(1, self%n_mu+1, self%dir>0)
    self%mue = merge(self%n_mu, size(short_wmu), self%dir>0)
    self%natom = merge(size(self%atoms), 0, associated(self%atoms))
    is_nlte = maxiter > 0 .or. ilambd > 0 .or. read_departures

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
            trans%line(j)%ired = fs
            trans%line(j)%iblue = fe
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
            trans%cont(j)%ired = fs
            trans%cont(j)%iblue = fe
            j = j + 1
          end do
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
    logical :: skip, get_cntrbf, end_boundary, do_rotate_atmos
    integer :: i, j, i_mu, zs, ze, fs, fe, ng, i_set, nrad
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%execute', itimer=itimer)
    fs = self%fs
    fe = self%fe
    nx = self%patch%n(1)
    ny = self%patch%n(2)
    nz = self%patch%n(3)
    ng = self%ng(3)
    zs = 1 - ng
    ze = nz + ng

    self%nu(fs:fe) => self%sp%nu(fs:fe)
    do_rotate_atmos = rotate_atmos .and. self%do_final
    get_cntrbf = self%do_final .and. save_cntrbf
    self%get_rates = self%get_rates .and. [(.not. self%atoms(i)%data%converged, i=1, self%natom)]
    end_boundary = self%patch%boundaries%is_set(bits%zl) .and. (self%dir<0)
    end_boundary = self%patch%boundaries%is_set(bits%zu) .and. (self%dir>0) .or. end_boundary
    if (.not. end_boundary .and. associated(self%dnstream)) then
      end_boundary = self%dnstream%is_finished()
    endif

    call allocate_arrays()

    if (do_rotate_atmos) then ! Final output spectrum from ray-aligned atmosphere
      do i_mu=self%mus, self%mue
        wmu = long_wmu(i_mu)
        vec = long_vec(:,i_mu)
        comp => self%m3d_parent%rotated_composition(i_mu)
        call self%get_vproj(vproj, vec, do_rotate_atmos, i_mu)
        do i_set=1, self%n_set
          call self%background(comp, i_mu, i_set, do_rotate_atmos)
          call self%back_to_work()
          if (i_set /= 2) call self%foreground(comp, i_mu, i_set, do_rotate_atmos, wmu, vproj)
          call self%eta_to_Snu()
          call self%set_upstream_boundary(i_mu, i_set, ng)
          call self%radiative_transfer(comp, vec, do_rotate_atmos)
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
          wmu = long_wmu(i_mu)
          vec = long_vec(:,i_mu)
          call self%get_vproj(vproj, vec, do_rotate_atmos, i_mu)
          call self%back_to_work()
          if (i_set /= 2) call self%foreground(comp, i_mu, i_set, do_rotate_atmos, wmu, vproj)
          call self%eta_to_Snu()
          call self%set_upstream_boundary(i_mu, i_set, ng)
          call self%radiative_transfer(comp, vec, do_rotate_atmos)
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
          wmu = short_wmu(i_mu)
          vec = short_vec(:,i_mu)
          call self%get_vproj(vproj, vec, do_rotate_atmos, i_mu)
          if (self%do_full_sp .and. fixed_Jnu > 0) then
            call self%set_fixed_Jnu()
          else
            call self%back_to_work()
            call self%foreground(comp, i_mu, i_set, do_rotate_atmos, wmu, vproj)
            call self%eta_to_Snu()
            call self%set_upstream_boundary(i_mu, i_set, ng)
            call self%radiative_transfer(comp, vec, do_rotate_atmos)
          endif
          call self%set_downstream_boundary(i_mu, i_set)
          call self%aniso_integrations(comp, i_mu, i_set, do_rotate_atmos, wmu, vproj)
        end do
        call self%iso_integrations()
      end do
    endif

    if (end_boundary) then
      if (self%do_final .and. self%n_mu>0) then
        call self%save_intensity()
        if (self%i_nu == -1) then
          call self%print_linelist(comp)
        endif
      endif
      call self%deallocate_mem()
    endif

    call deallocate_arrays()

    if (self%has_atom .and. self%iter > ilambd-1) then
      if (all([(self%atoms(i)%data%converged, i=1, self%natom)]) .and. premature_finish) then
        call self%finish()
      endif
    endif

    call trace%end(itimer)

  contains

    !-----------------------------------------------------------------------------
    subroutine allocate_arrays()
      allocate(back(fs:fe,nx,ny,nz,0:3), source=0.)
      allocate(work(fs:fe,nx,ny,zs:ze,3), source=0.)
      allocate(self%twohnu3_c2(fs:fe))
      allocate(vproj(nx,ny,nz))
      self%twohnu3_c2 = twoh_c2 * self%nu**3
      call mem_register(back,      label='back')
      call mem_register(work,      label='work')
      if (.not. allocated(self%Jnu)) then
        allocate(self%Jnu(fs:fe,nx,ny,nz,self%n_set), source=0.)
        call mem_register(self%Jnu, label='Jnu')
      endif
      do i=1, self%natom
        trans => self%trans(i)
        if (self%bound_bound) then
          call os%allocate(trans%adamp, [nx,ny,nz,size(trans%line)], label='adamp', zero=.true.)
        endif
        if (self%get_rates(i)) then
          nrad = size(trans%line) + size(trans%cont)
          call os%allocate(trans%gamma_ij, [nx,ny,nz,nrad], label='gamma_ij', zero=.true.)
          call os%allocate(trans%gamma_ji, [nx,ny,nz,nrad], label='gamma_ji', zero=.true.)
        endif
      end do
      if (any(self%get_rates)) then
        allocate(self%Psi(fs:fe,nx,ny,nz), source=0.)
        call mem_register(self%Psi, label='Psi')
        if (self%bound_free) then
          allocate(self%PsiJ(fs:fe,nx,ny,nz), source=0.)
          call mem_register(self%PsiJ, label='PsiJ')
        endif
      endif
      if (get_cntrbf) then
        allocate(cntrbf_data(fs:fe,nx,ny,nz,2,self%mus:self%mue))
        call mem_register(cntrbf_data, 'cntrbf_data')
        allocate(cntrbf_tmp(fs:fe,nx,ny,nz,2))
        call mem_register(cntrbf_tmp, 'cntrbf_tmp')
      endif
    end subroutine allocate_arrays

    !-----------------------------------------------------------------------------
    subroutine deallocate_arrays()
      call mem_unregister(work     , label='work')
      call mem_unregister(back     , label='back')
      deallocate(work, back)
      deallocate(self%twohnu3_c2)
      do i=1, self%natom
        trans => self%trans(i)
        if (self%bound_bound) then
          call os%deallocate(trans%adamp, label='adamp')
        endif
        if (self%get_rates(i)) then
          call os%deallocate(trans%gamma_ij, label='gamma_ij')
          call os%deallocate(trans%gamma_ji, label='gamma_ji')
        endif
      end do
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
        call mem_unregister(cntrbf_data, label='cntrbf_data')
        deallocate(cntrbf_data)
        call mem_unregister(cntrbf_tmp, label='cntrbf_tmp')
        deallocate(cntrbf_tmp)
      endif
    end subroutine deallocate_arrays
  END SUBROUTINE execute

  !===============================================================================
  SUBROUTINE background(self, comp, i_mu, i_set, do_rotate_atmos)
    class(rt_patch_t) :: self
    class(composition_type), pointer :: comp
    class(atom_data_t), pointer :: adata
    integer, intent(in) :: i_mu, i_set
    logical, intent(in) :: do_rotate_atmos
    real(8) :: vec(3)
    logical :: include_absmet, do_sample_opac, do_make_eos
    integer :: i
    !-----------------------------------------------------------------------------
    include_absmet = (self%do_full_sp .and. (self%i_nu > self%n_nu .or. .not. decouple_continuum)) &
                   .or. (self%do_out_sp .and. absmet_in_output)
    do_sample_opac = (self%do_full_sp .and. sample_back_opac) &
                   .or. (self%do_out_sp .and. sample_back_opac_in_output)
    do_make_eos = make_eos .and. self%do_out_sp

    if (do_rotate_atmos .and. i_mu /= 0) then
      vec = long_vec(:,i_mu)
    endif

    if (i_set==1) then
      call self%lte_background_opacities(comp, include_absmet, do_sample_opac, do_make_eos)
    endif
    call self%lte_background_emission(comp, i_set, vec, include_absmet, do_rotate_atmos)
    call self%convert_to_work_units()
    call self%reset_array(self%Jnu(:,:,:,:,i_set:i_set))
    call self%reset_array(work)
    
    do i=1, self%natom
      associate(atom=>self%atoms(i), trans=>self%trans(i))
        adata => atom%data
        if (do_rotate_atmos .and. i_mu /= 0) then
          adata => atom%rotated_data(i_mu)
        endif

        if (self%bound_free) then
          call self%add_bound_free_opacities(atom, trans, adata)
        endif
        if (self%bound_bound .and. i_set==1) then
          call self%bound_bound_damping(comp, atom, trans, adata)
        endif
      end associate
    end do
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
    logical :: use_linelist
    integer :: i
    !-----------------------------------------------------------------------------
    use_linelist = self%do_final .and. comp%ll_initialised .and. linelist_in_output .and. .not. sample_linelist
    
    do i=1, self%natom
      associate(atom=>self%atoms(i), trans=>self%trans(i))
        adata => atom%data
        if (do_rotate_atmos .and. i_mu /= 0) then
          adata => atom%rotated_data(i_mu)
        endif
        if (self%bound_bound) then
          call self%add_bound_bound_opacities(comp, atom, trans, adata, wmu, vproj, LTE=i_set==3)
        endif
      end associate
    end do

    if (use_linelist) then
      departures => null()
      if (associated(self%m3d_parent%departures) .and. i_set==1) then
        i = merge(i_mu+1, 1, do_rotate_atmos)
        departures => self%m3d_parent%departures(:,i)
      endif
      associate(xhat=>work(:,:,:,1:nz,1), ehat=>work(:,:,:,1:nz,2))
        call comp%add_linelist_opacities(self%nu, xhat, ehat, dlam, vproj, departures=departures)
      end associate
    endif
  END SUBROUTINE foreground

  !===============================================================================
  SUBROUTINE radiative_transfer(self, comp, vec, do_rotate_atmos)
    class(rt_patch_t) :: self
    class(composition_type) :: comp
    real(8), intent(in) :: vec(3)
    logical, intent(in) :: do_rotate_atmos
    !-----------------------------------------------------------------------------
    if (self%do_final) then
      if (.not. (one_dimensional .or. column_by_column .or. do_rotate_atmos)) then
        call self%rotate_array(vec, work(:,:,:,1:nz,1:2))
      endif

      select case(long_solver)
        case('linear')
        !if (one_dimensional .and. rstar > 2) then
        !  call self%long_char_linear_spherical(vec)
        !else
        call self%long_char_linear(vec)
        !endif
        case('quadratic')
        call self%long_char_quadratic(vec)
        case('quadratic2')
        call self%long_char_quadratic2(vec)
        case('quadratic3')
        call self%long_char_quadratic3(comp, vec)
        case default
        stop 'Unknown long_solver, use "linear" or "quadratic"'
      end select
    else
      if (one_dimensional .or. column_by_column) then
        call self%long_char_linear(vec)
      else
        call self%short_char(vec)
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
    get_Jnu = ((maxiter > 0) .or. (ilambd > 0)) .and. .not. self%do_final
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
          call self%bound_bound_rates(comp, atom, trans, adata, wmu, vproj)
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

    if (Jnu_to_file) then
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
        call self%save4(work(:,:,:,1:nz,1), sav_label, key=self%key(1:3))
      endif
    endif
    if (save_Snu) then
      if (i_set==1) then
        write(sav_label, '(A, I2.2)') 'Snu_', i_mu
        call self%save4(work(:,:,:,1:nz,2), sav_label, key=self%key(1:3))
      endif
    endif
    if (save_Inu) then
      if (i_set==1) then
        write(sav_label, '(A, I2.2)') 'Inu_', i_mu
        call self%save4(work(:,:,:,1:nz,3), sav_label, key=self%key(1:3))
      elseif (i_set==2) then
        write(sav_label, '(A, I2.2)') 'Cnt_', i_mu
        call self%save4(work(:,:,:,1:nz,3), sav_label, key=self%key(1:3))
      endif
    endif
  END SUBROUTINE save_debug_arrays


  !===============================================================================
  SUBROUTINE lte_background_opacities (self, comp, include_absmet, do_sample_opac, do_make_eos)
    class(rt_patch_t) :: self
    class(composition_type) :: comp
    real(4) :: h_kT, hv_kT, bnu
    real(8), dimension(:), allocatable, target :: tab_nu
    logical :: include_absmet, do_sample_opac, do_make_eos
    integer :: fs, fe
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%lte_background_opacities', itimer=itimer)
    fs = self%fs
    fe = self%fe

    do z=1, nz
      do y=1, ny
        do x=1, nx
          ! computing boltzmann factor
          h_kT = h_k / real(comp%temp(x,y,z))
          do f=fs, fe
            hv_kT = self%nu(f) * h_kT
            back(f,x,y,z,3) = exp(-hv_kT)
          end do
        end do
      end do
    end do

    if (.not. allocated(self%opac_table)) then
      if (do_make_eos) then
        call self%generate_opacity_table(self%nu(fs:fe), include_absmet, comp%ll_initialised)
      elseif (do_sample_opac) then
        if (include_absmet .or. sample_linelist) then
          tab_nu = self%nu(fs:fe)
        else
          call comp%get_table_frequencies(self%nu(fs:fe), tab_nu)
        endif
        call self%generate_opacity_table(tab_nu, include_absmet, include_llist=sample_linelist)
      endif
    endif

    ! getting background LTE absorption
    if (do_make_eos .and. self%iter==0) then
      call self%save_opac_table()
    elseif (do_sample_opac) then
      call self%interpolate_opacity_table(comp%temp, comp%ne, comp%chi500)
    else
      call comp%absorption(self%nu, back(:,:,:,:,0), &
        beta=back(:,:,:,:,3), include_absmet=include_absmet, sorted=self%sp%sorted)
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

    if (self%iter > 0 .or. (self%do_final .and. ilambd > 0)) then
      Jnu_to_Bnu = .false.
      if (Jnu_to_file) then
        Jnu_label = merge('Jnu_out', 'Jnu    ', self%do_out_sp)
        allocate(Jup_mem(fs:fe,nx,ny,nz,self%up_twin%n_set))
        allocate(Jdn_mem(fs:fe,nx,ny,nz,self%dn_twin%n_set))
        call self%up_twin%load5(Jnu_label, key=self%up_twin%key(1:3), mem=Jup_mem)
        call self%dn_twin%load5(Jnu_label, key=self%dn_twin%key(1:3), mem=Jdn_mem)
        Jup => Jup_mem
        Jdn => Jdn_mem
      else
        Jup => self%up_twin%Jnu
        Jdn => self%dn_twin%Jnu
      endif
      do z=1, nz
        do y=1, ny
          do x=1, nx
            do f=fs, fe
              self%Jnu(f,x,y,z,i_set) = 0.5*(Jup(f,x,y,z,i_set) + Jdn(f,x,y,z,i_set))
            end do
          end do
        end do
      end do
      if (do_rotate_atmos .and. .not. (one_dimensional .or. column_by_column)) then
        call self%rotate_array(vec, self%Jnu(:,:,:,:,i_set:i_set))
      endif
    elseif (.not. Jnu_to_Bnu) then
      call self%reset_array(self%Jnu(:,:,:,:,i_set:i_set))
    endif

    do z=1, nz
      do y=1, ny
        do x=1, nx
          do f=fs, fe
            associate(Jnu=>self%Jnu(f,x,y,z,i_set), chi=>back(f,x,y,z,0), eta=>back(f,x,y,z,2), beta=>back(f,x,y,z,3))
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

    call comp%scattering(self%nu, back(:,:,:,:,2), self%Jnu(:,:,:,:,i_set), include_absmet=include_absmet)
    call trace%end(itimer)
  END SUBROUTINE

  !===============================================================================
  SUBROUTINE generate_opacity_table (self, nu, include_absmet, include_llist)
    class(rt_patch_t), target :: self
    class(composition_type), pointer :: comp
    real(8), dimension(:), pointer, intent(in) :: nu
    real(4), dimension(:,:,:,:), pointer :: all_chi
    real(4), dimension(:,:,:,:), allocatable :: all_eta, all_beta
    logical, intent(in) :: include_absmet, include_llist
    real(4) :: over_chi500
    integer :: i,j,k, ntemp, nne, nnu
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%generate_opacity_table', itimer=itimer)
    comp => self%m3d_parent%tab_comp
    nnu = size(nu)
    ntemp = size(comp%temp,1)
    nne   = size(comp%temp,2)
    allocate(self%opac_table(nnu, ntemp, nne), source=0.)
    allocate(self%opac_table_nu, source=nu)
    call mem_register(self%opac_table, label='opac_table')

    all_chi(1:nnu, 1:ntemp, 1:nne, 1:1) => self%opac_table
    call comp%absorption(nu, all_chi, include_absmet=include_absmet, sorted=self%sp%sorted)

    if (include_llist) then
      ! absorption returns the actual chi, while add_linelist_opacities
      ! adds the working versions of chi and eta
      allocate(all_eta, all_beta, mold=all_chi)
      do k=1, nne
        do j=1, ntemp
          associate(chi=>all_chi(:,j,k,1), eta=>all_eta(:,j,k,1), beta=>all_beta(:,j,k,1), &
                    xhat=>all_chi(:,j,k,1), ehat=>all_eta(:,j,k,1))
            beta = exp(-h_k*nu / comp%temp(j,k,1))
            ehat = chi / (1 - beta)
            xhat = chi + ehat*beta
          end associate
        end do
      end do
      call comp%add_linelist_opacities(nu, all_chi, all_eta, dlam)
      do k=1, nne
        do j=1, ntemp
          associate(chi=>all_chi(:,j,k,1), eta=>all_eta(:,j,k,1), beta=>all_beta(:,j,k,1), &
                    xhat=>all_chi(:,j,k,1), ehat=>all_eta(:,j,k,1))
            if (eta > 0) then
              ehat = eta / beta
            endif
            chi = xhat - eta
          end associate
        end do
      end do
    endif

    do k=1, nne
      do j=1, ntemp
        over_chi500 = 1 / comp%chi500(j,k,1)
        self%opac_table(:,j,k) = self%opac_table(:,j,k) * over_chi500
      end do
    end do
    call trace%end(itimer)
  END SUBROUTINE generate_opacity_table

  !===============================================================================
  SUBROUTINE save_opac_table (self)
    class(rt_patch_t), target :: self
    class(composition_type), pointer :: comp
    real(4), dimension(:,:,:), allocatable :: eos
    real(4), dimension(:,:,:,:), pointer :: view
    real(8), dimension(:), pointer :: nu
    real(4) :: over_rho
    integer :: i,j,k, ntemp, nrho, nnu
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%save_opac_table', itimer=itimer)
    comp => self%m3d_parent%tab_comp
    nu => self%opac_table_nu
    ntemp = size(self%opac_table, 2)
    nrho = size(self%opac_table, 3)
    nnu = size(nu)

    ! save absorption part first (which has been normalised by chi500)
    do k=1, nrho
      do j=1, ntemp
        over_rho = comp%chi500(j,k,1) / comp%rho(j,k,1)
        self%opac_table(:,j,k) = self%opac_table(:,j,k) * over_rho
      end do
    end do
    call self%save3(self%opac_table, 'opa')

    ! then save scattering part
    self%opac_table = 0.0
    view(1:nnu, 1:ntemp, 1:nrho, 1:1) => self%opac_table
    call comp%electron_scattering(view)
    call comp%rayleigh_scattering(cc / nu * 1e8, view)
    do k=1, nrho
      do j=1, ntemp
        over_rho = 1 / comp%rho(j,k,1)
        self%opac_table(:,j,k) = self%opac_table(:,j,k) * over_rho
      end do
    end do
    call self%save3(self%opac_table, 'scat')
    call trace%end(itimer)
  END SUBROUTINE save_opac_table

  !===============================================================================
  SUBROUTINE interpolate_opacity_table (self, temp, ne, chi500)
    class(rt_patch_t), target :: self
    class(composition_type), pointer :: comp
    real(8), dimension(:,:,:), intent(in) :: temp, ne
    real(4), dimension(:,:,:), intent(in) :: chi500
    integer, dimension(:), allocatable :: inu
    real(4), dimension(:), allocatable :: wnu
    real(8), dimension(:), pointer :: tab_nu
    real(8) :: w1, w2, e, t, mint, mine, kappa
    real(8) :: dx, dy, a, b, c, d, fi1, fi2, fj1, fj2, ab, cd
    integer :: fs, fe, i, j, k
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%interpolate_opacity_table', itimer=itimer)
    comp => self%m3d_parent%tab_comp
    dx = 1. / (comp%temp(2,1,1) - comp%temp(1,1,1))
    dy = 1. / (log(comp%ne(1,2,1) / comp%ne(1,1,1)))
    mint = comp%temp(1,1,1)
    mine = log(comp%ne(1,1,1))
    tab_nu => self%opac_table_nu

    fs = self%fs
    fe = self%fe

    allocate(inu(fs:fe), wnu(fs:fe))
    inu(fs:fe) = [(searchsorted(tab_nu, self%sp%nu(f)), f=fs, fe)]
    wnu(fs:fe) = [((self%sp%nu(f) - tab_nu(inu(f))) / (tab_nu(inu(f)+1) - tab_nu(inu(f))), f=fs, fe)]

    do z=1, nz
      do y=1, ny
        do x=1, nx
          t = temp(x,y,z)
          e = log(ne(x,y,z))
          fi1 = (t - mint) * dx + 1
          fj1 = (e - mine) * dy + 1
          i = floor(fi1)
          j = floor(fj1)
          fi1 = fi1 - i
          fj1 = fj1 - j
          fi2 = 1 - fi1
          fj2 = 1 - fj1
          do f=fs, fe
            w2 = wnu(f)
            w1 = 1 - w2
            k = inu(f)
            a = w1*self%opac_table(k, i,   j)   + w2*self%opac_table(k+1, i,   j)
            b = w1*self%opac_table(k, i+1, j)   + w2*self%opac_table(k+1, i+1, j)
            c = w1*self%opac_table(k, i,   j+1) + w2*self%opac_table(k+1, i,   j+1)
            d = w1*self%opac_table(k, i+1, j+1) + w2*self%opac_table(k+1, i+1, j+1)
            ab    = fi2*a  + fi1*b
            cd    = fi2*c  + fi1*d
            kappa = fj2*ab + fj1*cd
            back(f,x,y,z,0) = kappa*chi500(x,y,z)
          end do
        end do
      end do
    end do
    call trace%end(itimer)
  END SUBROUTINE interpolate_opacity_table

  !===============================================================================
  SUBROUTINE reset_array (self, array)
    class(rt_patch_t), target :: self
    real(4), dimension(:,:,:,:,:), intent(inout) :: array
    integer :: i
    !-----------------------------------------------------------------------------
    do i=lbound(array, 5), ubound(array,5)
      do z=lbound(array, 4), ubound(array,4)
        do y=1, ny
          do x=1, nx
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
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%add_bound_free_opacities', itimer=itimer)
    
    ! Summing contributions from all bound-free transitions
    do i=1, size(trans%cont)
      cont => trans%cont(i)
      if (.not. cont%active) cycle
      fs = cont%ired
      fe = cont%iblue

      if (allocated(adata%b)) then
        do z=1, nz
          do y=1, ny
            do x=1, nx
              ni  = adata%nstar(x,y,z,cont%i) * adata%b(x,y,z,cont%i)
              tmp = adata%nstar(x,y,z,cont%i) * adata%b(x,y,z,cont%j)
              do f=fs, fe
                alpha = cont%alpha(f)
                back(f,x,y,z,1) = back(f,x,y,z,1) + alpha * ni  ! chi
                back(f,x,y,z,2) = back(f,x,y,z,2) + alpha * tmp ! eta
              end do
            end do
          end do
        end do
      else ! b = 1
        do z=1, nz
          do y=1, ny
            do x=1, nx
              ni  = adata%nstar(x,y,z,cont%i)
              tmp = adata%nstar(x,y,z,cont%i)
              do f=fs, fe
                alpha = cont%alpha(f)
                back(f,x,y,z,1) = back(f,x,y,z,1) + alpha * ni  ! chi
                back(f,x,y,z,2) = back(f,x,y,z,2) + alpha * tmp ! eta
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
  SUBROUTINE add_bound_bound_opacities (self, comp, atom, trans, adata, wmu, vproj, LTE)
    class(rt_patch_t), target :: self
    class(composition_type) :: comp
    class(atom_type) :: atom
    class(atom_trans_t) :: trans
    class(atom_data_t), target :: adata
    class(line_type), pointer :: line
    real(4), contiguous, dimension(:,:,:) :: vproj
    real(4), dimension(:,:,:), allocatable :: wphi, temp, ne, fo, vbar2h, tothi, t3nhe
    real(4), dimension(:), allocatable :: vphi
    real(4) :: wmu, idf, vpr, a, idf_4pi, ni, nj, Bji_nj, Bij_ni, widf, wp, amin, amax
    real(8) :: q
    real(4) :: lower_off, upper_off, ekt, wj, chi_fac, eta_fac, bi, bj, Rij
    integer :: i, ih, fs, fe
    logical :: LTE, do_LTE, get_wphi
    integer, save :: itimer=0, itimer1=0
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
      fs = line%ired
      fe = line%iblue

      lower_off = line%lower_ev - atom%ev(line%i)
      upper_off = line%upper_ev - atom%ev(line%j)

      ! Adding bound-bound opacities
      do z=1, nz
        do y=1, ny
          do x=1, nx
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

              associate(xhat=>work(fs:fe,x,y,z,1), ehat=>work(fs:fe,x,y,z,2))
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
    call trace%end(itimer)
  END SUBROUTINE add_bound_bound_opacities

  !===============================================================================
  SUBROUTINE bound_bound_rates (self, comp, atom, trans, adata, wmu, vproj)
    class(rt_patch_t), target :: self
    class(composition_type) :: comp
    class(atom_type) :: atom
    class(atom_trans_t) :: trans
    class(atom_data_t), target :: adata
    class(line_type), pointer :: line
    real(4), contiguous, dimension(:,:,:) :: vproj
    real(4), dimension(:,:,:), allocatable :: temp, ne, fo, vbar2h, tothi, t3nhe
    real(4), dimension(:), allocatable :: wmunu
    real(4) :: wmu, idf, vpr, a, idf_4pi, ni, nj, chi_fac, phi, vel, vgt
    real(4) :: vvv, a2, sum_I, sum_P, sum_A
    real(4) :: g_ij, g_ji, q
    real(4) :: lower_off, upper_off, ekt
    integer :: i, l, ih, fs, fe
    integer, save :: itimer=0
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
      fs = line%ired
      fe = line%iblue
      wmunu(fs:fe) = wmu * line%wnu(fs:fe) * line%Aij
      lower_off = line%lower_ev - atom%ev(line%i)
      upper_off = line%upper_ev - atom%ev(line%j)

      do z=1,nz
        do y=1,ny
          do x=1,nx
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
                associate(Inu=>work(f,x,y,z,3), Psi=>self%Psi(f,x,y,z))
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
                associate(Inu=>work(f,x,y,z,3), Psi=>self%Psi(f,x,y,z))
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
                associate(Inu=>work(f,x,y,z,3), Psi=>self%Psi(f,x,y,z))
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
    !-----------------------------------------------------------------------------
    call trace%begin (trim(self%kind)//'%bound_free_rates', itimer=itimer)
    allocate(wmunu(self%fs:self%fe))
    allocate(w_alpha(self%fs:self%fe))
    ! Computing bound-free radiative rates based on half integrated mean intensity
    do j=1, size(trans%cont)
      cont => trans%cont(j)
      if (.not. cont%active) cycle
      i = j + size(trans%line)
      fs = cont%ired
      fe = cont%iblue
      wmunu(fs:fe) = 0.5 * cont%wnu(fs:fe)
      w_alpha(fs:fe) = self%twohnu3_c2(fs:fe) * cont%alpha(fs:fe) * wmunu(fs:fe)

      do z=1,nz
        do y=1,ny
          do x=1,nx
            sum_abs = 0.0
            sum_ems = 0.0
            sum_pre = 0.0

            do f=fs, fe
              associate(Jnu=>self%Jnu(f,x,y,z,1), beta=>back(f,x,y,z,3), Psi=>self%PsiJ(f,x,y,z))
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
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%set_upstream_boundary', itimer=itimer)
    z = merge(lbound(work,4), ubound(work,4), self%dir > 0)

    if (self%boundaries%is_set(bits%zu)) then
      work(:,:,:,z,:) = 0.
    elseif (self%boundaries%is_set(bits%zl)) then
      xxold = [(i*self%ds(3), i=1, 3)]
      xxnew = [(i*self%ds(3), i=1-ng, 0)]
      do i=1, 2
        do y=1,  ny
          do x=1,  nx
            do f=self%fs, self%fe
              call pchip_mono(xxold, work(f,x,y,z+1:z+3,i), &
                xxnew, work(f,x,y,z:z,i), sorted=.true.)
              work(f,x,y,z,i) = max(work(f,x,y,z,i), 0.)
            end do
          end do
        end do
      end do
      do y=1,  ny
        do x=1,  nx
          do f=self%fs, self%fe
            work(f,x,y,z,3) = work(f,x,y,z,2)
          end do
        end do
      end do
    else
      do y=1,  ny
        do x=1,  nx
          do f=self%fs, self%fe
            work(f,x,y,z,1) = self%mem(f,x,y,i_mu,i_set,1)
            work(f,x,y,z,2) = self%mem(f,x,y,i_mu,i_set,2)
            work(f,x,y,z,3) = self%mem(f,x,y,i_mu,i_set,3)
          end do
        end do
      end do
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
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%integrate_intensity', itimer=itimer)
    wmu2 = 2*wmu
    do z=1, nz
      do y=1, ny
        do x=1, nx
          do f=self%fs, self%fe
            self%Jnu(f,x,y,z,i_set) = self%Jnu(f,x,y,z,i_set) + wmu2 * work(f,x,y,z,3)
          end do
        end do
      end do
    end do
    if (save_PsiJ) then
      do z=1, nz
        do y=1, ny
          do x=1, nx
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
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%set_downstream_boundary', itimer=itimer)

    if (self%print_update_info) then
      p = self%pcheck
      if (check_NaN) then
        nan_check = work(self%fcheck, p(1), p(2), p(3),1:3)
        if (any(nan_check/=nan_check)) then
          print *, 'NaN in rt_patch', nan_check
        endif
      endif
      if (i_mu==self%mus .and. i_set==1) then
        memcheck(1:3) = work(self%fcheck, p(1), p(2), p(3),:)
        !memcheck(2:3) = memcheck(2:3) * self%twohnu3_c2(self%fcheck)
        associate(beta=>back(self%fcheck, p(1), p(2), p(3), 3))
          memcheck(2:3) = memcheck(2:3) * (1 - beta) / beta
        end associate
      endif
    endif

    if (.not. associated(self%mem)) then
      call self%allocate_mem()
    endif

    z = merge(self%patch%n(3), 1, self%dir > 0)
    do y=1, ny
      do x=1, nx
        self%mem(:,x,y,i_mu,i_set,1) = work(:,x,y,z,1)
        self%mem(:,x,y,i_mu,i_set,2) = work(:,x,y,z,2)
        self%mem(:,x,y,i_mu,i_set,3) = work(:,x,y,z,3)
      end do
    end do
    call trace%end(itimer)
  END SUBROUTINE set_downstream_boundary

  !===============================================================================
  SUBROUTINE convert_to_work_units(self)
    class(rt_patch_t) :: self
    real(4), dimension(:,:,:,:), pointer :: chi, eta, beta
    integer :: fs, fe
    !-----------------------------------------------------------------------------
    fs = self%fs
    fe = self%fe
    do z=1, nz
      do y=1, ny
        do x=1, nx
          do f=fs, fe
            associate(chi=>back(f,x,y,z,0), eta=>back(f,x,y,z,2), beta=>back(f,x,y,z,3), &
                      xhat=>back(f,x,y,z,1), ehat=>back(f,x,y,z,2))
              xhat = chi + eta
              if (eta > 0) then
                ehat = eta / beta
              endif
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
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%back_to_work', itimer=itimer)
    ! copy back opacities into angle dependent opacities
    do z=1, nz
      do y=1, ny
        do x=1, nx
          do f=self%fs, self%fe
            associate(xhat=>work(f,x,y,z,1), ehat=>work(f,x,y,z,2))
              xhat = back(f,x,y,z,1)
              ehat = back(f,x,y,z,2)
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
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%eta_to_Snu', itimer=itimer)
    ! Converting eta into source function Snu=eta/chi
    do z=1, nz
      do y=1, ny
        do x=1, nx
          do f=self%fs, self%fe
            associate(chi=>work(f,x,y,z,1), eta=>work(f,x,y,z,2), beta=>back(f,x,y,z,3), &
                      xhat=>work(f,x,y,z,1), ehat=>work(f,x,y,z,2))
              eta = ehat * beta
              chi = xhat - eta
              work(f,x,y,z,2) = eta / chi
            end associate
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
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%Snu_to_eta', itimer=itimer)
    do z=1, nz
      do y=1, ny
        do x=1, nx
          do f=self%fs, self%fe
            associate(chi=>work(f,x,y,z,1), eta=>work(f,x,y,z,2), beta=>back(f,x,y,z,3), &
                      xhat=>work(f,x,y,z,1), ehat=>work(f,x,y,z,2))
              eta = work(f,x,y,z,2) * chi
              ehat = eta / beta
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
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin('net_opac_t%rotate_array', itimer=itimer)
    allocate(shifted, mold=array(:,:,:,1,1))
    delta = vec(:2) / (vec(3) * self%ds(:2)) ! shift in units of grid cells

    ! bot2pivot is the distance between the bottom of the patch and the pivot plane
    bottom = self%position(3) - 0.5*self%size(3)
    bot2pivot = pivot_z - bottom

    do z=1, nz
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
          do y=1, ny
            y1 = modulo(y+up2-1, ny)+1
            y2 = modulo(y+up2,   ny)+1

            do x=1, nx
              x1 = modulo(x+up1-1, nx)+1
              x2 = modulo(x+up1,   nx)+1

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
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%set_fixed_Jnu', itimer=itimer)
    do z=1, nz
      do y=1, ny
        do x=1, nx
          do f=self%fs, self%fe
            associate(beta=>back(f,x,y,z,3), Inu=>work(f,x,y,z,3))
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
          do z=1, nz
            do y=1, ny
              do x=1, nx
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
  SUBROUTINE short_char (self, vec)
    class(rt_patch_t) :: self
    type(ray_type) :: ray, amr_ray
    real(8), intent(in) :: vec(3)
    real(4) :: Inuup, xxup, ssup, xxmin, ssmin, Inumin
    real(4) :: dsup1, dsup2, dsup05, xx, ss, Inu, u0, Psi, kx, ky
    real(8) :: ds(3), dt, ex ! important
    real(4), dimension(4), target :: xkernel1, ykernel1, xkernel2, ykernel2
    real(4), dimension(:), pointer, contiguous :: xkern, ykern
    real(4), dimension(:,:), allocatable :: work_up, work_min
    real(4), dimension(:,:,:,:), allocatable :: row_red
    integer :: xup_t(4,nx), yup_t(4,ny)
    integer :: zs, ze, xup, yup, zup, x1, y1, i, ii, jj, upo1, upo2
    logical :: inner, amr_boundary, rebuild_tables, get_rates
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%short_char', itimer=itimer)
    allocate(work_up(self%fs:self%fe,3))
    allocate(work_min(self%fs:self%fe,3))
    allocate(row_red(self%fs:self%fe, 3, nx, ny))
    get_rates = any(self%get_rates)

    zs = merge(1, nz, self%dir > 0)
    ze = merge(nz, 1, self%dir > 0)

    call ray%set_stencil(self%ds * scaling%l, vec)
    xkernel1 = real(ray%xkern)
    ykernel1 = real(ray%ykern)
    dsup1 = 0.5*real(ray%dsup)

    amr_boundary = amr/='' .and. self%dir < 0
    if (amr_boundary) then
      ds = self%ds
      if (associated(self%upstream)) ds = self%upstream%ds
      call amr_ray%set_stencil(ds * scaling%l, vec)
      xkernel2 = real(amr_ray%xkern)
      ykernel2 = real(amr_ray%ykern)
      dsup2 = 0.5*real(amr_ray%dsup)
    endif

    rebuild_tables = .true.
    do z=zs, ze, self%dir
      zup = z - self%dir

      if (amr_boundary .and. z==zs) then
        xkern => xkernel2
        ykern => ykernel2
        dsup05 = dsup2
        upo1 = amr_ray%upo1 - 2
        upo2 = amr_ray%upo2 - 2
      else
        xkern => xkernel1
        ykern => ykernel1
        dsup05 = dsup1
        upo1 = ray%upo1 - 2
        upo2 = ray%upo2 - 2
      endif

      if (rebuild_tables) then
        do x=1,nx; do ii=1,4
          xup_t(ii,x) = modulo(ii + x + upo1, nx) + 1
        end do; end do
        do y=1,ny; do jj=1,4
          yup_t(jj,y) = modulo(jj + y + upo2, ny) + 1
        end do; end do
        rebuild_tables = amr_boundary .and. z==zs
      endif

      !--- ROW PASS: x-reduce the whole upstream plane ----------------------------
      row_red = 0.0
      do yup = 1, ny
        do x = 1, nx
          do ii = 1, 4
            xup = xup_t(ii,x)
            kx = xkern(ii)
            do i=1, 3
              do f = self%fs, self%fe
                row_red(f,i,x,yup) = row_red(f,i,x,yup) + kx*work(f,xup,yup,zup,i)
              end do
            end do
          end do
        end do
      end do

      do y=1, ny
        do x=1, nx
          work_up = 0.0
          work_min = huge(1.0)

          ! 4x4 kernel sum
          do jj=1,4
            yup = yup_t(jj,y)
            ky = ykern(jj)
            do i=1, 3
              do f=self%fs,self%fe
                work_up(f,i) = work_up(f,i) + ky*row_red(f,i,x,yup)
              end do
            end do
          end do

          ! 2x2 min over the inner cells
          do jj=2,3
            do ii=2,3
              xup = xup_t(ii,x)
              yup = yup_t(jj,y)
              do i=1, 3
                do f=self%fs,self%fe
                  work_min(f,i) = min(work_min(f,i), work(f,xup,yup,zup,i))
                end do
              end do
            end do
          end do

          work_up = max(work_up, work_min)

          do f=self%fs, self%fe
            xx = work(f,x,y,z,1)
            ss = work(f,x,y,z,2)
            xxup  = work_up(f,1)
            ssup  = work_up(f,2)
            Inuup = work_up(f,3)

            dt = dsup05 * ( xxup + xx ) ! upwind optical pathlength
            if (dt < 80) then
              ex = exp(-dt)
              if ( dt > 1d-3 ) then
                u0 = 1 - ex
                Psi = 1 - u0 / dt
              else
                Psi = dt * (0.5 - 0.16666667*dt)
                u0 = dt * (1 - Psi)
              endif

              Inu = Inuup*ex + (u0 - Psi)*ssup + Psi*ss
            else
              !ex = 0
              u0 = 1 / dt
              Psi = 1 - u0
              Inu = u0*ssup + Psi*ss
            endif

            work(f,x,y,z,3) = Inu
            if (get_rates) self%Psi(f,x,y,z) = Psi / work(f,x,y,z,1)
          end do ! f-loop
        end do
      end do
    end do ! z-loop
    call trace%end(itimer)
  END SUBROUTINE short_char

  !===============================================================================
  SUBROUTINE long_char_linear_spherical (self, vec)
    class(rt_patch_t) :: self
    real(8) :: vec(3)
    real(4) :: ds, Inuup, xxup, ssup, dsup05, dsup1, dsup2, xx, ss, Inu, u0, Psi
    real(8) :: dt, ex ! important
    real(4) :: mu, p, lower, ds1, ds2, r1, r2, s1, s2
    integer :: dx, dy, zs, ze, zup
    logical :: amr_boundary
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%long_char_linear_spherical', itimer=itimer)

    p = rstar * sqrt(1.d0 - vec(3)**2) ! impact parameter
    lower = (self%position(3) - 0.5*self%size(3))*scaling%l + rorigin
    ds1 = self%ds(3) * scaling%l

    amr_boundary = amr/='' .and. self%dir < 0 .and. associated(self%upstream)
    if (amr_boundary) then
      ds2 = self%upstream%ds(3) * scaling%l
    endif

    zs = merge(1, nz, self%dir > 0)
    ze = merge(nz, 1, self%dir > 0)

    do z=zs, ze, self%dir
      zup = z - self%dir
      r2 = lower + z*ds1

      if (amr_boundary) then
        r1 = r2 - self%dir*ds2
        amr_boundary = .false.
      else
        r1 = r2 - self%dir*ds1
      endif

      if (r2 < p .or. r1 < p) then
        work(:,:,:,z,3) = 0
        cycle
      endif

      s1 = sqrt(r1**2 - p**2)
      s2 = sqrt(r2**2 - p**2)
      dsup05 = 0.5d0 * abs(s2-s1)

      do y=1, ny
        do x=1, nx
          do f=self%fs, self%fe
            xxup = work(f,x,y,zup,1)
            ssup = work(f,x,y,zup,2)
            xx = work(f,x,y,z,1)
            ss = work(f,x,y,z,2)
            dt = dsup05 * ( xxup + xx ) ! upwind optical pathlength
            dt = max(dt, 0.d0)
            if (dt < 80) then
              Inuup = work(f,x,y,zup,3)
              ex = exp(-dt)
              if ( dt > 1d-3 ) then
                u0 = 1 - ex
                Psi = 1 - u0 / dt
              else
                Psi = dt * (0.5 - 0.16666667*dt)
                u0 = dt * (1 - Psi)
              endif

              Inu = Inuup*ex + (u0 - Psi)*ssup + Psi*ss
            else
              !ex = 0
              u0 = 1 / dt
              Psi = 1 - u0
              Inu = u0*ssup + Psi*ss
            endif

            work(f,x,y,z,3) = Inu
            if (any(self%get_rates)) then
              self%Psi(f,x,y,z) = Psi / xx
            endif
          end do ! f-loop
        end do
      end do
    end do ! z-loop
    call trace%end(itimer)
  END SUBROUTINE long_char_linear_spherical

  !===============================================================================
  SUBROUTINE long_char_linear (self, vec)
    class(rt_patch_t) :: self
    real(8) :: vec(3)
    real(4) :: ds, Inuup, xxup, ssup, dsup05, dsup1, dsup2, xx, ss, Inu, u0, Psi
    real(8) :: dt, ex ! important
    integer :: dx, dy, zs, ze, zup
    logical :: amr_boundary
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%long_char_linear', itimer=itimer)

    ds = self%ds(3) * scaling%l
    dsup1 = 0.5*ds/abs(vec(3))

    amr_boundary = amr/='' .and. self%dir < 0
    if (amr_boundary) then
      if (associated(self%upstream)) ds = self%upstream%ds(3) * scaling%l
      dsup2 = 0.5*ds/abs(vec(3))
    endif

    zs = merge(1, nz, self%dir > 0)
    ze = merge(nz, 1, self%dir > 0)

    do z=zs, ze, self%dir
      zup = z - self%dir

      if (amr_boundary) then
        dsup05 = dsup2
        amr_boundary = .false.
      else
        dsup05 = dsup1
      endif

      do y=1, ny
        do x=1, nx
          do f=self%fs, self%fe
            xxup = work(f,x,y,zup,1)
            ssup = work(f,x,y,zup,2)
            xx = work(f,x,y,z,1)
            ss = work(f,x,y,z,2)
            dt = dsup05 * ( xxup + xx ) ! upwind optical pathlength
            dt = max(dt, 0.d0)
            if (dt < 80) then
              Inuup = work(f,x,y,zup,3)
              ex = exp(-dt)
              if ( dt > 1d-3 ) then
                u0 = 1 - ex
                Psi = 1 - u0 / dt
              else
                Psi = dt * (0.5 - 0.16666667*dt)
                u0 = dt * (1 - Psi)
              endif

              Inu = Inuup*ex + (u0 - Psi)*ssup + Psi*ss
            else
              !ex = 0
              u0 = 1 / dt
              Psi = 1 - u0
              Inu = u0*ssup + Psi*ss
            endif

            work(f,x,y,z,3) = Inu
            if (any(self%get_rates)) then
              !Psi = Psi + Inuup*ex*exp(-dsup05*xx*ss)*dsup05*xx
              self%Psi(f,x,y,z) = Psi / xx
            endif
          end do ! f-loop
        end do
      end do
    end do ! z-loop
    call trace%end(itimer)
  END SUBROUTINE long_char_linear

  !===============================================================================
  SUBROUTINE long_char_quadratic (self, vec)
    class(rt_patch_t) :: self
    real(8) :: vec(3)
    real(4) :: ds, dsup, rdsup, dsup1, dsup2
    real(4) :: xx, ss, Inu, Psi
    real(4) :: Inuup, xxup, ssup
    real(4) :: xxdn, ssdn
    real(4) :: t1, t2, t3, dtup, dtdn, dt2, dt3, U0, U1, U2, a, b, c
    real(8) :: dt, ex ! important
    integer :: dx, dy, zs, ze, z0, zup, zdn
    logical :: edge
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%long_char_quadratic', itimer=itimer)

    ds = self%ds(3) * scaling%l
    dsup = ds / abs(vec(3))
    rdsup = 1./ dsup
    dsup1 = 0.5*dsup
    dsup2 = 1./6.*dsup*dsup

    zs = merge(1, nz, self%dir > 0)
    ze = merge(nz, 1, self%dir > 0)

    do z=zs, ze, self%dir
      edge = z==ze
      z0 = merge(z-1, z, edge)
      zup = z0 - self%dir
      zdn = z0 + self%dir
      do y=1, ny
        do x=1, nx
          do f=self%fs, self%fe
            xx = work(f,x,y,z0,1)
            xxup = work(f,x,y,zup,1)
            xxdn = work(f,x,y,zdn,1)
            Inuup = work(f,x,y,z-self%dir,3)

            ! one sided dchi/ds at upwind point
            t1 = (xx - xxup) * rdsup
            ! one sided dchi/ds at dnwind point
            t3 = (xxdn - xx) * rdsup

            t2 = 0.5*(t1 + t3)
            t2 = (sign(1.0, t1) + sign(1.0, t3)) * min(abs(t1), abs(t3), abs(0.5*t2))

            dtup = dsup1*(xxup + xx) + dsup2*(t1-t2)
            dtdn = dsup1*(xxdn + xx) + dsup2*(t2-t3)

            dt = merge(dtdn, dtup, edge)
            dt2 = dt*dt
            ex=exp(-dt)

            ! Compute the U functions, use Taylor expansion for small dt
            ! to avoid accuracy issues with 1-(number slightly smaller than one)
            if ( dt > 1e-3 ) then
              U0 = 1.0d0 - ex
              U1 = dt - U0
              U2 = dt2 - 2.0d0 * U1
            else
              dt3 = dt2*dt
              U1 = half * dt2 - one_sixth * dt3
              U0 = dt - U1
              U2 = one_third * dt3
            end if

            ss = work(f,x,y,z0,2)
            ssdn = work(f,x,y,zdn,2)
            ssup = work(f,x,y,zup,2)

            ! one sided dS/dt at upwind point
            t1 = (ss - ssup) / dtup
            t3 = (ssdn - ss) / dtdn

            t2 = 0.5*(t1 + t3)
            t2 = (sign(1.0, t1) + sign(1.0, t3)) * min(abs(t1), abs(t3), abs(0.5*t2))

            if (edge) then
              a = (t3-t2) / dt
              b = 2*t2 - t3
              c = ss
            else
              a = (t2-t1) / dt
              b = 2*t1 - t2
              c = ssup
            endif

            ! intensity at center point
            Inu = Inuup*ex + a*U2 + b*U1 + c*U0
            work(f,x,y,z,3) = Inu
          end do ! f-loop
        end do
      end do
    end do ! z-loop
    call trace%end(itimer)
  END SUBROUTINE long_char_quadratic

  !===============================================================================
  SUBROUTINE long_char_quadratic2 (self, vec)
    class(rt_patch_t) :: self
    real(8) :: vec(3)
    real(4) :: ds, dsup, rdsup, dsup1, dsup2
    real(4) :: xx, ss, Inu, Psi
    real(4) :: Inuup, xxup, ssup
    real(4) :: xxdn, ssdn
    real(4) :: t1, t2, t3, dtup, dtdn, dt2, dt3, U0, U1, U2, a, b, c
    real(8) :: dt, ex ! important
    integer :: dx, dy, zs, ze, z0, zup, zdn
    logical :: edge
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%long_char_quadratic2', itimer=itimer)

    ds = self%ds(3) * scaling%l
    dsup = ds / abs(vec(3))
    rdsup = 1./ dsup
    dsup1 = 0.5*dsup
    dsup2 = 1./6.*dsup*dsup

    zs = merge(1, nz, self%dir > 0)
    ze = merge(nz, 1, self%dir > 0)

    do z=zs, ze, self%dir
      edge = z==ze
      z0 = merge(z-1, z, edge)
      zup = z0 - self%dir
      zdn = z0 + self%dir
      do y=1, ny
        do x=1, nx
          do f=self%fs, self%fe
            xx = work(f,x,y,z0,1)
            xxup = work(f,x,y,zup,1)
            xxdn = work(f,x,y,zdn,1)
            Inuup = work(f,x,y,z-self%dir,3)

            ! one sided dchi/ds at upwind point
            t1 = (xx - xxup) * rdsup
            ! one sided dchi/ds at dnwind point
            t3 = (xxdn - xx) * rdsup

            t2 = 0.5*(t1 + t3)
            t2 = (sign(1.0, t1) + sign(1.0, t3)) * min(abs(t1), abs(t3), abs(0.5*t2))

            dtup = dsup1*(xxup + xx) + dsup2*(t1-t2)
            dtdn = dsup1*(xxdn + xx) + dsup2*(t2-t3)

            dt = merge(dtdn, dtup, edge)
            dt2 = dt*dt
            ex=exp(-dt)

            ! Compute the U functions, use Taylor expansion for small dt
            ! to avoid accuracy issues with 1-(number slightly smaller than one)
            if ( dt > 1e-3 ) then
              U0 = 1.0d0 - ex
              U1 = dt - U0
              U2 = dt2 - 2.0d0 * U1
            else
              dt3 = dt2*dt
              U1 = half * dt2 - one_sixth * dt3
              U0 = dt - U1
              U2 = one_third * dt3
            end if

            ss = work(f,x,y,z0,2)
            ssdn = work(f,x,y,zdn,2)
            ssup = work(f,x,y,zup,2)

            ! one sided dS/dt at upwind point
            t1 = ss - ssup
            t3 = ssdn - ss

            t2 = 0.5*(t1 + t3)
            t2 = (sign(1.0, t1) + sign(1.0, t3)) * min(abs(t1), abs(t3), abs(0.5*t2))

            if (edge) then
              t3 = t3 / dtdn
              t2 = t2 * rdsup / xx
              a = (t3-t2) / dtdn
              b = 2*t2 - t3
              c = ss
            else
              t1 = t1 / dtup
              t2 = t2 * rdsup / xx
              a = (t2-t1) / dtup
              b = 2*t1 - t2
              c = ssup
            endif

            ! intensity at center point
            Inu = Inuup*ex + a*U2 + b*U1 + c*U0
            work(f,x,y,z,3) = Inu
          end do ! f-loop
        end do
      end do
    end do ! z-loop
    call trace%end(itimer)
  END SUBROUTINE long_char_quadratic2

  !===============================================================================
  SUBROUTINE long_char_quadratic3 (self, comp, vec)
    class(rt_patch_t) :: self
    class(composition_type) :: comp
    real(8), dimension(:,:,:), allocatable :: chi500, dchi500, Bnu500, dBnu500
    real(8), dimension(:) :: xxold(3), xxnew(2)
    real(8), dimension(:), allocatable :: ds3
    real(8) :: vec(3), hnu
    real(4) :: ds, dsup, rdsup, dsup1, dsup2
    real(4) :: xx, ss, Inu, Psi, beta
    real(4) :: Inuup, xxup, ssup
    real(4) :: xxdn, ssdn
    real(4) :: xx500, xxup500, dt500
    real(4) :: ss500, ssup500, dS500, dSup500
    real(4) :: t1, t2, t3, dtup, dtdn, U0, U1, U2, a, b, c, dt2, dt3, tt2
    real(8) :: dt, ex ! important
    integer :: dx, dy, zs, ze, z0, zup, zdn
    integer :: i, nz_up
    logical :: edge
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%long_char_quadratic3', itimer=itimer)

    ds = self%ds(3) * scaling%l
    dsup = ds / abs(vec(3))
    dsup1 = 0.5*dsup
    dsup2 = 1./12.*dsup*dsup

    zs = merge(1, nz, self%dir > 0)
    ze = merge(nz, 1, self%dir > 0)

    allocate(ds3(0:nz+1))
    allocate(chi500(nx,ny,-1:nz+1), dchi500(nx,ny,0:nz))
    allocate(Bnu500(nx,ny,-1:nz+1), dBnu500(nx,ny,0:nz))
    ds3(1:nz) = self%ds(3)
    chi500(:,:,1:nz) = comp%chi500

    hnu = h_k * CC / air2vac(ref_lam) * 1e8
    Bnu500(:,:,1:nz) = exp(-hnu / comp%temp)
    Bnu500(:,:,1:nz) = Bnu500(:,:,1:nz) / (1 - Bnu500(:,:,1:nz))

    if (associated(self%upstream)) then
      nz_up = size(self%upstream%composition%chi500, 3)
      ds3(0) = self%upstream%ds(3)
      chi500(:,:,-1:0) = self%upstream%composition%chi500(:,:,nz_up-1:nz_up)
      Bnu500(:,:,-1:0) = exp(-hnu / self%upstream%composition%temp(:,:,nz_up-1:nz_up))
      Bnu500(:,:,-1:0) = Bnu500(:,:,-1:0) / (1 - Bnu500(:,:,-1:0))
    else
      ds3(0) = self%ds(3)
      xxold = [(i*self%ds(3), i= 1, 3)]
      xxnew = [(i*self%ds(3), i=-1, 0)]
      do y=1, ny
        do x=1, nx
          call pchip_mono(xxold, chi500(x,y,1:3), xxnew, chi500(x,y,-1:0), sorted=.true.)
          call pchip_mono(xxold, Bnu500(x,y,1:3), xxnew, Bnu500(x,y,-1:0), sorted=.true.)
        end do
      end do
    endif
    if (associated(self%dnstream)) then
      ds3(nz+1) = self%dnstream%ds(3)
      chi500(:,:,nz+1) = self%dnstream%m3d_parent%composition%chi500(:,:,1)
      Bnu500(:,:,nz+1) = exp(-hnu / self%dnstream%m3d_parent%composition%temp(:,:,1))
      Bnu500(:,:,nz+1) = Bnu500(:,:,nz+1) / (1 - Bnu500(:,:,nz+1))
    else
      ds3(nz+1) = self%ds(3)
      xxold = [(i*self%ds(3), i=1, 3)]
      xxnew(1:1) = [(i*self%ds(3), i=4, 4)]
      do y=1, ny
        do x=1, nx
          call pchip_mono(xxold, chi500(x,y,nz-2:nz), xxnew(1:1), chi500(x,y,nz+1:nz+1), sorted=.true.)
          call pchip_mono(xxold, Bnu500(x,y,nz-2:nz), xxnew(1:1), Bnu500(x,y,nz+1:nz+1), sorted=.true.)
        end do
      end do
    endif

    ds3 = ds3 * (scaling%l / abs(vec(3)))

    do z=0, nz
      do y=1, ny
        do x=1, nx
          xx = chi500(x,y,z)
          xxup = chi500(x,y,z-1)
          xxdn = chi500(x,y,z+1)
          t1 = (xx - xxup) / ds3(z)
          t3 = (xxdn - xx) / ds3(z+1)
          t2 = 0.5*(t1 + t3)
          t2 = (sign(1.0, t1) + sign(1.0, t3)) * min(abs(t1), abs(t3), abs(0.5*t2))
          dchi500(x,y,z) = t2

          ss = Bnu500(x,y,z)
          ssup = Bnu500(x,y,z-1)
          ssdn = Bnu500(x,y,z+1)
          t1 = (ss - ssup) / ds3(z)
          t3 = (ssdn - ss) / ds3(z+1)
          t2 = 0.5*(t1 + t3)
          t2 = (sign(1.0, t1) + sign(1.0, t3)) * min(abs(t1), abs(t3), abs(0.5*t2))
          dBnu500(x,y,z) = t2
        end do
      end do
    end do

    do z=zs, ze, self%dir
      zup = z - self%dir
      zdn = z + self%dir
      do y=1, ny
        do x=1, nx
          xx500   = chi500(x,y,z)
          xxup500 = chi500(x,y,zup)
          ss500   = Bnu500(x,y,z)
          ssup500 = Bnu500(x,y,zup)
          dS500   = dBnu500(x,y,z)
          dSup500 = dBnu500(x,y,zup)
          t1      = dchi500(x,y,zup)
          t2      = dchi500(x,y,z)

          dt500 = dsup1*(xxup500 + xx500) + dsup2*(t1-t2)
          dt500 = 0.5*dt500
          xxup500 = dt500 / xxup500
          xx500   = dt500 / xx500
          dS500   = (ss500 - ssup500) / dS500
          dSup500 = (ss500 - ssup500) / dSup500

          do f=self%fs, self%fe
            xx   = work(f,x,y,z,1)
            xxup = work(f,x,y,zup,1)
            ss   = work(f,x,y,z,2)
            ssup = work(f,x,y,zup,2)
            Inuup = work(f,x,y,z-self%dir,3)

            dt = xxup*xxup500 + xx*xx500
            ex=exp(-dt)
            ! Compute the U functions, use Taylor expansion for small dt
            ! to avoid accuracy issues with 1-(number slightly smaller than one)
            !dt2 = dt*dt
            !if ( dt > 1e-3 ) then
            !  U0 = 1.0d0 - ex
            !  U1 = dt - U0
            !  U2 = dt2 - 2.0d0 * U1
            !else
            !  dt3 = dt2*dt
            !  U1 = half * dt2 - one_sixth * dt3
            !  U0 = dt - U1
            !  U2 = one_third * dt3
            !end if
            !
            !t1 = ss - ssup
            !t2 = t1 / (xx*dS500)
            !t1 = t1 / dt
            !
            !a = (t2-t1) / dt
            !b = 2*t1 - t2
            !c = ssup
            !
            !! intensity at center point
            !Inu = Inuup*ex + a*U2 + b*U1 + c*U0
            !work(f,x,y,z,3) = Inu

            if (dt > 1e-3) then
              U0 = 1. - ex
              U1 = U0 - dt*ex
            else
              dt2 = 0.5*dt*dt
              U1 = dt2 * (1. - one_sixth*dt)
              U0 = dt - U1 + dt2*dt
            endif

            t1 = (ss - ssup) / dt
            t2 = (ss - ssup) / (xx*dS500)
            tt2 = 2*(t2 - t1) / dt

            ! intensity at center point
            Inu = (Inuup - ssup)*ex - U0*t2 + U1*tt2 + ss
            work(f,x,y,z,3) = Inu
          end do ! f-loop
        end do
      end do
    end do ! z-loop
    call trace%end(itimer)
  END SUBROUTINE long_char_quadratic3

  !===============================================================================
  SUBROUTINE get_vproj (self, vproj, vec, rotate, i_mu)
    class(rt_patch_t) :: self
    real(4), dimension(:,:,:) :: vproj
    real(8), dimension(:) :: vec
    real(4), dimension(:,:,:), pointer :: vx, vy, vz
    real(4), dimension(:,:,:,:), contiguous, pointer :: mem
    integer :: i_mu
    logical :: rotate
    integer, save :: itimer
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch_t%get_vproj', itimer=itimer)
    if (rotate) then
      mem => self%m3d_parent%mem(:,:,:,:,self%it,1+i_mu)
    else
      mem => self%m3d_parent%mem(:,:,:,:,self%it,1)
    endif

    vx => mem(:,:,:,self%i_vx)
    vy => mem(:,:,:,self%i_vy)
    vz => mem(:,:,:,self%i_vz)

    if (.not. use_vel(1)) vx = 0
    if (.not. use_vel(2)) vy = 0
    if (.not. use_vel(3)) vz = 0

    vproj = vec(1)*vx + vec(2)*vy + vec(3)*vz
    call trace%end(itimer)
  END SUBROUTINE get_vproj

  !===============================================================================
  SUBROUTINE save_contribution_function (self, i_set, i_mu)
    class(rt_patch_t) :: self
    integer :: i, j, i_set, i_mu, ff, im
    real(4) :: net_chi, net_src, lin_chi, lin_eta
    real(4) :: cnt_inu, cnt_chi, cnt_eta
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch_t%save_contribution_function', itimer=itimer)
    do z=1,nz
      do y=1,ny
        do x=1,nx
          do f=self%fs,self%fe
            associate(cf => cntrbf_tmp(f,x,y,z,:), cf_out => cntrbf_data(f,x,y,z,:,i_mu))
              if (i_set == 1) then ! NLTE set
                net_chi = work(f,x,y,z,1)
                net_src = work(f,x,y,z,2)
                cf(1) = net_chi
                cf(2) = merge(net_src*net_chi, net_src, depression_cntrbf)
              elseif (i_set == 2) then ! continuum set
                cnt_chi = work(f,x,y,z,1)
                cnt_eta = work(f,x,y,z,2)*cnt_chi
                cnt_inu = work(f,x,y,z,3)

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
      call self%save6(cntrbf_data, 'cntrbf')
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
    if (ilambd > 0) then
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
    ! ----------------------------------------------------------------------------
    call self%m3d_patch_t%finish()
    if (associated(self%upstream)) then
      if (associated(self%upstream%mem)) then
        call self%upstream%deallocate_mem()
      endif
    endif
    if (ilambd > 0 .and. (self%do_full_sp .or. self%do_final)) then
      if (self%up_twin%is_finished() .and. self%dn_twin%is_finished()) then
        if (Jnu_to_file) then
          if (.not. save_Jnu) then
            Jnu_label = merge('Jnu_out', 'Jnu    ', self%do_out_sp)
            call self%up_twin%remove(Jnu_label, key=self%up_twin%key(1:3))
            call self%dn_twin%remove(Jnu_label, key=self%dn_twin%key(1:3))
          endif
        else
          call mem_unregister(self%up_twin%Jnu, label='Jnu')
          call mem_unregister(self%dn_twin%Jnu, label='Jnu')
          deallocate(self%up_twin%Jnu)
          deallocate(self%dn_twin%Jnu)
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
      nx = self%n(1)
      ny = self%n(2)
      tmp => self%mem
      self%mem(fs:fe,1:nx,1:ny,self%mus:self%mue,1:self%n_set,1:3) => tmp
    endif
  END SUBROUTINE allocate_mem

  !===============================================================================
  SUBROUTINE save_intensity (self)
    class(rt_patch_t) :: self
    real(4), dimension(:,:,:,:,:), allocatable :: ie, i3
    integer :: fs, fe
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch_t%save_intensity', itimer=itimer)
    fs = self%fs
    fe = self%fe
    allocate(i3(fs:fe,nx,ny, self%mus:self%mue, self%n_set), source=0.)
    allocate(ie(fs:fe, 1, 1, self%mus:self%mue, self%n_set), source=0.)

    do f=fs, fe
      i3(f,:,:,:,:) = self%twohnu3_c2(f) * self%mem(f,:,:,:,:,3)
    end do

    do y=1, ny
      do x=1, nx
        ie(:,1,1,:,:) = ie(:,1,1,:,:) + i3(:,x,y,:,:)
      end do
    end do
    ie = ie / real(nx*ny)

    call self%save5(ie, 'ie')

    if (save_resolved) then
      call self%save5(i3, 'i3')
    endif
    call trace%end(itimer)
  END SUBROUTINE save_intensity

  !===============================================================================
  SUBROUTINE copy_to (self, tgt)
    class(rt_patch_t):: self
    class(patch_t):: tgt
    integer :: i, i_set, i_mu
    !-----------------------------------------------------------------------------
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
