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
  USE interpol_mod, only: pchip_mono
  USE handyfuncs_mod
  USE atom_mod, only: atom_type, atom_data_t, atom_trans_t, atomic_spectrum_type
  USE ray_patch_mod, only: ray_patch_t, output_patch_window, output_patch_range
  USE m3d_halo_mod, only: lock_mem, unlock_mem
  USE m3d_patch_mod, only: use_vel, m3d_patch_t, point_opac_set_t, &
    dlam, n_nu_out, &
    premature_finish, m3d_patch_pre_update=>pre_update, save_resolved, &
    column_by_column, rotate_atmos, save_cntrbf, depression_cntrbf, &
    decouple_continuum, initial_Jnu,  &
    Jnu_to_file, preconditioning, scattering_preconditioning, scat_op_max, &
    make_eos, amr, save_Inu, save_chi, save_Snu, &
    fixed_Jnu, n_nu, save_Jnu, &
    convlim, converged, flux_convergence, flux_conv_signal, &
    global_flux_err, global_flux_pos, global_flux_lambda, flux_update_counter, &
    maser_chi_floor, maser_dtau_max, maser_Inu_max, suppress_masers, set_file
  USE mpi_mod, only: mpi
  USE rt_params_mod, only: rt_params_t, rt_iter, rt_out, ray_offset, pivot_z, &
                           save_linelist
  USE ray_mod, only: one_dimensional, ray_type, rstar, rorigin
  USE trans_mod, only: line_type, cont_type
  USE bbtrans_mod, only: tolerance
  USE qlg, only: qlg_bot_profile      ! the one Voigt kernel; see linelist_mod's kernel note
  USE linelist_mod, only: ll_add_opacity => add_opacity, &
    ll_add_opacity_equi => add_opacity_equi, voigt_peak
  USE element_mod, only: element_type
  USE molecule_mod, only: molec_type
  USE ion_mod, only: particle_type
  USE bftrans_mod, only: H_occupation_probability
  USE hprofl_mod, only: hprofile
  USE m3d_debug_mod
  USE shared_mod
  USE composition_mod, only: composition_type, ref_lam, &
    departure_pointer_t, read_departures, &
    line_scattering
  USE osmet_mod, only: compute_epsilon_stage
  USE spectrum_mod, only: spectrum_type
  USE optab_mod, only: sample_optab_type
  USE omp_lock_mod, only: lock_t
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer :: n_negative_chi = 0                 ! diagnostic: cells with chi<=0 in work_to_transfer_units

  !-----------------------------------------------------------------------------
  ! What a set *is*, as opposed to where it sits.  Historically the meaning of a
  ! set was carried by its index alone (set 1 = main, 2 = continuum, 3 = forced
  ! LTE), which stops working the moment the set axis becomes variable-length --
  ! and is already wrong for the ALI branch below, whose single set is the main
  ! one and has no continuum companion.  These kinds are recorded per set and
  ! written to 'sets.txt' so the reader looks the layout up instead of guessing
  ! it from the component count.
  !-----------------------------------------------------------------------------
  integer, parameter, public :: set_continuum=0    ! background only; skips foreground
  integer, parameter, public :: set_nlte=1         ! the main spectrum, with departures
  integer, parameter, public :: set_lte=2          ! LTE: the main spectrum of an LTE
                                                   ! run, or the forced-LTE companion
                                                   ! of an NLTE one
  character(len=9), dimension(0:2), parameter :: set_kind_name = &
    ['continuum', 'nlte     ', 'lte      ']
  integer, public, save :: output_n_set=0
  integer, public, save :: output_i_set_main=0, output_i_set_cnt=0
  integer, allocatable, public, save :: output_set_kind(:)
  real(8), allocatable, public, save :: output_set_da(:,:) ! (atomic number, set)
  logical, allocatable, public, save :: output_set_element(:,:) ! explicit rows, including da=0

  !-----------------------------------------------------------------------------
  ! Output patches use IDs -n_nu_out:-1, internal patches use 1:n_nu.  n_nu_out
  ! (m3d_patch_mod, set in composition_pre_init) equals n_nu for a contiguous output
  ! grid; a window_file grid gets a whole number of patches per window.
  !-----------------------------------------------------------------------------
  public :: verify_output_patch_coverage

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
    real(4), dimension(:,:,:,:), allocatable :: Psi
    ! PsiJ is the angle-averaged diagonal operator and Jali the full-sphere
    ! accelerated mean intensity of the continuum-scattering ALI.  Both carry a
    ! set index (the output-spectrum path has n_set>1) and both are persisted
    ! exactly like Jnu, because the OAB update is a two-term recurrence that
    ! needs the previous iteration's values at the twin-combine point.
    ! Jali_valid: per set, has a previous Jali been stored yet (else the
    ! recurrence is seeded with the plain formal solution).
    ! do_scat_precond: is the correction active for this patch.
    ! compute_psi: must the formal solvers build Psi at all (rates OR ALI).
    real(4), dimension(:,:,:,:,:), allocatable :: PsiJ, Jali
    logical, allocatable :: Jali_valid(:)
    logical :: do_scat_precond = .false.
    logical :: compute_psi = .false.
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
    integer, dimension(:), allocatable :: set_kind   ! set_continuum / set_nlte / set_lte
    integer :: i_set_main=1     ! the set the run's own spectrum comes from
    integer :: i_set_cnt=0      ! the continuum set, or 0 where none is computed
    logical :: has_nlte_set=.false.   ! is there a set carrying departures?
    logical :: bound_bound, bound_free
    logical, dimension(:), allocatable :: get_rates
    integer :: nx, ny, nz
    real(4), dimension(:,:,:,:,:), allocatable :: work
    ! Direction-independent background, per (f,x,y,z), rebuilt once per i_set in
    ! background().  Filled in this order:
    !   lte_background_opacities  chi_cont, sigma  (from the table or comp%opacity)
    !                             beta = exp(-h*nu/kT), clamped off zero
    !   lte_background_emission   ehat = chi_cont*Bnu + sigma*(Jnu - Bnu)   [= eta]
    !   convert_to_work_units     xhat = chi_cont + ehat, then ehat -> ehat/beta
    !   add_bound_free_opacities  accumulates into xhat and ehat
    ! so xhat/ehat hold the emissivity eta in *physical* units until
    ! convert_to_work_units runs, and the working pair afterwards.  back_to_work
    ! then seeds work(...,1:2) from xhat/ehat for each direction.
    real(4), dimension(:,:,:,:), allocatable :: chi_cont  ! continuum extinction, total
    real(4), dimension(:,:,:,:), allocatable :: sigma     ! coherent-scattering part of chi_cont (<= it)
    real(4), dimension(:,:,:,:), allocatable :: beta      ! LTE Boltzmann factor exp(-h*nu/kT)
    real(4), dimension(:,:,:,:), allocatable :: xhat      ! working extinction (see above)
    real(4), dimension(:,:,:,:), allocatable :: ehat      ! emissivity / working source (see above)
    ! line_scattering: per-(f,x,y,z) accumulator of the (1-eps)-weighted LTE
    ! line opacity built during foreground, and (direct path only) the per-stage
    ! scattering fraction 1-eps on the patch frequency grid (index 1=neutral,
    ! 2=ionized), refreshed per nu-chunk
    ! line_scattering scratch: chi_scat is the (1-eps)-weighted line opacity for
    ! the current direction, dJnu the (Jnu - Bnu) the background left behind, and
    ! do_line_scat says whether this pass filled them (it follows the per-pass
    ! include_linelist / sample_linelist, so it can differ between rt_iter and
    ! rt_out while the buffers stay allocated)
    logical :: do_line_scat=.false.
    real(4), dimension(:,:,:,:),   allocatable :: chi_scat, dJnu
    real(4), dimension(:,:,:,:,:), allocatable :: one_m_eps
    real(4), dimension(:,:,:,:,:,:), allocatable :: cntrbf_data
    real(4), dimension(:,:,:,:,:,:), allocatable :: cntrbf_cnt   ! per-mu continuum stash
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
    procedure :: save_set_manifest
    procedure :: init_set_layout
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
    procedure :: work_to_transfer_units
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
    procedure :: free_ll_point_table
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
  ! NB: work, the background arrays, cntrbf_data and the patch
  ! dimensions nx,ny,nz are
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
      call self%init_set_layout(is_nlte)
    else
      ! For the ALI we do not need to compute the continuum or the LTE set
      self%n_set = 1
      ! ALI / scattering iteration: the single set is the main one, no continuum
      allocate(self%set_kind(1))
      self%set_kind(1) = merge(set_nlte, set_lte, is_nlte)
      self%i_set_cnt   = 0
      self%i_set_main  = 1
    endif
    self%has_nlte_set  = is_nlte

    ! Store (mu,set) as the generic variable axis and the three transfer fields
    ! as its generic w axis. ray_patch pack/unpack then traverses the complete
    ! boundary buffer; the set axis is not hidden in task_t's time slot.
    call self%set_nv_and_nw(self%n_mu*self%n_set, 3)
    self%nt = 1

    call self%set_freq_indices()
    call self%init_patch_trans()
    allocate(self%get_rates(self%natom))
    self%get_rates = self%do_full_sp .and. self%i_nu>0

    self%ng(3) = 1
    call self%post_pre_init()

    call trace%end(itimer)
  END SUBROUTINE rt_patch_pre_init

  !===============================================================================
  !> Derive the historical base layout, or read an abundance-ladder manifest.
  !> The process-global copy is also consumed once by each composition instance.
  !===============================================================================
  SUBROUTINE init_set_layout (self, is_nlte)
    class(rt_patch_t) :: self
    logical, intent(in) :: is_nlte
    integer :: unit, ios, number, max_set, z, i, ncnt, nmain
    character(len=1024) :: line
    character(len=16) :: kind, element, da_word
    real(8) :: da
    logical, allocatable :: seen(:), has_offset(:)
    logical, save :: initialised=.false.
    !-----------------------------------------------------------------------------
    !$omp critical (cr_init_set_layout)
    if (.not. initialised) then
      if (len_trim(set_file) == 0) then
        output_n_set = merge(3, 2, is_nlte)
        allocate(output_set_kind(output_n_set), source=0)
        allocate(output_set_da(92,output_n_set), source=0d0)
        allocate(output_set_element(92,output_n_set), source=.false.)
        output_set_kind(1) = set_continuum
        output_set_kind(2) = merge(set_nlte, set_lte, is_nlte)
        if (is_nlte) output_set_kind(3) = set_lte
      else
        if (index(trim(set_file), trim(io%run_dir)) == 1) then
          print '(a)', 'ERROR: output_params set_file must be outside the run directory'
          error stop 'invalid set_file location'
        endif
        open(newunit=unit, file=trim(set_file), status='old', action='read', iostat=ios)
        if (ios /= 0) then
          print '(2a)', 'ERROR: cannot open set_file: ', trim(set_file)
          error stop 'cannot open set_file'
        endif
        max_set = 0
        do
          read(unit, '(A)', iostat=ios) line
          if (ios /= 0) exit
          if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
          read(line, *, iostat=ios) number, kind, element, da_word
          if (ios /= 0 .or. number < 1) error stop 'malformed row in set_file'
          max_set = max(max_set, number)
        enddo
        if (max_set == 0) error stop 'set_file contains no sets'
        rewind(unit)
        allocate(output_set_kind(max_set), source=-1)
        allocate(output_set_da(92,max_set), source=0d0)
        allocate(output_set_element(92,max_set), source=.false.)
        allocate(seen(max_set), source=.false.)
        do
          read(unit, '(A)', iostat=ios) line
          if (ios /= 0) exit
          if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
          read(line, *, iostat=ios) number, kind, element, da_word
          if (ios /= 0 .or. number < 1 .or. number > max_set) &
            error stop 'malformed row in set_file'
          select case (trim(kind))
          case ('continuum'); i = set_continuum
          case ('nlte');      i = set_nlte
          case ('lte');       i = set_lte
          case default
            print '(2a)', 'ERROR: unknown kind in set_file: ', trim(kind)
            error stop 'invalid set kind'
          end select
          if (seen(number) .and. output_set_kind(number) /= i) &
            error stop 'one set number has conflicting kinds'
          output_set_kind(number) = i
          seen(number) = .true.
          if (trim(element) == '-') then
            if (trim(da_word) /= '-') error stop 'set_file has da without element'
          else
            z = atomic_number(element)
            if (z == 0) then
              print '(2a)', 'ERROR: unknown element in set_file: ', trim(element)
              error stop 'invalid ladder element'
            endif
            read(da_word, *, iostat=ios) da
            if (ios /= 0) error stop 'invalid abundance offset in set_file'
            if (output_set_element(z,number)) &
              error stop 'duplicate element in one set_file set'
            output_set_da(z,number) = da
            output_set_element(z,number) = .true.
          endif
        enddo
        close(unit)
        if (.not. all(seen)) error stop 'set numbers must be contiguous from one'
        output_n_set = max_set
        deallocate(seen)
      endif

      allocate(has_offset(output_n_set))
      has_offset = any(output_set_element, dim=1)
      ncnt = count(output_set_kind == set_continuum .and. .not. has_offset)
      nmain = count(output_set_kind == merge(set_nlte,set_lte,is_nlte) .and. .not. has_offset)
      if (ncnt /= 1) error stop 'set layout needs exactly one unperturbed continuum set'
      if (nmain /= 1) error stop 'set layout needs exactly one unperturbed main set'
      if (.not. is_nlte .and. any(output_set_kind == set_nlte)) &
        error stop 'nlte set is illegal without departure coefficients'
      if (any(output_set_kind < set_continuum) .or. any(output_set_kind > set_lte)) &
        error stop 'invalid set kind'
      output_i_set_cnt = findloc(output_set_kind == set_continuum .and. .not. has_offset, .true., dim=1)
      output_i_set_main = findloc(output_set_kind == merge(set_nlte,set_lte,is_nlte) .and. &
                                  .not. has_offset, .true., dim=1)
      deallocate(has_offset)
      initialised = .true.
    endif
    !$omp end critical (cr_init_set_layout)

    self%n_set = output_n_set
    allocate(self%set_kind(self%n_set), source=output_set_kind)
    self%i_set_cnt = output_i_set_cnt
    self%i_set_main = output_i_set_main
    if (self%n_set > merge(3,2,is_nlte)) then
      if (rt_out%sample_linelist) &
        error stop 'abundance ladders require rt_output_params sample_linelist=.false.'
    endif
  contains
    integer function atomic_number(symbol) result(number)
      character(len=*), intent(in) :: symbol
      character(len=3), parameter :: ids(92) = [ &
        'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ', &
        'Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr ', &
        'Mn ','Fe ','Co ','Ni ','Cu ','Zn ','Ga ','Ge ','As ','Se ','Br ','Kr ', &
        'Rb ','Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ', &
        'In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ', &
        'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ', &
        'Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po ', &
        'At ','Rn ','Fr ','Ra ','Ac ','Th ','Pa ','U  ']
      character(len=3) :: key
      key = adjustl(symbol)
      number = findloc(ids, key, dim=1)
    end function atomic_number
  END SUBROUTINE init_set_layout

  !===============================================================================
  !> Write 'sets.txt': one row per (set, element), naming what each slot on the
  !> saved set axis of 'ie' and 'flux' actually contains.  The element and offset
  !> columns are placeholders until per-set abundance offsets exist; a set with no
  !> offsets is one row with '-' in both.  Written once per run, by the first
  !> output patch to reach pre_init.
  !===============================================================================
  SUBROUTINE save_set_manifest (self)
    class(rt_patch_t) :: self
    integer :: lunit, i_set, z
    ! left-justified placeholder columns; a bare '-' under an A7 edit descriptor
    ! would right-justify and misalign against the header
    character(len=7) :: c_element, c_da
    logical, save :: first_time = .true.
    !-----------------------------------------------------------------------------
    !$omp critical (cr_save_set_manifest)
    if (first_time .and. mpi%master) then
      first_time = .false.
      c_element = '-'
      c_da      = '-'
      open(newunit=lunit, file=trim(io%run_dir)//'sets.txt', form='FORMATTED')
      write(lunit, '(A)') '# set  kind       element  da'
      do i_set=1, self%n_set
        c_element = '-'
        c_da = '-'
        if (any(output_set_element(:,i_set))) then
          do z=1,92
            if (.not. output_set_element(z,i_set)) cycle
            write(c_element, '(A)') trim(element_symbol(z))
            write(c_da, '(SP,F7.3)') output_set_da(z,i_set)
            write(lunit, '(I5, 2X, A9, 2X, A7, 2X, A7)') &
              i_set, set_kind_name(self%set_kind(i_set)), c_element, c_da
          enddo
        else
          write(lunit, '(I5, 2X, A9, 2X, A7, 2X, A7)') &
            i_set, set_kind_name(self%set_kind(i_set)), c_element, c_da
        endif
      end do
      close(lunit)
    endif
    !$omp end critical (cr_save_set_manifest)
  contains
    function element_symbol(z) result(symbol)
      integer, intent(in) :: z
      character(len=2) :: symbol
      character(len=2), parameter :: ids(92) = [ &
        'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar', &
        'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
        'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe', &
        'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf', &
        'Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
        'Pa','U ']
      symbol = ids(z)
    end function element_symbol
  END SUBROUTINE save_set_manifest

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
    ! The path covers opacity only (ll_add_opacity); there is no SoA rate kernel,
    ! so enabling it needs one written first -- otherwise its lines carry opacity
    ! but no radiative rates. bound_bound_rates stops rather than allow that.
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
    if (self%do_out_sp) call self%save_set_manifest()
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
    ! Continuum-scattering ALI.  It applies to the internal NLTE iteration
    ! (do_full_sp, where Psi is already built for the rates) and to the
    ! output-spectrum scattering loop (do_out_sp).  The do_final single pass has
    ! no iteration to accelerate and is left untouched.  The maxiter test is the
    ! same one get_Jnu uses: without it there is no Jnu to accelerate.
    self%do_scat_precond = scattering_preconditioning &
      .and. .not. self%do_final &
      .and. (any(self%get_rates) .or. self%do_out_sp) &
      .and. ((rt_iter%maxiter > 0) .or. (rt_out%maxiter > 0))
    ! Psi must be built for the rates OR for the scattering ALI, so the formal
    ! solvers gate its computation on compute_psi rather than on get_rates.
    self%compute_psi = any(self%get_rates) .or. self%do_scat_precond
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
          if (i_set /= self%i_set_cnt) &
            call self%foreground(comp, i_mu, i_set, do_rotate_atmos, wmu, vproj)
          call self%work_to_transfer_units()
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
          if (i_set /= self%i_set_cnt) &
            call self%foreground(comp, i_mu, i_set, do_rotate_atmos, wmu, vproj)
          call self%work_to_transfer_units()
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
            call self%work_to_transfer_units()
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
    if (self%do_final) call self%free_ll_point_table()

    if (end_boundary) then
      if (self%do_final .and. self%n_mu>0) then
        call self%save_intensity()
        if (save_linelist) call self%print_linelist(comp, self%n_nu)
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
      ! background arrays: see the component declarations for how they are filled
      allocate(self%chi_cont(fs:fe,self%nx,self%ny,self%nz), source=0.)
      allocate(self%sigma   (fs:fe,self%nx,self%ny,self%nz), source=0.)
      allocate(self%beta    (fs:fe,self%nx,self%ny,self%nz), source=0.)
      allocate(self%xhat    (fs:fe,self%nx,self%ny,self%nz), source=0.)
      allocate(self%ehat    (fs:fe,self%nx,self%ny,self%nz), source=0.)
      allocate(self%work(fs:fe,self%nx,self%ny,zs:ze,3), source=0.)
      allocate(self%twohnu3_c2(fs:fe))
      allocate(vproj(self%nx,self%ny,self%nz))
      self%b4_ptr => null()  ! departure coefficients may have changed since last execute
      self%twohnu3_c2 = twoh_c2 * self%nu**3
      call mem_register(self%chi_cont,  label='self%chi_cont')
      call mem_register(self%sigma,     label='self%sigma')
      call mem_register(self%beta,      label='self%beta')
      call mem_register(self%xhat,      label='self%xhat')
      call mem_register(self%ehat,      label='self%ehat')
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
        ! Jnu is zeroed for a fresh angle accumulation right after the
        ! background, so the foreground cannot read the previous iteration's
        ! mean intensity from it.  Keep (Jnu - Bnu) here instead, filled once
        ! per set where the continuum scattering term uses the same numbers.
        allocate(self%dJnu(fs:fe,self%nx,self%ny,self%nz), source=0.)
        call mem_register(self%dJnu, label='self%dJnu')
        ! one_m_eps is filled per foreground for the detailed line path; the
        ! sampled path bakes epsilon into the table instead
        allocate(self%one_m_eps(fs:fe,self%nx,self%ny,self%nz,2), source=0.)
        call mem_register(self%one_m_eps, label='self%one_m_eps')
      endif
      if (self%compute_psi) then
        allocate(self%Psi(fs:fe,self%nx,self%ny,self%nz), source=0.)
        call mem_register(self%Psi, label='Psi')
        ! PsiJ is a within-iteration angle accumulation for bound-free, but the
        ! scattering ALI needs it to survive to the next iteration's twin
        ! combine, like Jnu.  Allocate-if-absent here plus the per-set reset in
        ! background() covers both uses; deallocate_arrays frees it only when
        ! Jnu_to_file, where persistence goes through the file instead.
        if (self%bound_free .or. self%do_scat_precond) then
          if (.not. allocated(self%PsiJ)) then
            allocate(self%PsiJ(fs:fe,self%nx,self%ny,self%nz,self%n_set), source=0.)
            call mem_register(self%PsiJ, label='PsiJ')
          endif
        endif
        if (self%do_scat_precond .and. .not. allocated(self%Jali)) then
          allocate(self%Jali(fs:fe,self%nx,self%ny,self%nz,self%n_set), source=0.)
          call mem_register(self%Jali, label='Jali')
        endif
        ! Per-set seed flags persist for the patch's lifetime, so each set seeds
        ! the recurrence with its first formal solution exactly once.
        if (self%do_scat_precond .and. .not. allocated(self%Jali_valid)) then
          allocate(self%Jali_valid(self%n_set), source=.false.)
        endif
      endif
      if (get_cntrbf) then
        allocate(self%cntrbf_data(fs:fe,self%nx,self%ny,self%nz,2,self%mus:self%mue))
        call mem_register(self%cntrbf_data, 'self%cntrbf_data')
        ! (chi, eta, Inu) of the continuum set, per mu.  Its predecessor carried no
        ! mu axis, which was a latent bug: see save_contribution_function.
        allocate(self%cntrbf_cnt(fs:fe,self%nx,self%ny,self%nz,3,self%mus:self%mue))
        call mem_register(self%cntrbf_cnt, 'self%cntrbf_cnt')
      endif
    end subroutine allocate_arrays

    !-----------------------------------------------------------------------------
    subroutine deallocate_arrays()
      call mem_unregister(self%work     , label='self%work')
      call mem_unregister(self%chi_cont , label='self%chi_cont')
      call mem_unregister(self%sigma    , label='self%sigma')
      call mem_unregister(self%beta     , label='self%beta')
      call mem_unregister(self%xhat     , label='self%xhat')
      call mem_unregister(self%ehat     , label='self%ehat')
      deallocate(self%work)
      deallocate(self%chi_cont, self%sigma, self%beta, self%xhat, self%ehat)
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
        call mem_unregister(self%dJnu, label='self%dJnu')
        deallocate(self%chi_scat, self%one_m_eps, self%dJnu)
      endif
      if (allocated(self%b4)) deallocate(self%b4)
      if (allocated(self%Jnu) .and. Jnu_to_file) then
        call mem_unregister(self%Jnu, label='Jnu')
        deallocate(self%Jnu)
      endif
      if (self%compute_psi) then
        call mem_unregister(self%Psi, label='Psi')
        deallocate(self%Psi)
        ! PsiJ and Jali mirror the Jnu lifecycle: when the scattering ALI carries
        ! them to the next iteration's twin combine through memory they must stay
        ! alive across the execute; with Jnu_to_file persistence goes through the
        ! file instead, and without the ALI they are freed exactly as before.
        if (allocated(self%PsiJ) .and. .not. (self%do_scat_precond .and. .not. Jnu_to_file)) then
          call mem_unregister(self%PsiJ, label='PsiJ')
          deallocate(self%PsiJ)
        endif
      endif
      if (allocated(self%Jali) .and. Jnu_to_file) then
        call mem_unregister(self%Jali, label='Jali')
        deallocate(self%Jali)
      endif
      if (get_cntrbf) then
        call mem_unregister(self%cntrbf_data, label='self%cntrbf_data')
        deallocate(self%cntrbf_data)
        call mem_unregister(self%cntrbf_cnt, label='self%cntrbf_cnt')
        deallocate(self%cntrbf_cnt)
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

    ! i_set==1 here means "the first set iterated", not "the main set": the LTE
    ! background is set-independent and is computed once, then reused by the rest.
    if (i_set==1) then
      call self%lte_background_opacities(comp, include_absmet, do_sample_opac, do_make_eos)
    endif
    call self%lte_background_emission(comp, i_set, vec, include_absmet, do_rotate_atmos)
    call self%convert_to_work_units()
    call self%reset_array(self%Jnu(:,:,:,:,i_set:i_set))
    ! PsiJ now persists across iterations for the scattering ALI, so zero it for a
    ! fresh angle accumulation here -- after lte_background_emission has already
    ! read the previous iteration's twin PsiJ for the OAB correction.
    if (allocated(self%PsiJ)) self%PsiJ(:,:,:,:,i_set) = 0.
    call self%reset_array(self%work)
    
    do i=1, self%natom
      associate(atom=>self%atoms(i), trans=>self%trans(i))
        adata => atom%data
        if (do_rotate_atmos .and. i_mu /= 0) then
          adata => atom%rotated_data(i_mu)
        endif

        if (self%bound_free) then
          call self%add_bound_free_opacities(comp, atom, trans, adata, i_set)
        endif
        ! adamp is only allocated where the old-formalism loops still need it
        ! (atoms with hydrogen lines)
        ! i_set==1: computed once on the first set and reused (see above)
        if (self%bound_bound .and. i_set==1 .and. allocated(trans%adamp)) then
          call self%bound_bound_damping(comp, atom, trans, adata)
        endif
      end associate
    end do

    ! log-space tolerance reference for the SoA line kernels; the continuum xhat
    ! (work units, incl. bound-free) is fixed per i_set, so compute it once here and
    ! reuse it for all mu directions
    if (allocated(self%logref) .and. i_set /= self%i_set_cnt) then
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
              b = self%xhat(f,x,y,z)
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
    logical :: use_linelist, sample_ll, do_scat, use_store
    integer :: i, i_seg
    real(dp), allocatable :: window_lam(:)
    !-----------------------------------------------------------------------------
    ! detailed line list vs sampled LTE line table, per pass (rt_iter or rt_out)
    rtp => pass_params(self)
    ! An output patch on a window_file grid visits the lines of its whole window, however
    ! many patches that window is split into; unallocated, window_lam is absent.
    if (self%do_out_sp .and. self%i_nu < 0 .and. self%sp%n_segments > 0) then
      i_seg = output_patch_window(self%sp, self%i_nu, n_nu_out)
      window_lam = self%sp%lam([self%sp%segment_first(i_seg), self%sp%segment_last(i_seg)])
    endif
    use_linelist = comp%ll_initialised .and. rtp%include_linelist .and. .not. rtp%sample_linelist
    sample_ll = rtp%sample_linelist .and. associated(self%ll_opac)
    if (sample_ll) sample_ll = self%ll_opac%initialised

    ! line_scattering: accumulate the (1-eps)-weighted LTE line opacity, then
    ! redirect that fraction of the line source from Bnu to Jnu (see below)
    do_scat = line_scattering .and. (use_linelist .or. sample_ll)
    if (do_scat) self%chi_scat = 0.

    do i=1, self%natom
      associate(atom=>self%atoms(i), trans=>self%trans(i))
        adata => atom%data
        if (do_rotate_atmos .and. i_mu /= 0) then
          adata => atom%rotated_data(i_mu)
        endif
        if (self%bound_bound) then
          ! force LTE only for an LTE *companion* set: the main set of a run with
          ! no departures is already LTE and must not take this branch
          call self%add_bound_bound_opacities(comp, atom, trans, adata, wmu, vproj, &
                                              LTE=self%set_kind(i_set)==set_lte .and. &
                                                  self%has_nlte_set, &
                                              use_ll=self%use_ll(i), i_set=i_set)
        endif
      end associate
    end do

    ! NLTE departure coefficients for this (i_set, i_mu), shared by both the full
    ! detailed line list and the sampled-mode detailed departure-species pass
    departures => null()
    if (associated(self%m3d_parent%departures) .and. &
        self%set_kind(i_set)==set_nlte) then
      i = merge(i_mu+1, 1, do_rotate_atmos)
      departures => self%m3d_parent%departures(:,i)
    endif

    ! The stored line opacity of a 1D output patch holds for any direction as long
    ! as the line profiles are unshifted; a moving atmosphere recomputes every call.
    use_store = use_linelist .and. self%do_out_sp .and. &
                associated(self%m3d_parent%ll_point_tables)
    if (use_store) use_store = all(vproj == 0.)

    if (use_store) then
      call stored_linelist()
    elseif (use_linelist) then
      call detailed_linelist()
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
            departures=departures, only_departures=.true., i_set=i_set, &
            window_lam=window_lam)
        end associate
      endif
    endif

    ! line_scattering: the LTE lines were added as pure thermal absorbers (source
    ! Bnu).  The scattering fraction (1-eps) of their opacity now has its source
    ! redirected from Bnu to the mean intensity, Snu_line = eps*Bnu + (1-eps)*Jnu,
    ! by adding Delta_eta = chi_scat*(Jnu - Bnu) to the emissivity without
    ! altering chi.  That is done in work_to_transfer_units, which already forms
    ! eta and so absorbs the whole correction into two extra flops per point --
    ! no separate pass over the work arrays, and no division by beta.
    self%do_line_scat = do_scat

  contains

    !-----------------------------------------------------------------------------
    ! Full detailed line list (non-sampled mode): every species, with departures.
    ! chi_lines/eta_lines, when present, receive the increments added to the work
    ! arrays.
    subroutine detailed_linelist(chi_lines, eta_lines)
      real(4), dimension(:,:,:,:), optional :: chi_lines, eta_lines
      !---------------------------------------------------------------------------
      associate(xhat=>self%work(:,:,:,1:self%nz,1), ehat=>self%work(:,:,:,1:self%nz,2))
        if (do_scat) then
          ! per-stage scattering fraction on the detailed-path frequency grid
          call compute_epsilon_stage(1, self%sub_sp%nu, comp%temp, comp%ne, &
            self%one_m_eps(:,:,:,:,1))
          call compute_epsilon_stage(2, self%sub_sp%nu, comp%temp, comp%ne, &
            self%one_m_eps(:,:,:,:,2))
          call comp%add_linelist_opacities2(self%sub_sp, xhat, ehat, dlam, vproj, &
            departures=departures, one_m_eps=self%one_m_eps, chi_scat=self%chi_scat, &
            beta=self%beta, i_set=i_set, window_lam=window_lam, &
            chi_lines=chi_lines, eta_lines=eta_lines)
        else
          call comp%add_linelist_opacities2(self%sub_sp, xhat, ehat, dlam, vproj, &
                                            departures=departures, i_set=i_set, &
                                            window_lam=window_lam, &
                                            chi_lines=chi_lines, eta_lines=eta_lines)
        endif
      end associate
    end subroutine detailed_linelist

    !-----------------------------------------------------------------------------
    ! The detailed line list of a 1D output patch, computed once per set and then
    ! copied.  Temperature, densities and departures stay fixed through the output
    ! scattering iterations and the final pass.  The one input that moves is the
    ! kernels' wing-cutoff reference: work xhat carries sigma*(Jnu-Bnu), so later
    ! passes would cut line wings against a slightly different continuum.  The store
    ! keeps the cutoff of the first pass, where Jnu=Bnu; the difference is below
    ! the kernel tolerance.  The first foreground of the (i_nu, set) fills the store
    ! while its i_mu twins and the final patch wait on the lock.  Adding the stored
    ! increments repeats the kernel's own last step.
    subroutine stored_linelist()
      type(point_opac_set_t), pointer :: store
      type(lock_t), pointer :: lock
      logical :: filled_here
      integer :: x, y, z, f
      !---------------------------------------------------------------------------
      if (.not. one_dimensional) &
        error stop 'rt_patch%foreground: the stored line-list opacity requires a 1D atmosphere'
      lock => self%m3d_parent%ll_point_locks(self%i_nu)
      filled_here = .false.
      call lock%set('ll_point_table')
      associate(table => self%m3d_parent%ll_point_tables(self%i_nu))
        if (.not. allocated(table%set)) allocate(table%set(self%n_set))
        store => table%set(i_set)
      end associate
      if (.not. store%initialised) then
        allocate(store%chi(self%fs:self%fe,self%nx,self%ny,self%nz), &
                 store%eta(self%fs:self%fe,self%nx,self%ny,self%nz))
        call mem_register(store%chi, label='ll_point_table')
        call mem_register(store%eta, label='ll_point_table')
        call detailed_linelist(store%chi, store%eta)
        if (do_scat) then
          allocate(store%scat, source=self%chi_scat)
          call mem_register(store%scat, label='ll_point_table')
        endif
        store%initialised = .true.
        filled_here = .true.
      endif
      call lock%unset('ll_point_table')
      if (filled_here) return

      do z=1, self%nz
        do y=1, self%ny
          do x=1, self%nx
            do f=self%fs, self%fe
              self%work(f,x,y,z,2) = self%work(f,x,y,z,2) + store%eta(f,x,y,z)
              self%work(f,x,y,z,1) = self%work(f,x,y,z,1) + store%chi(f,x,y,z)
            end do
          end do
        end do
      end do
      if (do_scat) self%chi_scat = store%scat
    end subroutine stored_linelist

  END SUBROUTINE foreground

  !===============================================================================
  !> Free the stored line-list opacity of this patch's output i_nu.  The final pass
  !> is its last reader: the i_mu twins have finished their iterations before it runs.
  !===============================================================================
  SUBROUTINE free_ll_point_table (self)
    class(rt_patch_t) :: self
    integer :: i
    !-----------------------------------------------------------------------------
    if (.not. associated(self%m3d_parent%ll_point_tables)) return
    call self%m3d_parent%ll_point_locks(self%i_nu)%set('ll_point_table')
    associate(table => self%m3d_parent%ll_point_tables(self%i_nu))
      if (allocated(table%set)) then
        do i=1, size(table%set)
          associate(store => table%set(i))
            if (allocated(store%chi)) then
              call mem_unregister(store%chi, label='ll_point_table')
              call mem_unregister(store%eta, label='ll_point_table')
            endif
            if (allocated(store%scat)) call mem_unregister(store%scat, label='ll_point_table')
          end associate
        end do
        deallocate(table%set)
      endif
    end associate
    call self%m3d_parent%ll_point_locks(self%i_nu)%unset('ll_point_table')
  END SUBROUTINE free_ll_point_table

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
          case('hermite')
            call self%short_char_hermite(comp, vec)
          case default
            stop 'Unknown rt_iteration_params solver, use "linear" or "hermite"'
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
    get_PsiJ = (any(self%get_rates) .and. self%bound_free) .or. self%do_scat_precond

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

    ! Twin sub-patches only exist for do_full_sp when rt_iter%maxiter>0, and for
    ! the output-spectrum scattering loop when rt_out%maxiter>0 (see
    ! composition_patch_mod.f90 create_sub_patches calls) -- saving Jnu outside
    ! those conditions leaves it with no reader.
    !
    ! The i_mu=[1,2] hemisphere pair created for rt_out%maxiter>0 is do_out_sp
    ! and *not* do_final, so the second test has to be do_out_sp: those are the
    ! patches that both read (each other, at the next iteration) and are read by
    ! the do_final ray.  Testing do_final here was dead code -- iso_integrations
    ! is only reached from the non-final branch of execute() -- and left the
    ! Jnu_out files unwritten, so the readers below died on a missing file.
    if (Jnu_to_file .and. ((self%do_full_sp .and. rt_iter%maxiter > 0) .or. &
                           (self%do_out_sp .and. rt_out%maxiter > 0))) then
      Jnu_label = merge('Jnu_out', 'Jnu    ', self%do_out_sp)
      call self%save5(self%Jnu, trim(Jnu_label), key=self%key(1:3), add2meta=save_Jnu)
    endif
    ! Scattering ALI: persist the hemisphere PsiJ (twin-combined at the next
    ! iteration's OAB step) and the accelerated Jali, mirroring the Jnu channel.
    ! iso_integrations runs once per set and saves the whole array, so the last
    ! call of an execute holds the complete state; the reader only ever sees the
    ! writer's last save.  Seed flags are raised per set in the OAB block.
    !
    ! add2meta=.false.: these are scratch, deleted again in finish(), and the
    ! readers below pass mem= (never alloc=), so nothing consults the meta file.
    ! Indexing them would append a line per patch per set per iteration, under
    ! the global io lock, describing files that no longer exist by run's end.
    if (self%do_scat_precond .and. Jnu_to_file) then
      call self%save5(self%PsiJ, 'PsiJ', key=self%key(1:3), add2meta=.false.)
      call self%save5(self%Jali, 'Jali', key=self%key(1:3), add2meta=.false.)
    endif
  END SUBROUTINE iso_integrations

  !===============================================================================
  SUBROUTINE save_debug_arrays(self, i_set, i_mu)
    class(rt_patch_t) :: self
    integer, intent(in) :: i_set, i_mu
    character(6) :: sav_label
    !-----------------------------------------------------------------------------
    if (save_chi) then
      if (i_set==self%i_set_main) then
        write(sav_label, '(A, I2.2)') 'chi_', i_mu
        call self%save4(self%work(:,:,:,1:self%nz,1), sav_label, key=self%key(1:3))
      endif
    endif
    if (save_Snu) then
      if (i_set==self%i_set_main) then
        write(sav_label, '(A, I2.2)') 'Snu_', i_mu
        call self%save4(self%work(:,:,:,1:self%nz,2), sav_label, key=self%key(1:3))
      endif
    endif
    if (save_Inu) then
      if (i_set==self%i_set_main) then
        write(sav_label, '(A, I2.2)') 'Inu_', i_mu
        call self%save4(self%work(:,:,:,1:self%nz,3), sav_label, key=self%key(1:3))
      elseif (i_set==self%i_set_cnt) then
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
            self%beta(f,x,y,z) = max(exp(-hv_kT), tiny(1.0))
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
            ! shift into the foreground rather than the velocity-blind background.
            call self%generate_opacity_table(tab_sp, include_absmet, include_llist=.false., &
              build_deriv=.true., disjoint=.true.)
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
        self%chi_cont, self%sigma, self%sp%lam, self%fs, self%fe)
    else
      call comp%opacity(self%sub_sp, self%chi_cont, self%sigma, &
        include_absmet=include_absmet, include_llist=.false., &
        beta=self%beta, sorted=self%sp%sorted)
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
    real(4), dimension(:,:,:,:,:), allocatable, target :: PsiJup_mem, PsiJdn_mem, Jali_mem
    real(4), dimension(:,:,:,:,:), pointer, contiguous :: PsiJup, PsiJdn
    real(4), dimension(:,:,:,:), allocatable :: Jali_old
    real(8) :: op, Jfs, Jali_new, op_raw, op_max, amp_max
    integer :: n_neg
    real(4) :: Bnu, wscat
    character(8) :: Jnu_label
    logical :: Jnu_to_Bnu, do_scat_now
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
      !-------------------------------------------------------------------------
      ! Continuum-scattering ALI (Olson-Auer-Buchler local operator).
      !
      ! self%Jnu now holds the full-sphere formal-solution mean intensity
      ! J_FS = 0.5*(Jup+Jdn).  The source it feeds is
      !   eta = chi_tot*Bnu + sigma_scat*(J - Bnu),
      ! so J appears on both sides and plain Lambda-iteration converges only at
      ! the scattering albedo per iteration.  Replace J_FS with
      !   Jali = (J_FS - op*Jali_old)/(1 - op),   op = sigma_scat*PsiJ_full,
      ! the standard (1 - Lambda* a)^-1 acceleration.  Psi is stored as
      ! Lambda_diag/chi (rt_solvers_mod), so the dimensionless self-coupling is
      ! sigma_scat*PsiJ -- multiply by the scattering *opacity*, mirroring the
      ! bound-free `phi = PsiJ * cont%alpha`.  PsiJ_full = 0.5*(PsiJup+PsiJdn)
      ! because each twin integrates only its own hemisphere; unlike bound-free,
      ! which consumes PsiJ linearly, this update divides by (1-op) and so must
      ! see the full sphere.
      !
      ! op <= 1 identically: sigma is the scattering part of chi_cont
      ! by construction on both the exact and the sampled path, and
      ! Lambda_diag <= 1.  So scat_op_max is not there to rescue an out-of-range
      ! ratio; it bounds the amplification 1/(1-op) where op legitimately
      ! approaches 1 (see its declaration in m3d_patch_mod).  The difference
      ! in the numerator cancels heavily as op -> 1, so it is formed in real(8).
      ! Where there is no scattering op = 0 and Jali = J_FS exactly.
      !-------------------------------------------------------------------------
      ! Both twins must be carrying a PsiJ.  do_scat_precond follows
      ! any(get_rates), which drops as soon as the atom converges, and the twins
      ! reach that point on their own schedules -- so a twin can have freed its
      ! PsiJ while this one still wants it.  Skip the correction for this
      ! iteration rather than dereferencing a dead array.  This is a local
      ! decision only: self%do_scat_precond stays up so PsiJ keeps accumulating
      ! and the recurrence can resume as soon as the twin is back.
      do_scat_now = self%do_scat_precond
      if (do_scat_now .and. .not. Jnu_to_file) then
        do_scat_now = allocated(self%up_twin%PsiJ) .and. allocated(self%dn_twin%PsiJ)
      endif
      if (do_scat_now) then
        allocate(Jali_old(fs:fe,self%nx,self%ny,self%nz))
        if (Jnu_to_file) then
          allocate(PsiJup_mem(fs:fe,self%nx,self%ny,self%nz,self%up_twin%n_set))
          allocate(PsiJdn_mem(fs:fe,self%nx,self%ny,self%nz,self%dn_twin%n_set))
          call self%up_twin%load5('PsiJ', key=self%up_twin%key(1:3), mem=PsiJup_mem)
          call self%dn_twin%load5('PsiJ', key=self%dn_twin%key(1:3), mem=PsiJdn_mem)
          PsiJup => PsiJup_mem
          PsiJdn => PsiJdn_mem
        else
          PsiJup => self%up_twin%PsiJ
          PsiJdn => self%dn_twin%PsiJ
        endif
        if (self%Jali_valid(i_set)) then
          if (Jnu_to_file) then
            allocate(Jali_mem(fs:fe,self%nx,self%ny,self%nz,self%n_set))
            call self%load5('Jali', key=self%key(1:3), mem=Jali_mem)
            Jali_old = Jali_mem(:,:,:,:,i_set)
          else
            Jali_old = self%Jali(:,:,:,:,i_set)
          endif
        else
          Jali_old = self%Jnu(:,:,:,:,i_set)   ! seed the recurrence with J_FS
        endif
        op_raw = 0
        op_max = 0
        amp_max = 0
        n_neg = 0
        do z=1, self%nz
          do y=1, self%ny
            do x=1, self%nx
              do f=fs, fe
                op = self%sigma(f,x,y,z) * 0.5d0 * &
                     (dble(PsiJup(f,x,y,z,i_set)) + dble(PsiJdn(f,x,y,z,i_set)))
                op_raw = max(op_raw, op)
                op = max(0.0d0, min(op, scat_op_max))
                op_max = max(op_max, op)
                Jfs = self%Jnu(f,x,y,z,i_set)
                Jali_new = (Jfs - op*dble(Jali_old(f,x,y,z))) / (1.0d0 - op)
                if (Jali_new < 0d0) then
                  n_neg = n_neg + 1
                  ! A negative mean intensity is unphysical and feeds a negative
                  ! source function straight into the next formal solution.  It
                  ! means the local operator over-corrected; fall back to the
                  ! unaccelerated value for this point rather than propagate it.
                  Jali_new = Jfs
                endif
                if (Jfs > 0d0) amp_max = max(amp_max, Jali_new/Jfs)
                self%Jnu(f,x,y,z,i_set) = real(Jali_new, 4)
              end do
            end do
          end do
        end do
        if (diag_pops) then
          !$omp critical (diag_scat_print)
          print '(A,I8,A,I3,A,F9.5,A,F9.5,A,ES11.3,A,I9,A,I9)', &
            'DIAGS ', self%id, ' set=', i_set, ' op_raw_max=', op_raw, &
            ' op_max=', op_max, ' amp_max=', amp_max, &
            ' Jneg=', n_neg, '/', (fe-fs+1)*self%nx*self%ny*self%nz
          !$omp end critical (diag_scat_print)
        endif
        self%Jali(:,:,:,:,i_set) = self%Jnu(:,:,:,:,i_set)
        self%Jali_valid(i_set) = .true.
        if (allocated(PsiJup_mem)) deallocate(PsiJup_mem, PsiJdn_mem)
        if (allocated(Jali_mem)) deallocate(Jali_mem)
        deallocate(Jali_old)
      endif
    elseif (.not. Jnu_to_Bnu) then
      call self%reset_array(self%Jnu(:,:,:,:,i_set:i_set))
    endif

    !---------------------------------------------------------------------------
    ! The LTE background source, including coherent scattering.
    !
    ! sigma is the scattering part *of* chi_cont, so
    !   eta = chi_true*Bnu + sigma*Jnu = chi_cont*Bnu + sigma*(Jnu - Bnu),
    ! and the second form is the one to use: every optical depth stays exactly
    ! as it is and only the departure of J from Planck drives the scattering
    ! term.  Because sigma is built together with chi_cont -- summed after
    ! interpolation on the sampled path, accumulated side by side on the exact
    ! one -- the effective single-scattering albedo sigma/chi_cont is <= 1 by
    ! construction, which is what makes this plain Lambda-iteration a
    ! contraction rather than an amplifier.
    !
    ! Jnu is left holding (Jnu - Bnu) unconditionally.  The foreground's
    ! line-scattering term needs the same quantity, but Jnu is zeroed for a
    ! fresh angle accumulation before it runs, so keep a copy in dJnu.
    !
    ! diag_no_scat_loop severs the scattering Lambda-loop while holding every
    ! optical depth fixed, by zeroing the weight rather than branching inside
    ! the loop -- see m3d_debug_mod.
    !---------------------------------------------------------------------------
    wscat = 1.0
    if (diag_no_scat_loop) wscat = 0.0
    do z=1, self%nz
      do y=1, self%ny
        do x=1, self%nx
          do f=fs, fe
            associate(Jnu=>self%Jnu(f,x,y,z,i_set), chi=>self%chi_cont(f,x,y,z), &
                      eta=>self%ehat(f,x,y,z), beta=>self%beta(f,x,y,z), &
                      sigma=>self%sigma(f,x,y,z))
              Bnu = beta / (1 - beta)
              if (Jnu_to_Bnu) then
                Jnu = Bnu
              endif
              Jnu = Jnu - Bnu
              eta = chi*Bnu + wscat*sigma*Jnu
            end associate
            if (allocated(self%dJnu)) self%dJnu(f,x,y,z) = wscat*self%Jnu(f,x,y,z,i_set)
          end do
        end do
      end do
    end do
    call trace%end(itimer)
  END SUBROUTINE

  !===============================================================================
  SUBROUTINE generate_opacity_table (self, sp, include_absmet, include_llist, build_deriv, &
                                     disjoint)
    class(rt_patch_t), target :: self
    class(composition_type), pointer :: comp
    type(spectrum_type) :: sp
    real(4), dimension(:,:,:,:), pointer :: all_chi, all_scat
    logical, intent(in) :: include_absmet, include_llist, build_deriv
    ! Both components are always filled.  disjoint (the sampling path) keeps
    ! val = true absorption and scat = scattering, each normalised by chi500, so
    ! that sample_interpolate can sum them after interpolation and the sampled
    ! albedo stays <= 1.  Without it (the make_eos path) val is the total and
    ! scat stays in raw 1/cm, which is what save_opac_table writes out.
    logical, intent(in), optional :: disjoint
    logical :: do_disjoint
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
    do_disjoint = .false.
    if (present(disjoint)) do_disjoint = disjoint
    allocate(self%opac%scat(nnu, ntemp, nne), source=0.)
    call mem_register(self%opac%scat, label='opac_table')
    all_scat(1:nnu, 1:ntemp, 1:nne, 1:1) => self%opac%scat
    call comp%opacity(sp, all_chi, all_scat, include_absmet=include_absmet, &
      include_llist=include_llist, dlam=dlam, sorted=self%sp%sorted, &
      disjoint=do_disjoint)

    do k=1, nne
      do j=1, ntemp
        over_chi500 = 1 / comp%chi500(j,k,1)
        self%opac%val(:,j,k) = self%opac%val(:,j,k) * over_chi500
        ! the sampled table carries both components in the same units; the
        ! make_eos table leaves scat in 1/cm for save_opac_table
        if (do_disjoint) self%opac%scat(:,j,k) = self%opac%scat(:,j,k) * over_chi500
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

    ! Continuum reference for add_linelist's opacity threshold.  The pipeline
    ! hands back the working form directly -- continuum_opacities returns
    ! xhat = ehat = chi/(1-beta) -- so there is nothing to convert here; only the
    ! scattering opacity is dropped, since a threshold reference wants the
    ! opacity the lines are compared against, not the source function.
    all_chi(1:nnu, 1:ntemp, 1:nne, 1:1) => self%ll_opac%val
    allocate(all_eta, all_beta, xhat_cont, mold=all_chi)
    call comp%boltzmann_factor(psp, all_beta)
    call comp%continuum_opacities(psp, all_chi, beta=all_beta, sorted=self%sp%sorted)
    if (include_absmet) &
      call comp%add_absmet_opacities(psp, all_chi, all_beta)
    all_eta = all_chi       ! in LTE the continuum's two working forms coincide
    xhat_cont = all_chi     ! snapshot of the continuum working xhat

    ! LTE line list (no departures/vproj); tab_comp is masked to non-departure
    ! species, so only LTE species contribute here
    if (line_scattering) then
      ! scattering split: precompute the per-stage scattering fraction 1-eps at
      ! the node (T,ne) on the (padded) table grid, accumulate the
      ! (1-eps)-weighted line opacity into chi_scat, and store it as the table's
      ! second component (physical 1/cm, normalised the same way as the total)
      allocate(one_m_eps(nnu, ntemp, nne, 1, 2))
      allocate(chi_scat(nnu, ntemp, nne, 1), source=0.)
      call compute_epsilon_stage(1, psp%nu, comp%temp, comp%ne, one_m_eps(:,:,:,:,1))
      call compute_epsilon_stage(2, psp%nu, comp%temp, comp%ne, one_m_eps(:,:,:,:,2))
      call comp%add_linelist_opacities2(psp, all_chi, all_eta, dlam, &
        one_m_eps=one_m_eps, chi_scat=chi_scat, beta=all_beta)
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
    real(4) :: over_rho
    integer :: i,j,k, ntemp, nrho
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%save_opac_table', itimer=itimer)
    comp => self%m3d_parent%tab_comp
    ntemp = size(self%opac%val, 2)
    nrho = size(self%opac%val, 3)

    ! save absorption part first (which has been normalised by chi500)
    do k=1, nrho
      do j=1, ntemp
        over_rho = comp%chi500(j,k,1) / comp%rho(j,k,1)
        self%opac%val(:,j,k) = self%opac%val(:,j,k) * over_rho
      end do
    end do
    call self%save3(self%opac%val, 'opa')

    ! then save the scattering part, which generate_opacity_table left in 1/cm
    do k=1, nrho
      do j=1, ntemp
        over_rho = 1 / comp%rho(j,k,1)
        self%opac%scat(:,j,k) = self%opac%scat(:,j,k) * over_rho
      end do
    end do
    call self%save3(self%opac%scat, 'scat')
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
  SUBROUTINE add_bound_free_opacities (self, comp, atom, trans, adata, i_set)
    class(rt_patch_t), target :: self
    class(composition_type) :: comp
    class(atom_type) :: atom
    class(atom_trans_t) :: trans
    class(atom_data_t), target :: adata
    integer, intent(in) :: i_set
    class(cont_type), pointer :: cont
    real(4) :: ni, nj, tmp, alpha
    integer :: i, fs, fe
    integer, save :: itimer=0
    integer :: x, y, z, f
    real(8), allocatable :: pop_ratio(:,:,:)
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%add_bound_free_opacities', itimer=itimer)
    allocate(pop_ratio, mold=adata%totn)
    call comp%population_ratio(atom%id, i_set, pop_ratio)
    
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
              ni  = adata%nstar(x,y,z,cont%i) * adata%b(x,y,z,cont%i) * pop_ratio(x,y,z)
              tmp = adata%nstar(x,y,z,cont%i) * adata%b(x,y,z,cont%j) * pop_ratio(x,y,z)
              do f=fs, fe
                alpha = cont%alpha(f)
                self%xhat(f,x,y,z) = self%xhat(f,x,y,z) + alpha * ni
                self%ehat(f,x,y,z) = self%ehat(f,x,y,z) + alpha * tmp
              end do
            end do
          end do
        end do
      else ! b = 1
        do z=1, self%nz
          do y=1, self%ny
            do x=1, self%nx
              ni  = adata%nstar(x,y,z,cont%i) * pop_ratio(x,y,z)
              tmp = adata%nstar(x,y,z,cont%i) * pop_ratio(x,y,z)
              do f=fs, fe
                alpha = cont%alpha(f)
                self%xhat(f,x,y,z) = self%xhat(f,x,y,z) + alpha * ni
                self%ehat(f,x,y,z) = self%ehat(f,x,y,z) + alpha * tmp
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
  SUBROUTINE add_bound_bound_opacities (self, comp, atom, trans, adata, wmu, vproj, LTE, use_ll, &
                                        i_set)
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
    integer, intent(in) :: i_set
    real(8), allocatable :: pop_ratio(:,:,:)
    integer, save :: itimer=0, itimer1=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin (trim(self%kind)//'%add_bound_bound_opacities', itimer=itimer)
    ! Copying continuum opacities to angle dependent opacities

    allocate(vphi(self%fs:self%fe))
    allocate(wphi, mold=vproj)
    get_wphi = self%iter == 0 .and. .not. self%do_out_sp
    do_LTE = LTE .or. .not. allocated(adata%b)
    allocate(pop_ratio, mold=adata%totn)
    call comp%population_ratio(atom%id, i_set, pop_ratio)

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
            ni = adata%nstar(x,y,z,line%i) * pop_ratio(x,y,z)
            nj = adata%nstar(x,y,z,line%j) * pop_ratio(x,y,z)

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
                  phi = phi + qlg_bot_profile(a, line%q(fs:fe), vpr, idf) * line%wt(ih)
                end do
              else
                vpr = vproj(x,y,z)
                a = trans%adamp(x,y,z,i)
                phi = qlg_bot_profile(a, line%q(fs:fe), vpr, idf)
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
      call get_ll_env(comp, atom, trans, adata, e_kT4, ne4, tothi4, sp4, idf3, vb3, &
                      lnbase, i_set)
      associate(xhat=>self%work(:,:,:,1:self%nz,1), ehat=>self%work(:,:,:,1:self%nz,2), &
                cont=>self%xhat)
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
                         idf3, vb3, lnbase, i_set)
    class(composition_type), target :: comp
    class(atom_type) :: atom
    class(atom_trans_t) :: trans
    class(atom_data_t) :: adata
    real(4), dimension(:,:,:),   allocatable, intent(out) :: e_kT, ne, tothi, stark_pe
    real(4), dimension(:,:,:),   allocatable, intent(out) :: idf3, vb3
    real(4), dimension(:,:,:,:), allocatable, intent(out) :: lnbase
    integer, intent(in), optional :: i_set
    class(element_type), pointer :: element
    class(molec_type), pointer :: mol
    integer :: m(3), np, s
    integer :: set_index
    real(8), dimension(:,:,:), pointer :: set_n
    real(8) :: set_scale
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%get_ll_env', itimer=itimer)
    m = shape(adata%inv_df)
    set_index = 0
    if (present(i_set)) set_index = i_set
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
          call comp%set_population(ion%id, set_index, set_n, set_scale)
          if (associated(set_n)) then
            lnbase(:,:,:,s) = log(max(set_n,tiny(1d0)))
          else
            lnbase(:,:,:,s) = log(max(ion%n,tiny(1d0))) + real(log(set_scale),4)
          endif
          lnbase(:,:,:,s) = lnbase(:,:,:,s) - ion%part &
                          + real(comp%dxi) * (s-1) * e_kT + log(idf3)
        end associate
      end do
    else
      call comp%molecule_by_id(atom%id, mol, verbose=.false.)
      if (associated(mol)) then
        if (allocated(mol%n)) then
          call comp%set_population(mol%id, set_index, set_n, set_scale)
          if (associated(set_n)) then
            lnbase(:,:,:,1) = log(max(set_n,tiny(1d0))) - mol%part + log(idf3)
          else
            lnbase(:,:,:,1) = log(max(mol%n,tiny(1d0))) - mol%part + log(idf3) &
                            + real(log(set_scale),4)
          endif
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

    ! use_ll routes the non-hydrogen lines to the SoA kernel, which covers opacity
    ! only: there is no rate counterpart, so those lines would end up with opacity
    ! and no radiative rates.  pre_init holds use_ll at .false. for that reason;
    ! stop here rather than let a silently wrong rate matrix through.
    if (use_ll) error stop 'bound_bound_rates: use_ll has no rate kernel (see pre_init)'

    do i=1, size(trans%line)
      line => trans%line(i)
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
                    vgt = qlg_bot_profile(a2, q, vpr, idf)
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
                  phi = qlg_bot_profile(a, line%q(f), vpr, idf)
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
              associate(Jnu=>self%Jnu(f,x,y,z,self%i_set_main), beta=>self%beta(f,x,y,z), &
                        Psi=>self%PsiJ(f,x,y,z,self%i_set_main))
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
              self%PsiJ(f,x,y,z,i_set) = self%PsiJ(f,x,y,z,i_set) + wmu2 * self%Psi(f,x,y,z)
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
      if (i_mu==self%mus .and. i_set==self%i_set_main) then   ! sample the spectrum set once
        memcheck(1:3) = self%work(self%fcheck, p(1), p(2), p(3),:)
        !memcheck(2:3) = memcheck(2:3) * self%twohnu3_c2(self%fcheck)
        associate(beta=>self%beta(self%fcheck, p(1), p(2), p(3)))
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
  ! ehat enters holding the physical emissivity eta (from lte_background_emission)
  ! and leaves holding the working form eta/beta.  xhat must therefore be formed
  ! first, while ehat is still eta -- the two statements do not commute.
  SUBROUTINE convert_to_work_units(self)
    class(rt_patch_t) :: self
    integer :: fs, fe
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    fs = self%fs
    fe = self%fe
    do z=1, self%nz
      do y=1, self%ny
        do x=1, self%nx
          do f=fs, fe
            self%xhat(f,x,y,z) = self%chi_cont(f,x,y,z) + self%ehat(f,x,y,z)
            ! beta is clamped nonzero (see lte_background_opacities)
            self%ehat(f,x,y,z) = self%ehat(f,x,y,z) / self%beta(f,x,y,z)
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
    ! seed the angle-dependent working arrays from the direction-independent background
    do z=1, self%nz
      do y=1, self%ny
        do x=1, self%nx
          do f=self%fs, self%fe
            self%work(f,x,y,z,1) = self%xhat(f,x,y,z)
            self%work(f,x,y,z,2) = self%ehat(f,x,y,z)
          end do
        end do
      end do
    end do
    call trace%end(itimer)
  END SUBROUTINE back_to_work

  !===============================================================================
  !===============================================================================
  ! Close the accumulation: turn the working pair (xhat, ehat) into the (chi, Snu)
  ! pair the formal solvers consume, and apply the line-scattering source on the
  ! way.
  !
  !   chi = xhat - ehat*beta                    physical extinction, scattering
  !                                             opacity included (the lines went
  !                                             in as absorbers, so chi is whole)
  !   eta = ehat*beta + chi_scat*(Jnu - Bnu)    redirect the (1-eps) fraction of
  !                                             the line source from Bnu to Jnu
  !   Snu = eta/chi
  !===============================================================================
  SUBROUTINE work_to_transfer_units (self)
    class(rt_patch_t) :: self
    logical :: do_scat
    integer, save :: itimer=0
    integer :: x, y, z, f, i
    real(4) :: xhat, ehat, beta, eta, chi, floor
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%work_to_transfer_units', itimer=itimer)
    do_scat = self%do_line_scat
    do z=1, self%nz
      do y=1, self%ny
        do x=1, self%nx
          do f=self%fs, self%fe
            xhat = self%work(f,x,y,z,1)
            ehat = self%work(f,x,y,z,2)
            beta = self%beta(f,x,y,z)
            eta  = ehat * beta
            chi  = xhat - eta
            if (do_scat) eta = eta + self%chi_scat(f,x,y,z)*self%dJnu(f,x,y,z)
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
                    self%beta(f,x,y,z)
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
  END SUBROUTINE work_to_transfer_units

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
            associate(beta=>self%beta(f,x,y,z), Inu=>self%work(f,x,y,z,3))
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
    real(4) :: net_chi, net_src, net_eta, lin_chi, lin_eta
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    ! The continuum set now leads, so it is the one that stashes and the main set
    ! that combines -- the reverse of the historical order.  The stash keeps all
    ! three continuum quantities (chi, eta, Inu) rather than pre-combining them,
    ! so the main set still forms the differences chi-chi_cnt and eta-eta_cnt
    ! directly.  That matters: deep in the atmosphere those are catastrophically
    ! cancelling, and any algebraically-equivalent regrouping that multiplies out
    ! first loses most of the significant digits there.
    !
    ! The stash also carries a mu axis, which its predecessor (cntrbf_tmp) did
    ! not.  In the set-outer/mu-inner loop of `execute` that meant every mu was
    ! paired with only the *last* mu of the other set -- harmless for a static
    ! atmosphere, where chi and eta do not depend on direction, but wrong as soon
    ! as vproj does (any snapshot with velocities, i.e. the default rotate_atmos
    ! = .false. path).  Fixed here rather than reproduced.
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch_t%save_contribution_function', itimer=itimer)
    do z=1,self%nz
      do y=1,self%ny
        do x=1,self%nx
          do f=self%fs,self%fe
            associate(cnt => self%cntrbf_cnt(f,x,y,z,:,i_mu), &
                      cf_out => self%cntrbf_data(f,x,y,z,:,i_mu))
              if (i_set == self%i_set_cnt) then ! continuum set: stash
                cnt(1) = self%work(f,x,y,z,1)                ! chi_cnt
                cnt(2) = self%work(f,x,y,z,2)*cnt(1)         ! eta_cnt
                cnt(3) = self%work(f,x,y,z,3)                ! I_cnt
              elseif (i_set == self%i_set_main) then ! main set: combine
                net_chi = self%work(f,x,y,z,1)
                net_src = self%work(f,x,y,z,2)

                ! net_eta is a rounding point, not a convenience: folding the
                ! product into the subtraction lets the compiler contract it into
                ! an FMA, and the difference below is catastrophically cancelling
                ! deep in the atmosphere, where a single unrounded product moves
                ! the result by O(1) in relative terms
                net_eta = net_src*net_chi
                lin_chi = net_chi - cnt(1)
                lin_eta = net_eta - cnt(2)

                cf_out(1) = merge(cnt(3)*lin_chi - lin_eta, net_chi*net_src, &
                                  depression_cntrbf)
                cf_out(1) = cf_out(1) * self%twohnu3_c2(f)
                cf_out(2) = net_chi
              endif
            end associate
          end do
        end do
      end do
    end do
    if (i_set==self%i_set_main .and. i_mu==self%mue) then
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
          !-----------------------------------------------------------------
          ! The scattering ALI persists PsiJ and Jali alongside Jnu (see
          ! iso_integrations), but unlike Jnu they are pure iteration scratch
          ! -- there is no save_PsiJ -- so they go regardless of save_Jnu.
          ! The twins are exactly the writers: for do_full_sp this patch and
          ! its partner, for do_final the i_mu=[1,2] output pair.
          !
          ! Gate on the namelist flag rather than on the twins'
          ! do_scat_precond: that one is recomputed every execute and drops
          ! with get_rates as atoms converge, so a twin that wrote files
          ! earlier in the iteration may no longer admit to it here.  remove()
          ! is a no-op when the file is absent, so the loose gate is free.
          !-----------------------------------------------------------------
          if (scattering_preconditioning) then
            call self%up_twin%remove('PsiJ', key=self%up_twin%key(1:3))
            call self%dn_twin%remove('PsiJ', key=self%dn_twin%key(1:3))
            call self%up_twin%remove('Jali', key=self%up_twin%key(1:3))
            call self%dn_twin%remove('Jali', key=self%dn_twin%key(1:3))
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
      m = [self%nf, self%n(1), self%n(2), self%n_mu*self%n_set, 1, 3]
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
        sum(sum(dble(self%mem(fs:fe,:,:,i,self%i_set_main,3)), dim=3), dim=2) &
          / dble(self%nx * self%ny)
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

  !===============================================================================
  !> Verify that the n_nu_out output patches of a segmented (window_file) grid tile
  !> it exactly: in i_nu order each patch is non-empty, starts where the previous
  !> one ended (the windows are adjacent in grid index) and lies inside one window.
  !===============================================================================
  SUBROUTINE verify_output_patch_coverage(output_sp, n_nu_out)
    class(atomic_spectrum_type), intent(in) :: output_sp
    integer, intent(in) :: n_nu_out
    integer :: i_nu, i_seg, fs, fe, prev_fe
    !---------------------------------------------------------------------------
    if (output_sp%n_segments == 0) return

    if (mod(n_nu_out, output_sp%n_segments) /= 0) then
      write(stderr,'(A,I0,A,I0)') &
        'ERROR: n_nu_out=', n_nu_out, ' is not a multiple of n_segments=', output_sp%n_segments
      error stop 'n_nu_out and n_segments mismatch'
    endif

    prev_fe = 0
    do i_nu = -1, -n_nu_out, -1
      call output_patch_range(output_sp, i_nu, n_nu_out, fs, fe)
      i_seg = output_patch_window(output_sp, i_nu, n_nu_out)
      if (fs /= prev_fe + 1 .or. fe < fs .or. &
          fs < output_sp%segment_first(i_seg) .or. fe > output_sp%segment_last(i_seg)) then
        write(stderr,'(A,I0,A,I0,A,I0,A,I0)') 'ERROR: output patch i_nu=', i_nu, &
          ' covers [', fs, ',', fe, '] in window ', i_seg
        error stop 'output patches do not tile the windows'
      endif
      prev_fe = fe
    enddo
    if (prev_fe /= output_sp%nnu) then
      write(stderr,'(A,I0,A,I0)') 'ERROR: output patches end at ', prev_fe, ' of ', output_sp%nnu
      error stop 'incomplete segment coverage'
    endif

    write(stdout,'(A,I0,A,I0,A)') 'Verified: ', n_nu_out, ' output patches tile ', &
      output_sp%n_segments, ' windows'
  END SUBROUTINE verify_output_patch_coverage

END MODULE rt_patch_mod
