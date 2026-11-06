!===============================================================================
!> The RT slabs, which compute the opacities, formal solution, mean intensity and rates of one window, direction and depth slab.
!===============================================================================
MODULE rt_patch_mod
  USE os_mod
  USE io_mod
  USE scaling_mod
  USE trace_mod
  USE bits_mod
  USE const_mod, only: CC, h_k, twoh_c2, half, one_third, one_sixth, HZ_TO_AA, pi
  USE patch_mod, only: patch_t
  USE interpol_mod, only: pchip_mono
  USE handyfuncs_mod
  USE atom_mod, only: atom_type, atom_data_t, atom_trans_t, atomic_spectrum_type, radiative_cut
  USE ray_patch_mod, only: ray_patch_t, output_patch_window, output_patch_range, window_range
  USE m3d_halo_mod, only: m3d_halo_t, lock_mem, unlock_mem
  USE m3d_patch_mod, only: m3d_patch_t, point_opac_set_t, bg_opac_t, &
    dlam, out_patch_offset, &
    premature_finish, m3d_patch_pre_update=>pre_update, save_resolved, &
    column_by_column, rotate_atmos, save_cntrbf, depression_cntrbf, &
    decouple_continuum, initial_Jnu,  &
    rt_to_file, preconditioning, rate_formulation, &
    scattering_preconditioning, scat_op_max, &
    make_eos, amr, save_Inu, save_chi, save_Snu, &
    fixed_Jnu, n_windows, save_Jnu, &
    convlim, flux_convergence, flux_conv_signal, &
    global_flux_err, global_flux_pos, global_flux_lambda, flux_update_counter, &
    maser_chi_floor, maser_dtau_max, maser_Inu_max, suppress_masers, set_file
  USE mpi_mod, only: mpi
  USE rt_params_mod, only: rt_params_t, rt_iter, rt_out, ray_offset, pivot_z, &
                           save_linelist
  USE ray_mod, only: one_dimensional, ray_type
  USE trans_mod, only: line_type, cont_type
  USE qlg, only: qlg_bot_profile      ! the one Voigt kernel; see linelist_mod's kernel note
  USE bftrans_mod, only: H_occupation_probability
  USE hprofl_mod, only: hprofile, hline_t, hydrogen_line, hcell_t, hcell_prepare
  USE hvcs_mod, only: vcs_profile_t, vcs_supported
  USE m3d_debug_mod
  USE shared_mod
  USE composition_mod, only: composition_type, element_ids, &
    departure_pointer_t, read_departures, &
    line_scattering
  USE osmet_mod, only: compute_epsilon_stage
  USE spectrum_mod, only: spectrum_type
  USE optab_mod, only: sample_optab_type
  USE omp_lock_mod, only: lock_t
  USE console_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer :: n_negative_chi = 0                 !< diagnostic: slabs that reported a floored chi in work_to_transfer_units

  !-----------------------------------------------------------------------------
  ! The kind of a set, as opposed to its position.  The set axis has a variable
  ! length (an iteration slab has a single main set and no continuum set), so the
  ! kind is recorded per set and written to 'sets.txt'; readers look the layout
  ! up instead of inferring it from the set index or the number of sets.
  !-----------------------------------------------------------------------------
  integer, parameter, public :: set_continuum=0 !< background only; skips foreground
  integer, parameter, public :: set_nlte=1      !< the main spectrum, with departures
  integer, parameter, public :: set_lte=2       !< LTE: the main spectrum of an LTE run, or the LTE companion of an NLTE one
  character(len=9), dimension(0:2), parameter :: set_kind_name = &
    ['continuum', 'nlte     ', 'lte      ']                     !< names of the set kinds in 'sets.txt'
  integer, public, save :: output_n_set=0                       !< number of output sets
  integer, public, save :: output_i_set_main=0                  !< the output set the run's own spectrum comes from
  integer, public, save :: output_i_set_cnt=0                   !< the unperturbed continuum output set
  integer, allocatable, public, save :: output_set_kind(:)      !< kind of each output set
  real(8), allocatable, public, save :: output_set_da(:,:)      !< abundance offset [dex] per (atomic number, set)
  logical, allocatable, public, save :: output_set_element(:,:) !< set_file names the element for the set, including da=0

  !-----------------------------------------------------------------------------
  ! Output patches use IDs -n_windows_out:-1, internal patches use 1:n_windows.  n_windows_out
  ! and the patches of each output window (out_patch_offset) are set by set_n_windows_out
  ! in composition_patch_mod.
  !-----------------------------------------------------------------------------
  public :: verify_output_patch_coverage, window_shape, slab_memory, trans_copy_bytes, &
            point_table_bytes

  !===============================================================================
  !> What the memory of one slab depends on besides its cells: the frequency
  !> points, directions and sets of its window and the transitions of each atom
  !> that reach into it.  window_shape predicts it before the slab exists.
  !===============================================================================
  type, public :: slab_shape_t
    integer :: nf=0                      !< frequency points of the window
    integer :: n_mu=0                    !< directions of the slab
    integer :: n_set=1                   !< sets of the slab
    integer, allocatable :: nline(:)     !< lines of each atom that reach into the window
    integer, allocatable :: ncont(:)     !< continua of each atom that reach into the window
    integer, allocatable :: nhtab(:)     !< hydrogen lines with a tabulated profile
    logical, allocatable :: get_rates(:) !< the slab accumulates the rates of each atom
    logical :: bound_bound=.false.       !< some atom has lines
    logical :: bound_free=.false.        !< some atom has continua
    logical :: do_final=.false.          !< a slab of the final output pass
    logical :: do_out_sp=.false.         !< a slab working on the output spectrum
  end type slab_shape_t

  !===============================================================================
  !> The tabulated profiles of an atom's hydrogen lines (hcell_t%vcs), prepared once
  !> per cell and line and shared by the directions, sets and passes of one execute.
  !===============================================================================
  type :: hvcs_store_t
    integer, allocatable :: slot(:)                !< index of each line in p, 0 for a line without a table
    type(vcs_profile_t), allocatable :: p(:,:,:,:) !< prepared profiles per (x, y, z, slot)
  end type hvcs_store_t

  !===============================================================================
  !> The state of one frequency window, ray direction and depth slab.  A slab is
  !> not a task: the rt_column_t of its window and direction sweeps it in order
  !> (rt_column_mod), sets upstream to the slab swept before it, hands it the
  !> boundary buffer (mem) and marks the last slab of a complete sweep with
  !> column_end.  up_twin and dn_twin are set there too.
  !===============================================================================
  type, public, extends(ray_patch_t) :: rt_patch_t
    class(rt_patch_t), pointer :: upstream=>null()              !< the slab swept before this one
    class(rt_patch_t), pointer :: up_twin=>null()               !< itself if twinned; final pass: the output pair's i_mu=1 slab
    class(rt_patch_t), pointer :: dn_twin=>null()               !< its twin; final pass: the output pair's i_mu=2 slab
    logical :: Jnu_freed=.false.                                !< finish has freed the twins' Jnu
    logical :: column_end=.false.                               !< the last slab of a complete sweep of its column
    logical :: settled=.false.                                  !< atoms converged here and below (rt_patch_pre_update)
    integer :: hemisphere_iter=-1                               !< last iteration whose hemisphere this slab has completed
    integer :: combined_iter=-1                                 !< iteration of the full-sphere Jnu held here; -1 before the first
    logical :: PsiJ_combined=.false.                            !< the last twin combine averaged PsiJ too
    class(atom_trans_t), dimension(:), pointer :: trans=>null() !< lines and continua of each atom that reach into the window

    type(spectrum_type) :: sub_sp                     !< frequencies of the window
    real(4), dimension(:,:,:,:,:), allocatable :: Jnu !< mean intensity (occupation number) per (f, x, y, z, set)
    real(4), dimension(:,:,:,:), allocatable :: Psi   !< diagonal of the approximate operator of the current direction

    !> Per frequency and cell: the largest single contribution to the emissivity
    !> (index 1) and the sum of all the others (index 2).  Holding the largest
    !> apart is what lets a line's own contribution be taken out of the local
    !> source term without a cancelling subtraction.  Contribution 0 is everything
    !> that is not a model-atom line (continuum, bound-free and line list), seeded
    !> as one contribution in foreground.  It is stored without the
    !> stimulated-emission factor beta, which is common to all contributions at a
    !> frequency and cell and is applied where they are read.  The line-scattering
    !> emissivity is not tracked: it is added to eta in work_to_transfer_units and
    !> to the other contributions in the rates.
    real(4), allocatable :: split_opac(:,:,:,:,:)
    real(4), allocatable :: split_chi(:,:,:,:,:)  !< absorption held the same way, so chi - chi_t is never a subtraction
    integer, allocatable :: split_id(:,:,:,:)     !< contribution index of the largest emissivity
    integer, allocatable :: split_chi_id(:,:,:,:) !< contribution index of the largest absorption
    integer, allocatable :: split_offset(:)       !< offset of each atom's line indices among the contributions
    type(hvcs_store_t), allocatable :: hvcs(:)    !< prepared tabulated hydrogen profiles of each atom
    real(4), allocatable :: Irem(:,:,:,:)         !< formal solution without its own local term, Inu - Psi*eta_total
    real(4), allocatable :: escape(:,:,:,:)       !< escape factor 1 - Psi*chi of the local term
    logical :: split_rates=.false.                !< the contributions are tracked (rate_formulation='largest')

    !> Angle-averaged diagonal operator per (f, x, y, z, set).  Like Jali it
    !> persists between iterations as Jnu does, because the Olson-Auer-Buchler
    !> update of the continuum-scattering ALI is a two-term recurrence that needs
    !> the previous iteration's values at the twin combine.
    real(4), dimension(:,:,:,:,:), allocatable :: PsiJ
    real(4), dimension(:,:,:,:,:), allocatable :: Jali !< accelerated full-sphere mean intensity of the scattering ALI
    logical, allocatable :: Jali_valid(:)              !< per set, a previous Jali is stored (else the recurrence starts from Jnu)
    logical :: do_scat_precond = .false.               !< the continuum-scattering ALI is active for this slab
    logical :: compute_psi = .false.                   !< the formal solvers build Psi (for the rates or the ALI)

    !> Sampled background opacity table: the composition patch's shared table for
    !> this i_nu or, on the make_eos path, a private one (save_opac_table
    !> overwrites it in place).
    type(sample_optab_type), pointer :: opac => null()
    type(lock_t), pointer :: opac_lock => null() !< lock of the shared opac table
    !> Sampled LTE line-list table (sample_linelist): non-departure species only,
    !> interpolated into the foreground work arrays with a per-cell Doppler shift;
    !> the i_mu twins share one instance per i_nu, like opac.
    type(sample_optab_type), pointer :: ll_opac => null()
    type(lock_t), pointer :: ll_opac_lock => null() !< lock of the shared ll_opac table

    real(4), dimension(:), allocatable :: twohnu3_c2   !< 2 h nu^3/c^2 per frequency [erg/cm^2]
    real(8), dimension(:), pointer :: nu               !< frequencies of the window [Hz]
    integer :: mus                                     !< first direction of the slab
    integer :: mue                                     !< last direction of the slab
    integer :: n_set                                   !< number of sets
    integer :: n_mu                                    !< number of directions
    integer :: natom                                   !< number of model atoms
    integer, dimension(:), allocatable :: set_kind     !< set_continuum, set_nlte or set_lte
    integer :: i_set_main=1                            !< the set the run's own spectrum comes from
    integer :: i_set_cnt=0                             !< the continuum set, or 0 where none is computed
    logical :: has_nlte_set=.false.                    !< there is a set carrying departures
    logical :: bound_bound                             !< some atom has lines
    logical :: bound_free                              !< some atom has continua
    logical, dimension(:), allocatable :: get_rates    !< per atom, this slab accumulates its rates
    integer :: nx                                      !< cells in x
    integer :: ny                                      !< cells in y
    integer :: nz                                      !< cells in z
    real(4), dimension(:,:,:,:,:), allocatable :: work !< (f, x, y, z, k): k=1 extinction, 2 source, 3 intensity

    ! Direction-independent background, per (f,x,y,z), rebuilt once per i_set in
    ! background().  Filled in this order:
    !   lte_background_opacities  chi_cont, sigma  (from the table or comp%opacity)
    !                             beta = exp(-h*nu/kT), clamped off zero, and
    !                             Bnu = beta/(1-beta)
    !   lte_background_emission   ehat = chi_cont*Bnu + sigma*(Jnu - Bnu)   [= eta]
    !   convert_to_work_units     xhat = chi_cont + ehat, then ehat -> ehat/beta
    !   add_bound_free_opacities  accumulates into xhat and ehat
    ! so xhat/ehat hold the emissivity eta in *physical* units until
    ! convert_to_work_units runs, and the working pair afterwards.  back_to_work
    ! then seeds work(...,1:2) from xhat/ehat for each direction.
    real(4), dimension(:,:,:,:), allocatable :: chi_cont !< continuum extinction, total
    real(4), dimension(:,:,:,:), allocatable :: sigma    !< coherent-scattering part of chi_cont (<= it)
    real(4), dimension(:,:,:,:), allocatable :: beta     !< LTE Boltzmann factor exp(-h*nu/kT)
    real(4), dimension(:,:,:,:), allocatable :: Bnu      !< Planck function as occupation number, beta/(1-beta)
    real(4), dimension(:,:,:,:), allocatable :: xhat     !< working extinction, chi_cont + eta plus the bound-free terms
    real(4), dimension(:,:,:,:), allocatable :: ehat     !< emissivity eta, then eta/beta (working source)

    logical :: do_line_scat=.false.                             !< this pass filled chi_scat and dJnu
    real(4), dimension(:,:,:,:), allocatable :: chi_scat        !< (1-eps)-weighted LTE line opacity of the current direction
    real(4), dimension(:,:,:,:), allocatable :: dJnu            !< Jnu - Bnu as the background left it
    real(4), dimension(:,:,:,:,:), allocatable :: one_m_eps     !< scattering fraction 1-eps per stage (1 neutral, 2 ionised)
    real(4), dimension(:,:,:,:,:,:), allocatable :: cntrbf_data !< contribution function and extinction per direction
    real(4), dimension(:,:,:,:,:,:), allocatable :: cntrbf_cnt  !< chi, eta and Inu of the continuum set per direction
    real(8), dimension(:), allocatable :: flux_prev             !< previous iteration's emergent flux per frequency
    real(8) :: flux_global=0.d0                                 !< largest relative flux change so far in this round
  contains
    procedure :: pre_init => rt_patch_pre_init
    procedure :: init_patch_trans
    procedure :: pre_update => rt_patch_pre_update
    procedure :: execute
    procedure :: estimate_memory
    procedure :: check_flux_convergence
    procedure :: background
    procedure :: foreground
    procedure :: aniso_integrations
    procedure :: iso_integrations
    procedure :: has_twin
    procedure :: restore_combined_state
    procedure :: load_Jnu_set
    procedure :: combine_twin_hemispheres
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
    procedure :: long_char_quadratic
    procedure :: long_char_hermite
    procedure :: long_char_hermite2
    procedure :: short_char_hermite
    procedure :: get_vproj
    procedure :: free_point_tables
    procedure :: print_post_update_info
    procedure :: append_to_task_list
    procedure :: finish
    procedure :: allocate_mem
    procedure :: allocate_boundary
    procedure :: retain_exit_plane
    procedure :: restore_exit_plane
    procedure :: release_exit_plane
  end type rt_patch_t

  !===============================================================================
  !> A pointer to an RT slab, for arrays of slabs.
  !===============================================================================
  type, public :: rt_patch_ptr_t
    class(rt_patch_t), pointer :: p => null() !< the slab
  end type rt_patch_ptr_t

  !-----------------------------------------------------------------------------
  ! NB: work, the background arrays, cntrbf_data and the patch dimensions
  ! nx,ny,nz are components of rt_patch_t (owned per patch); the loop indices
  ! x,y,z,f are plain locals in each routine.  Keep them out of module
  ! threadprivate storage: on macOS/GCC that routes every access through emulated
  ! TLS (pthread_getspecific).  One thread per patch makes instance/local data
  ! both race-free and TLS-free.

  !===============================================================================
  !> Formal solvers, implemented in the rt_solvers_mod submodule (rt_solvers_mod.f90)
  !===============================================================================
  interface
    !=============================================================================
    !> Linear formal solution by short characteristics in a 3D slab, which builds
    !> Psi for the iteration pass.
    !=============================================================================
    module subroutine short_char (self, vec)
      class(rt_patch_t) :: self
      real(8), intent(in) :: vec(3) !< direction of the ray
    end subroutine short_char

    !=============================================================================
    !> Linear formal solution along the vertical columns of the slab (1D,
    !> column_by_column or a ray-aligned atmosphere), which builds Psi.
    !=============================================================================
    module subroutine long_char_linear (self, vec)
      class(rt_patch_t) :: self
      real(8) :: vec(3) !< direction of the ray
    end subroutine long_char_linear

    !=============================================================================
    !> Parabolic transfer moments over an interval of optical depth t.
    !=============================================================================
    pure module subroutine parab_moments (t, E, U0, U1, q)
      real(4), intent(in)  :: t  !< optical depth of the interval
      real(4), intent(out) :: E  !< exp(-t)
      real(4), intent(out) :: U0 !< 1 - E
      real(4), intent(out) :: U1 !< t - U0
      real(4), intent(out) :: q  !< (t^2 - 2*U1)/t
    end subroutine parab_moments

    !=============================================================================
    !> Quadratic formal solution along the vertical columns, with a limited
    !> source slope.
    !=============================================================================
    module subroutine long_char_quadratic (self, vec)
      class(rt_patch_t), target :: self
      real(8) :: vec(3) !< direction of the ray
    end subroutine long_char_quadratic

    !=============================================================================
    !> Cubic Hermite formal solution along the vertical columns, with the slopes of
    !> chi and Snu from the composition gradients by the chain rule.
    !=============================================================================
    module subroutine long_char_hermite (self, comp, vec, comp_up)
      class(rt_patch_t), target :: self
      class(composition_type) :: comp             !< composition of the slab
      class(composition_type), pointer :: comp_up !< composition of the upstream slab, null at the domain boundary
      real(8), intent(in) :: vec(3)               !< direction of the ray
    end subroutine long_char_hermite

    !=============================================================================
    !> long_char_hermite's optical depth with long_char_quadratic's source
    !> treatment, which needs only the Pe gradients.
    !=============================================================================
    module subroutine long_char_hermite2 (self, comp, vec, comp_up)
      class(rt_patch_t), target :: self
      class(composition_type) :: comp             !< composition of the slab
      class(composition_type), pointer :: comp_up !< composition of the upstream slab, null at the domain boundary
      real(8), intent(in) :: vec(3)               !< direction of the ray
    end subroutine long_char_hermite2

    !=============================================================================
    !> short_char with long_char_hermite's chain-rule slopes, which builds Psi.
    !=============================================================================
    module subroutine short_char_hermite (self, comp, vec)
      class(rt_patch_t) :: self
      class(composition_type) :: comp !< composition of the slab
      real(8), intent(in) :: vec(3)   !< direction of the ray
    end subroutine short_char_hermite
  end interface

contains

  !===============================================================================
  !> Leading registered RT allocations of this slab, in bytes, before its first
  !> execution; see slab_memory.  u_phase counts until the slab allocates it.
  !> tables are the point tables of its window and depth slab (point_table_bytes).
  !===============================================================================
  SUBROUTINE estimate_memory (self, work, retained, boundary, persistent, tables)
    class(rt_patch_t), intent(in) :: self
    real(8), intent(out) :: work       !< execution-local bytes
    real(8), intent(out) :: retained   !< bytes kept between executions
    real(8), intent(out) :: boundary   !< bytes of the boundary buffer
    real(8), intent(out) :: persistent !< bytes of u_phase, 0 once it is allocated
    real(8), intent(out) :: tables     !< bytes of the point tables
    type(slab_shape_t) :: shape
    integer :: i
    !-----------------------------------------------------------------------------
    shape%nf = self%nf
    shape%n_mu = self%n_mu
    shape%n_set = self%n_set
    shape%bound_bound = self%bound_bound
    shape%bound_free = self%bound_free
    shape%do_final = self%do_final
    shape%do_out_sp = self%do_out_sp
    allocate(shape%nline(self%natom), shape%ncont(self%natom), shape%nhtab(self%natom))
    do i=1, self%natom
      shape%nline(i) = size(self%trans(i)%line)
      shape%ncont(i) = size(self%trans(i)%cont)
      shape%nhtab(i) = count(table_line_slots(self%trans(i)%line) > 0)
    end do
    allocate(shape%get_rates, source=self%get_rates)
    call slab_memory (shape, self%n, self%ng, work, retained, boundary, persistent)
    if (associated(self%u_phase)) persistent = 0.d0
    tables = point_table_bytes(shape, self%n)
  END SUBROUTINE estimate_memory

  !===============================================================================
  !> Shape of the slab of window i_nu and direction i_mu that create_sub_patches
  !> would make from spectrum sp and atoms, as rt_patch_pre_init and
  !> init_patch_trans set it up.  The window counts n_windows and n_windows_out
  !> must be set.
  !===============================================================================
  FUNCTION window_shape (sp, atoms, i_nu, i_mu) RESULT (shape)
    class(atomic_spectrum_type), intent(in) :: sp    !< the spectrum the windows divide
    type(atom_type), pointer, intent(in) :: atoms(:) !< model atoms, may be unassociated
    integer, intent(in) :: i_nu                      !< window, negative for an output window
    integer, intent(in) :: i_mu                      !< direction, negative for the final pass
    type(slab_shape_t) :: shape
    logical :: is_nlte
    integer :: i, natom, fs, fe
    !-----------------------------------------------------------------------------
    shape%do_final = i_mu < 0
    shape%do_out_sp = i_nu < 0 .or. i_mu < 0
    shape%n_mu = merge(size(rt_out%wmu), size(rt_iter%wmu)/2, shape%do_final)
    if (shape%do_out_sp) then
      is_nlte = rt_iter%maxiter > 0 .or. rt_out%maxiter > 0 .or. read_departures
      call init_output_sets(is_nlte)
      shape%n_set = output_n_set
    endif
    call window_range(sp, i_nu, fs, fe)
    shape%nf = fe - fs + 1
    natom = 0
    if (associated(atoms)) natom = size(atoms)
    allocate(shape%nline(natom), shape%ncont(natom), shape%nhtab(natom))
    if (natom > 0) then
      shape%bound_bound = any([(associated(atoms(i)%trans%line), i=1, natom)])
      shape%bound_free = any([(associated(atoms(i)%trans%cont), i=1, natom)])
    endif
    do i=1, natom
      call count_window_trans(atoms(i), shape%bound_bound, shape%bound_free, fs, fe, &
                              shape%nline(i), shape%ncont(i), shape%nhtab(i))
    end do
    allocate(shape%get_rates(natom), source=(.not. shape%do_out_sp))
  END FUNCTION window_shape

  !===============================================================================
  !> Bytes of the transition copies (init_patch_trans) of a slab of this shape,
  !> registered when the slab is made.  A continuum usually spans many windows,
  !> so its copies grow with n_windows.
  !===============================================================================
  FUNCTION trans_copy_bytes (shape) RESULT (bytes)
    type(slab_shape_t), intent(in) :: shape !< shape of the slab
    real(8) :: bytes
    type(line_type) :: line
    type(cont_type) :: cont
    !-----------------------------------------------------------------------------
    bytes = (real(sum(shape%nline),8)*storage_size(line) + &
             real(sum(shape%ncont),8)*storage_size(cont))/8
  END FUNCTION trans_copy_bytes

  !===============================================================================
  !> Bytes of the 1D point tables a depth slab of n cells holds for the window of
  !> this shape: the LTE background (bg_point_tables, chi_cont and sigma) and, in
  !> the output pass, the detailed line opacity of every set with lines
  !> (ll_point_tables, chi, eta and with line scattering chi_scat).  The slabs of
  !> a window share them, so they count once per window and depth slab.  An
  !> output window frees its tables after its final pass; without output Lambda
  !> iterations they live only within that one execution and are not counted.  A
  !> moving 1D atmosphere never fills the line-opacity store, which is counted
  !> all the same.
  !===============================================================================
  FUNCTION point_table_bytes (shape, n) RESULT (bytes)
    type(slab_shape_t), intent(in) :: shape !< shape of the slab
    integer, intent(in) :: n(3)             !< cells of the depth slab
    real(8) :: bytes, field
    !-----------------------------------------------------------------------------
    bytes = 0.d0
    if (.not. one_dimensional) return
    if (shape%do_out_sp .and. rt_out%maxiter <= 0) return
    field = real(storage_size(0.0_4)/8,8)*product(real(n,8))*shape%nf
    if (.not. make_eos) bytes = 2*field
    if (shape%do_out_sp .and. rt_out%include_linelist .and. .not. rt_out%sample_linelist) &
      bytes = bytes + merge(3, 2, line_scattering)*field*count(output_set_kind /= set_continuum)
  END FUNCTION point_table_bytes

  !===============================================================================
  !> Leading registered RT allocations of a slab of n cells and ng ghost cells, in
  !> bytes, before the first execution.  Keep the array counts consistent with
  !> allocate_arrays and deallocate_arrays.  work is execution-local; retained
  !> survives between executions when rt_to_file is false; boundary is the
  !> buffer its column hands from slab to slab; persistent is the slab's
  !> u_phase.  The point tables (point_table_bytes), the sampled opacity tables
  !> and solver temporaries are excluded.
  !===============================================================================
  SUBROUTINE slab_memory (shape, n, ng, work, retained, boundary, persistent)
    type(slab_shape_t), intent(in) :: shape !< shape of the slab
    integer, intent(in) :: n(3)             !< cells of the slab
    integer, intent(in) :: ng(3)            !< ghost cells of the slab
    real(8), intent(out) :: work            !< execution-local bytes
    real(8), intent(out) :: retained        !< bytes kept between executions
    real(8), intent(out) :: boundary        !< bytes of the boundary buffer
    real(8), intent(out) :: persistent      !< bytes of u_phase
    real(8) :: word, cell, field, sets
    logical :: scat_ali, psi, rates
    type(vcs_profile_t) :: hvcs_sample
    integer :: i
    !-----------------------------------------------------------------------------
    word = real(storage_size(0.0_4)/8,8)
    cell = word*product(real(n,8))
    field = cell*shape%nf
    sets = field*shape%n_set
    ! chi_cont, sigma, beta, Bnu, xhat, ehat, and work(...,1:3) with its ghost zones
    work = 6*field + 3*word*shape%nf*n(1)*n(2)*(real(n(3),8)+2*ng(3))
    retained = 0.d0
    boundary = 3*word*shape%nf*n(1)*n(2)*shape%n_mu*shape%n_set
    ! u_phase, allocated by patch_t%init when the patch first executes
    persistent = word*product(real(n+2*ng,8))
    rates = any(shape%get_rates)

    if (rt_to_file) then
      work = work + sets
    else
      retained = retained + sets
    endif
    do i=1, size(shape%nline)
      if (shape%bound_bound) work = work + cell*shape%nline(i)
      ! hvcs, the prepared tabulated hydrogen profiles
      if (shape%bound_bound .and. allocated(shape%nhtab)) work = work + &
        product(real(n,8))*shape%nhtab(i)*storage_size(hvcs_sample)/8
      if (shape%get_rates(i)) work = work + 2*cell*(shape%nline(i) + shape%ncont(i))
    enddo
    ! split_opac(2), split_id, split_chi(2), split_chi_id, Irem and escape
    if (rate_formulation == 'largest' .and. rates) work = work + 8*field
    if (line_scattering) work = work + 4*field
    ! the do_scat_precond and compute_psi tests of execute
    scat_ali = scattering_preconditioning .and. .not. shape%do_final &
      .and. (rates .or. shape%do_out_sp) &
      .and. (rt_iter%maxiter > 0 .or. rt_out%maxiter > 0)
    psi = rates .or. scat_ali
    if (psi) then
      work = work + field
      if (shape%bound_free .or. scat_ali) then
        if (scat_ali .and. .not. rt_to_file) then
          retained = retained + sets
        else
          work = work + sets
        endif
      endif
      if (scat_ali) then
        if (rt_to_file) then
          work = work + sets
        else
          retained = retained + sets
        endif
      endif
    endif
    if (shape%do_final .and. save_cntrbf) work = work + 5*field*shape%n_mu
  END SUBROUTINE slab_memory

  !===============================================================================
  !> Sets the slab up before its first execution: its directions and sets, the
  !> layout of its boundary buffer, the frequencies and transitions of its window
  !> and the atoms whose rates it accumulates.
  !===============================================================================
  SUBROUTINE rt_patch_pre_init (self)
    class(rt_patch_t), target:: self
    integer, save :: itimer=0
    logical :: is_nlte
    ! ----------------------------------------------------------------------------
    call trace%begin('rt_patch_t%pre_init', itimer=itimer)
    self%kind = 'rt_patch'
    call self%m3d_patch_t%pre_init()

    self%dir = merge(1, -1, abs(self%i_mu)==1)
    self%n_windows = n_windows
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
  !> Give this output slab the process-global set layout.
  !===============================================================================
  SUBROUTINE init_set_layout (self, is_nlte)
    class(rt_patch_t) :: self
    logical, intent(in) :: is_nlte !< the run has departure coefficients
    !-----------------------------------------------------------------------------
    call init_output_sets(is_nlte)
    self%n_set = output_n_set
    allocate(self%set_kind(self%n_set), source=output_set_kind)
    self%i_set_cnt = output_i_set_cnt
    self%i_set_main = output_i_set_main
    if (self%n_set > merge(3,2,is_nlte)) then
      if (rt_out%sample_linelist) &
        error stop 'abundance ladders require rt_output_params sample_linelist=.false.'
    endif
  END SUBROUTINE init_set_layout

  !===============================================================================
  !> Sets up the output set layout once per process: a continuum set, the main
  !> set and, with departures, an LTE set or, with set_file, the sets of an
  !> abundance-ladder manifest.  Each composition instance also reads the
  !> process-global layout once.
  !===============================================================================
  SUBROUTINE init_output_sets (is_nlte)
    logical, intent(in) :: is_nlte !< the run has departure coefficients
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
          call console_error ('output_params set_file must be outside the run directory')
          error stop 'invalid set_file location'
        endif
        open(newunit=unit, file=trim(set_file), status='old', action='read', iostat=ios)
        if (ios /= 0) then
          call console_error ('cannot open set_file: '//trim(set_file))
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
            call console_error ('unknown kind in set_file: '//trim(kind))
            error stop 'invalid set kind'
          end select
          if (seen(number) .and. output_set_kind(number) /= i) &
            error stop 'one set number has conflicting kinds'
          output_set_kind(number) = i
          seen(number) = .true.
          if (trim(element) == '-') then
            if (trim(da_word) /= '-') error stop 'set_file has da without element'
          else
            z = findloc(element_ids, adjustl(element), dim=1)
            if (z == 0) then
              call console_error ('unknown element in set_file: '//trim(element))
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
  END SUBROUTINE init_output_sets

  !===============================================================================
  !> Write 'sets.txt': one row per (set, element), naming what each slot on the
  !> saved set axis of 'ie' and 'flux' contains, with the element's abundance
  !> offset; a set without offsets is one row with '-' in both columns.  Written
  !> once per run, by the first output slab to execute.
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
      open(newunit=lunit, file=trim(io%run_dir)//'sets.txt', form='FORMATTED')
      write(lunit, '(A)') '# set  kind       element  da'
      do i_set=1, self%n_set
        c_element = '-'
        c_da = '-'
        if (any(output_set_element(:,i_set))) then
          do z=1,92
            if (.not. output_set_element(z,i_set)) cycle
            write(c_element, '(A)') trim(element_ids(z))
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
  END SUBROUTINE save_set_manifest

  !===============================================================================
  !> Copies the active lines and continua of each atom that reach into the window
  !> to trans, with their frequency ranges clipped to the window.
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
        call count_window_trans(atom, self%bound_bound, self%bound_free, self%fs, self%fe, &
                                nline, ncont)
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
      end associate
    end do
  END SUBROUTINE init_patch_trans

  !===============================================================================
  !> Number of active lines and continua of atom that reach into [fs_in, fe_in].
  !===============================================================================
  SUBROUTINE count_window_trans (atom, bound_bound, bound_free, fs_in, fe_in, nline, ncont, nhtab)
    type(atom_type), intent(in) :: atom     !< model atom
    logical, intent(in) :: bound_bound      !< count the lines
    logical, intent(in) :: bound_free       !< count the continua
    integer, intent(in) :: fs_in            !< first frequency of the window
    integer, intent(in) :: fe_in            !< last frequency of the window
    integer, intent(out) :: nline           !< number of lines
    integer, intent(out) :: ncont           !< number of continua
    integer, intent(out), optional :: nhtab !< number of hydrogen lines with a tabulated profile
    integer :: i, fs, fe, nh
    logical :: skip
    !-----------------------------------------------------------------------------
    nline = 0
    ncont = 0
    nh = 0
    if (bound_bound) then
      do i=1, atom%nline
        associate (line => atom%trans%line(i))
          if (.not. line%active) cycle
          call line%check_range(fs_in, fe_in, fs, fe, skip)
          if (skip) cycle
          nline = nline + 1
          if (line%hydrogen) then
            if (vcs_supported(line%i, line%j)) nh = nh + 1
          endif
        end associate
      end do
    endif
    if (bound_free) then
      do i=1, atom%ncont
        associate (cont => atom%trans%cont(i))
          if (.not. cont%active) cycle
          call cont%check_range(fs_in, fe_in, fs, fe, skip)
          if (.not. skip) ncont = ncont + 1
        end associate
      end do
    endif
    if (present(nhtab)) nhtab = nh
  END SUBROUTINE count_window_trans

  !===============================================================================
  !> The index of each line among the hydrogen lines with a tabulated profile, 0
  !> for the other lines (hvcs_store_t%slot).  The level indices of a hydrogen
  !> line are its principal quantum numbers.
  !===============================================================================
  FUNCTION table_line_slots (lines) RESULT (slot)
    type(line_type), intent(in) :: lines(:) !< lines of the window
    integer :: slot(size(lines))
    integer :: i, k
    !-----------------------------------------------------------------------------
    k = 0
    slot = 0
    do i=1, size(lines)
      if (.not. lines(i)%hydrogen) cycle
      if (.not. vcs_supported(lines(i)%i, lines(i)%j)) cycle
      k = k + 1
      slot(i) = k
    end do
  END FUNCTION table_line_slots

  !===============================================================================
  !> With premature_finish, an iteration slab freezes (static_mem) once it is
  !> settled: it sweeps once more and finishes.  Both columns then freeze from
  !> the bottom, the up column the first slabs of its sweep and the down column
  !> the last, so the two twins of a slab freeze in the same iteration, and only
  !> the up column has slabs after its frozen ones, which receive the exit plane
  !> of the highest frozen slab (rt_column_t%execute).  Population convergence
  !> does not end the task while flux_convergence decides.
  !===============================================================================
  SUBROUTINE rt_patch_pre_update (self)
    class(rt_patch_t), target :: self
    !-----------------------------------------------------------------------------
    call m3d_patch_pre_update(self)
    if (self%is_finished() .or. .not. (self%has_atom .and. self%do_full_sp)) return
    self%static_mem = self%settled .and. premature_finish .and. .not. flux_convergence &
      .and. self%iter > rt_out%maxiter-1
  END SUBROUTINE rt_patch_pre_update

  !===============================================================================
  !> One sweep step of the slab.  For each set and direction it adds the
  !> background and foreground opacities, solves the transfer and hands the exit
  !> plane downstream, then integrates Jnu and the rates over the directions.  The
  !> final pass solves each direction once, with rotate_atmos from the ray-aligned
  !> atmosphere; the iteration and output-scattering passes cover their hemisphere
  !> and combine it with their twin's.
  !===============================================================================
  SUBROUTINE execute (self)
    class(rt_patch_t) :: self
    class(atom_trans_t), pointer :: trans
    class(composition_type), pointer :: comp
    real(4), dimension(:,:,:), allocatable :: vproj
    real(4) :: wmu
    real(8) :: vec(3)
    logical :: get_cntrbf, do_rotate_atmos, top_of_domain
    integer :: i, i_mu, zs, ze, fs, fe, ng, i_set, nrad
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
    ! Every solver admitted for the iteration pass keeps the formal solution's
    ! non-local part alongside Psi, so the Rybicki & Hummer rates need no further
    ! condition.  fixed_Jnu needs none either: it replaces the formal solution
    ! altogether, leaving Psi at the zero it is allocated with, and
    ! bound_bound_rates then takes the unpreconditioned rates.
    self%split_rates = rate_formulation == 'largest' .and. any(self%get_rates)
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
      if (self%has_twin() .and. self%iter > 0) call self%restore_combined_state()
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
      if (self%has_twin()) call self%combine_twin_hemispheres()
    endif

    if (self%do_final) then
      call self%free_point_tables()
      if (save_linelist) call self%collect_line_diagnostics()
    endif

    ! mem holds the emergent intensity at the end of a complete sweep; the
    ! column frees it afterwards
    if (self%column_end) then
      if (self%do_final .and. self%n_mu>0) then
        call self%save_intensity()
      elseif (self%do_full_sp .and. top_of_domain) then
        call self%check_flux_convergence()
      endif
    endif

    call deallocate_arrays()

    ! a frozen slab has made its last sweep (rt_patch_pre_update)
    if (self%static_mem) call self%finish()

    call trace%end(itimer)

  contains

    !-----------------------------------------------------------------------------
    subroutine allocate_arrays()
      ! background arrays: see the component declarations for how they are filled
      allocate(self%chi_cont(fs:fe,self%nx,self%ny,self%nz), source=0.)
      allocate(self%sigma   (fs:fe,self%nx,self%ny,self%nz), source=0.)
      allocate(self%beta    (fs:fe,self%nx,self%ny,self%nz), source=0.)
      allocate(self%Bnu     (fs:fe,self%nx,self%ny,self%nz), source=0.)
      allocate(self%xhat    (fs:fe,self%nx,self%ny,self%nz), source=0.)
      allocate(self%ehat    (fs:fe,self%nx,self%ny,self%nz), source=0.)
      allocate(self%work(fs:fe,self%nx,self%ny,zs:ze,3), source=0.)
      if (self%split_rates) then
        allocate(self%split_opac(fs:fe,self%nx,self%ny,self%nz,2), source=0.)
        allocate(self%split_id(fs:fe,self%nx,self%ny,self%nz), source=0)
        allocate(self%split_chi(fs:fe,self%nx,self%ny,self%nz,2), source=0.)
        allocate(self%split_chi_id(fs:fe,self%nx,self%ny,self%nz), source=0)
        allocate(self%Irem(fs:fe,self%nx,self%ny,self%nz), source=0.)
        allocate(self%escape(fs:fe,self%nx,self%ny,self%nz), source=1.)
        allocate(self%split_offset(self%natom))
        nrad = 0
        do i=1, self%natom
          self%split_offset(i) = nrad
          nrad = nrad + size(self%trans(i)%line)
        end do
        call mem_register(self%split_opac, label='split_opac')
        call mem_register(self%split_chi, label='split_chi')
        call mem_register(self%Irem, label='Irem')
        call mem_register(self%escape, label='escape')
      endif
      allocate(self%twohnu3_c2(fs:fe))
      allocate(vproj(self%nx,self%ny,self%nz))
      self%twohnu3_c2 = twoh_c2 * self%nu**3
      call mem_register(self%chi_cont,  label='self%chi_cont')
      call mem_register(self%sigma,     label='self%sigma')
      call mem_register(self%beta,      label='self%beta')
      call mem_register(self%Bnu,       label='self%Bnu')
      call mem_register(self%xhat,      label='self%xhat')
      call mem_register(self%ehat,      label='self%ehat')
      call mem_register(self%work,      label='self%work')
      if (.not. allocated(self%Jnu)) then
        allocate(self%Jnu(fs:fe,self%nx,self%ny,self%nz,self%n_set), source=0.)
        call mem_register(self%Jnu, label='Jnu')
      endif
      if (self%bound_bound) allocate(self%hvcs(self%natom))
      do i=1, self%natom
        trans => self%trans(i)
        if (self%bound_bound) then
          call os%allocate(trans%adamp, [self%nx,self%ny,self%nz,size(trans%line)], label='adamp', zero=.true.)
          associate (store => self%hvcs(i))
            allocate(store%slot(size(trans%line)), source=table_line_slots(trans%line))
            if (any(store%slot > 0)) then
              allocate(store%p(self%nx,self%ny,self%nz,maxval(store%slot)))
              call os%mem_register(storage_size(store%p), size(store%p), 'hvcs')
            endif
          end associate
        endif
        if (self%get_rates(i)) then
          nrad = size(trans%line) + size(trans%cont)
          call os%allocate(trans%gamma_ij, [self%nx,self%ny,self%nz,nrad], label='gamma_ij', zero=.true.)
          call os%allocate(trans%gamma_ji, [self%nx,self%ny,self%nz,nrad], label='gamma_ji', zero=.true.)
        endif
      end do
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
        ! rt_to_file, where persistence goes through the file instead.
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
        ! (chi, eta, Inu) of the continuum set, per direction (see
        ! save_contribution_function)
        allocate(self%cntrbf_cnt(fs:fe,self%nx,self%ny,self%nz,3,self%mus:self%mue))
        call mem_register(self%cntrbf_cnt, 'self%cntrbf_cnt')
      endif
    end subroutine allocate_arrays

    !-----------------------------------------------------------------------------
    subroutine deallocate_arrays()
      if (allocated(self%split_opac)) then
        call mem_unregister(self%split_opac, label='split_opac')
        call mem_unregister(self%split_chi, label='split_chi')
        call mem_unregister(self%Irem, label='Irem')
        call mem_unregister(self%escape, label='escape')
        deallocate(self%split_opac, self%split_id, self%split_chi, self%split_chi_id, &
                   self%split_offset, self%Irem, self%escape)
      endif
      call mem_unregister(self%work     , label='self%work')
      call mem_unregister(self%chi_cont , label='self%chi_cont')
      call mem_unregister(self%sigma    , label='self%sigma')
      call mem_unregister(self%beta     , label='self%beta')
      call mem_unregister(self%Bnu      , label='self%Bnu')
      call mem_unregister(self%xhat     , label='self%xhat')
      call mem_unregister(self%ehat     , label='self%ehat')
      deallocate(self%work)
      deallocate(self%chi_cont, self%sigma, self%beta, self%Bnu, self%xhat, self%ehat)
      deallocate(self%twohnu3_c2)
      do i=1, self%natom
        trans => self%trans(i)
        if (self%bound_bound .and. allocated(trans%adamp)) then
          call os%deallocate(trans%adamp, label='adamp')
        endif
        if (allocated(self%hvcs)) then
          if (allocated(self%hvcs(i)%p)) &
            call os%mem_register(storage_size(self%hvcs(i)%p), -size(self%hvcs(i)%p), 'hvcs')
        endif
        if (self%get_rates(i)) then
          call os%deallocate(trans%gamma_ij, label='gamma_ij')
          call os%deallocate(trans%gamma_ji, label='gamma_ji')
        endif
      end do
      if (allocated(self%hvcs)) deallocate(self%hvcs)
      if (allocated(self%chi_scat)) then
        call mem_unregister(self%chi_scat, label='self%chi_scat')
        call mem_unregister(self%one_m_eps, label='self%one_m_eps')
        call mem_unregister(self%dJnu, label='self%dJnu')
        deallocate(self%chi_scat, self%one_m_eps, self%dJnu)
      endif
      if (allocated(self%Jnu) .and. rt_to_file) then
        call mem_unregister(self%Jnu, label='Jnu')
        deallocate(self%Jnu)
      endif
      if (self%compute_psi) then
        call mem_unregister(self%Psi, label='Psi')
        deallocate(self%Psi)
        ! PsiJ and Jali mirror the Jnu lifecycle: when the scattering ALI carries
        ! them to the next iteration's twin combine through memory they must stay
        ! alive across the execute; with rt_to_file persistence goes through the
        ! file instead, and without the ALI they are freed after each execute.
        if (allocated(self%PsiJ) .and. .not. (self%do_scat_precond .and. .not. rt_to_file)) then
          call mem_unregister(self%PsiJ, label='PsiJ')
          deallocate(self%PsiJ)
        endif
      endif
      if (allocated(self%Jali) .and. rt_to_file) then
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
  !> The direction-independent opacity and emissivity of set i_set: the LTE
  !> background (computed for the first set and reused), its scattering source
  !> from the previous Jnu, the bound-free opacities of the atoms and the damping
  !> of their lines.  Also resets the angle sums of Jnu and PsiJ.
  !===============================================================================
  SUBROUTINE background(self, comp, i_mu, i_set, do_rotate_atmos)
    class(rt_patch_t) :: self
    class(composition_type), pointer :: comp !< composition, ray-aligned with do_rotate_atmos
    integer, intent(in) :: i_mu              !< direction of a ray-aligned atmosphere, else 0
    integer, intent(in) :: i_set             !< set
    logical, intent(in) :: do_rotate_atmos   !< the composition is ray-aligned
    class(atom_data_t), pointer :: adata
    type(rt_params_t), pointer :: rtp
    real(8) :: vec(3)
    logical :: include_absmet, do_sample_opac, do_make_eos
    integer :: i
    !-----------------------------------------------------------------------------
    rtp => pass_params(self)
    ! with a decoupled continuum, only the bound-free half of the internal
    ! spectrum (i_nu > n_windows) can carry the absmet background
    include_absmet = rtp%include_absmet .and. (self%do_out_sp .or. &
                     self%i_nu > self%n_windows .or. .not. decouple_continuum)
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
    call self%lte_background_emission(comp, i_set, vec, do_rotate_atmos)
    call self%convert_to_work_units()
    call self%reset_array(self%Jnu(:,:,:,:,i_set:i_set))
    ! PsiJ persists across iterations for the scattering ALI, so zero it for a
    ! fresh angle accumulation here, after lte_background_emission has read the
    ! previous iteration's combined PsiJ for the OAB correction.
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
        ! the line damping is computed for the first set and reused (see above)
        if (self%bound_bound .and. i_set==1 .and. allocated(trans%adamp)) then
          call self%bound_bound_damping(comp, atom, trans, adata)
        endif
      end associate
    end do
  END SUBROUTINE background

  !===============================================================================
  !> Adds the direction-dependent line opacities of set i_set to the work arrays:
  !> the line list (detailed, stored or sampled) and the lines of the model atoms,
  !> Doppler-shifted by vproj.  A continuum set has none.
  !===============================================================================
  SUBROUTINE foreground(self, comp, i_mu, i_set, do_rotate_atmos, wmu, vproj)
    class(rt_patch_t) :: self
    class(composition_type), pointer :: comp       !< composition, ray-aligned with do_rotate_atmos
    integer, intent(in) :: i_mu                    !< direction
    integer, intent(in) :: i_set                   !< set
    logical, intent(in) :: do_rotate_atmos         !< the composition is ray-aligned
    real(4), intent(in) :: wmu                     !< quadrature weight of the direction
    real(4), intent(in), dimension(:,:,:) :: vproj !< velocity along the ray [c]
    class(atom_data_t), pointer :: adata
    class(departure_pointer_t), pointer :: departures(:)
    type(rt_params_t), pointer :: rtp
    logical :: use_linelist, sample_ll, do_scat, use_store
    integer :: i, i_seg, x, y, z, f
    real(4), allocatable :: scat_lte(:,:,:,:)
    real(dp), allocatable :: window_lam(:)
    !-----------------------------------------------------------------------------
    ! A continuum set carries the background only, in the scattering iterations as
    ! well as in the final pass: its Jnu, and with it the scattering term of the
    ! continuum spectrum, must come from the line-free problem.
    if (self%set_kind(i_set) == set_continuum) then
      self%do_line_scat = .false.
      return
    endif
    ! detailed line list vs sampled LTE line table, per pass (rt_iter or rt_out)
    rtp => pass_params(self)
    ! An output patch on a window_file grid visits the lines of its whole window, however
    ! many patches that window is split into; unallocated, window_lam is absent.
    if (self%do_out_sp .and. self%i_nu < 0 .and. self%sp%n_segments > 0) then
      i_seg = output_patch_window(out_patch_offset, self%i_nu)
      allocate(window_lam(2))
      window_lam = self%sp%lam([self%sp%segment_first(i_seg), self%sp%segment_last(i_seg)])
    endif
    use_linelist = comp%ll_initialised .and. rtp%include_linelist .and. .not. rtp%sample_linelist
    sample_ll = rtp%sample_linelist .and. associated(self%ll_opac)
    if (sample_ll) sample_ll = self%ll_opac%initialised

    ! line_scattering: accumulate the (1-eps)-weighted LTE line opacity, then
    ! redirect that fraction of the line source from Bnu to Jnu (see below)
    do_scat = line_scattering .and. (use_linelist .or. sample_ll)
    if (do_scat) self%chi_scat = 0.

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
          call comp%add_linelist_opacities(self%sub_sp, xhat, ehat, dlam, vproj, &
            departures=departures, only_departures=.true., i_set=i_set, &
            window_lam=window_lam)
        end associate
      endif
      ! (c) the LTE companion set: the table leaves the departure-carrying species
      ! out, so their lines are added in detail with b=1, exactly as the detailed
      ! line list does for this set (including the scattering split)
      if (self%set_kind(i_set) == set_lte .and. self%has_nlte_set .and. &
          associated(self%m3d_parent%departures)) then
        i = merge(i_mu+1, 1, do_rotate_atmos)
        associate(xhat=>self%work(:,:,:,1:self%nz,1), ehat=>self%work(:,:,:,1:self%nz,2))
          if (do_scat) then
            allocate(scat_lte, mold=self%chi_scat)
            scat_lte = 0.
            call compute_epsilon_stage(1, self%sub_sp%nu, comp%temp, comp%ne, &
              self%one_m_eps(:,:,:,:,1))
            call compute_epsilon_stage(2, self%sub_sp%nu, comp%temp, comp%ne, &
              self%one_m_eps(:,:,:,:,2))
            call comp%add_linelist_opacities(self%sub_sp, xhat, ehat, dlam, vproj, &
              only_departures=.true., lte_of=self%m3d_parent%departures(:,i), &
              one_m_eps=self%one_m_eps, chi_scat=scat_lte, beta=self%beta, &
              i_set=i_set, window_lam=window_lam)
            self%chi_scat = self%chi_scat + scat_lte
            deallocate(scat_lte)
          else
            call comp%add_linelist_opacities(self%sub_sp, xhat, ehat, dlam, vproj, &
              only_departures=.true., lte_of=self%m3d_parent%departures(:,i), &
              i_set=i_set, window_lam=window_lam)
          endif
        end associate
      endif
    endif

    ! rate_formulation='largest' tracks each model-atom line's share of the
    ! emissivity separately, so that a line's own contribution can be taken back
    ! out of the source term exactly rather than by subtraction.  Everything that
    ! is not such a line (continuum, bound-free and the line list) enters as
    ! one contribution with index 0, which is why the line list is accumulated
    ! first and the model atoms only afterwards.  Seeding the contribution from
    ! the finished non-atom emissivity keeps it a single exactly stored float.
    if (self%split_rates) then
      do z=1, self%nz
        do y=1, self%ny
          do x=1, self%nx
            do f=self%fs, self%fe
              self%split_opac(f,x,y,z,1) = self%work(f,x,y,z,2)
              self%split_opac(f,x,y,z,2) = 0.
              self%split_id(f,x,y,z) = 0
              self%split_chi(f,x,y,z,1) = self%work(f,x,y,z,1)
              self%split_chi(f,x,y,z,2) = 0.
              self%split_chi_id(f,x,y,z) = 0
            end do
          end do
        end do
      end do
    endif

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
                                              i_set=i_set, iatom=i)
        endif
      end associate
    end do

    ! line_scattering: the LTE lines were added as pure thermal absorbers (source
    ! Bnu).  The scattering fraction (1-eps) of their opacity now has its source
    ! redirected from Bnu to the mean intensity, Snu_line = eps*Bnu + (1-eps)*Jnu,
    ! by adding Delta_eta = chi_scat*(Jnu - Bnu) to the emissivity without
    ! altering chi.  That is done in work_to_transfer_units, which already forms
    ! eta and so absorbs the whole correction into two extra flops per point,
    ! with no separate pass over the work arrays and no division by beta.
    self%do_line_scat = do_scat

  contains

    !=============================================================================
    !> Full detailed line list (non-sampled mode): every species, with departures.
    !> chi_lines/eta_lines, when present, receive the increments added to the work
    !> arrays.
    !=============================================================================
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
          call comp%add_linelist_opacities(self%sub_sp, xhat, ehat, dlam, vproj, &
            departures=departures, one_m_eps=self%one_m_eps, chi_scat=self%chi_scat, &
            beta=self%beta, i_set=i_set, window_lam=window_lam, &
            chi_lines=chi_lines, eta_lines=eta_lines)
        else
          call comp%add_linelist_opacities(self%sub_sp, xhat, ehat, dlam, vproj, &
                                           departures=departures, i_set=i_set, &
                                           window_lam=window_lam, &
                                           chi_lines=chi_lines, eta_lines=eta_lines)
        endif
      end associate
    end subroutine detailed_linelist

    !=============================================================================
    !> The detailed line list of a 1D output patch, computed once per set and then
    !> copied.  Without velocities it is isotropic, and temperature, densities and
    !> departures stay fixed through the output scattering iterations and the final
    !> pass, so every direction of every pass reuses it.  The one input that moves
    !> is the kernels' wing-cutoff reference: work xhat carries sigma*(Jnu-Bnu), so later
    !> passes would cut line wings against a slightly different continuum.  The store
    !> keeps the cutoff of the first pass, where Jnu=Bnu; the difference is below
    !> the kernel tolerance.  The first foreground of the (i_nu, set) fills the store
    !> while its i_mu twins and the final patch wait on the lock.  Adding the stored
    !> increments repeats the kernel's own last step.
    !=============================================================================
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
  !> Free the stored line-list and background opacity of this patch's output i_nu.
  !> The final pass is their last reader: the i_mu twins have finished their
  !> iterations before it runs.
  !===============================================================================
  SUBROUTINE free_point_tables (self)
    class(rt_patch_t) :: self
    integer :: i
    !-----------------------------------------------------------------------------
    if (associated(self%m3d_parent%bg_point_tables)) then
      call self%m3d_parent%bg_point_locks(self%i_nu)%set('bg_point_table')
      associate(store => self%m3d_parent%bg_point_tables(self%i_nu))
        if (store%initialised) then
          call mem_unregister(store%chi, label='bg_point_table')
          call mem_unregister(store%sigma, label='bg_point_table')
          deallocate(store%chi, store%sigma)
          store%initialised = .false.
        endif
      end associate
      call self%m3d_parent%bg_point_locks(self%i_nu)%unset('bg_point_table')
    endif
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
  END SUBROUTINE free_point_tables

  !===============================================================================
  !> Formal solution of the current direction with the solver of the pass
  !> (rt_output_params or rt_iteration_params solver).  The iteration pass admits
  !> only the solvers that build Psi.
  !===============================================================================
  SUBROUTINE radiative_transfer(self, comp, vec, do_rotate_atmos, i_mu)
    class(rt_patch_t) :: self
    class(composition_type) :: comp        !< composition of the slab
    real(8), intent(in) :: vec(3)          !< direction of the ray
    logical, intent(in) :: do_rotate_atmos !< the composition is ray-aligned
    integer, intent(in) :: i_mu            !< direction
    class(composition_type), pointer :: comp_up
    !-----------------------------------------------------------------------------
    ! composition of the real upstream neighbour patch, one plane before self's
    ! ghost boundary; long_char_hermite uses it in place of a same-patch proxy
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
        call self%long_char_linear(vec)
        case('quadratic')
        call self%long_char_quadratic(vec)
        case('hermite')
        ! the chain rule needs the composition gradients along the ray; in 3D
        ! that means ray-aligned (rotated) compositions, since the rotate_array
        ! shortcut above rotates only chi/Snu, not the thermodynamic state
        if (.not. (one_dimensional .or. column_by_column .or. do_rotate_atmos)) &
          error stop "rt_output_params solver='hermite' in 3D requires rotate_atmos=.true."
        call self%long_char_hermite(comp, vec, comp_up)
        case('hermite2')
        if (.not. (one_dimensional .or. column_by_column .or. do_rotate_atmos)) &
          error stop "rt_output_params solver='hermite2' in 3D requires rotate_atmos=.true."
        call self%long_char_hermite2(comp, vec, comp_up)
        case default
        error stop 'Unknown rt_output_params solver, use "linear", "quadratic", "hermite" or "hermite2"'
      end select
    else
      ! The iteration pass exists to build rates or the scattering ALI, and both
      ! need the approximate operator.  Only the solvers that fill Psi are
      ! admitted here; the quadratic and long-characteristic Hermite schemes do
      ! not, and with them preconditioning would silently degrade to a plain
      ! lambda iteration.  They remain available for the output pass.
      if (one_dimensional .or. column_by_column) then
        select case(rt_iter%solver)
          case('linear')
            call self%long_char_linear(vec)
          case default
            error stop 'rt_iteration_params solver: only "linear" builds Psi in 1D'
        end select
      else
        select case(rt_iter%solver)
          case('linear')
            call self%short_char(vec)
          case('hermite')
            call self%short_char_hermite(comp, vec)
          case default
            error stop 'rt_iteration_params solver: only "linear" and "hermite" build Psi in 3D'
        end select
      endif
    endif
  END SUBROUTINE radiative_transfer

  !===============================================================================
  !> Adds the current direction to the angle integrals: Jnu (and PsiJ) when there
  !> is an iteration, and the bound-bound rates of the atoms whose rates the slab
  !> accumulates.
  !===============================================================================
  SUBROUTINE aniso_integrations(self, comp, i_mu, i_set, do_rotate_atmos, wmu, vproj)
    class(rt_patch_t) :: self
    class(composition_type), pointer :: comp       !< composition, ray-aligned with do_rotate_atmos
    integer, intent(in) :: i_mu                    !< direction
    integer, intent(in) :: i_set                   !< set
    logical, intent(in) :: do_rotate_atmos         !< the composition is ray-aligned
    real(4), intent(in) :: wmu                     !< quadrature weight of the direction
    real(4), dimension(:,:,:), intent(in) :: vproj !< velocity along the ray [c]
    class(atom_data_t), pointer :: adata
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
          call self%bound_bound_rates(comp, atom, trans, adata, wmu, vproj, i)
        endif
      end associate
    end do
  END SUBROUTINE aniso_integrations

  !===============================================================================
  !> After the last direction of a set: the bound-free rates, which take the
  !> angle-averaged Jnu, and the transfer of the slab's rates to the atoms
  !> (integrate_rates).
  !===============================================================================
  SUBROUTINE iso_integrations(self)
    class(rt_patch_t) :: self
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
  END SUBROUTINE iso_integrations

  !===============================================================================
  !> Does this patch iterate as one of an i_mu hemisphere pair?  Such twins exist
  !> for do_full_sp when rt_iter%maxiter>0, and for the output-spectrum scattering
  !> loop when rt_out%maxiter>0 (see the create_sub_patches calls in
  !> composition_patch_mod.f90).  The do_final ray reads the output pair, but is
  !> not part of it.
  !===============================================================================
  LOGICAL FUNCTION has_twin (self)
    class(rt_patch_t) :: self
    !-----------------------------------------------------------------------------
    has_twin = .not. self%do_final .and. &
      ((self%do_full_sp .and. rt_iter%maxiter > 0) .or. &
       (self%do_out_sp .and. rt_out%maxiter > 0))
  END FUNCTION has_twin

  !===============================================================================
  !> Each twin integrates Jnu (and PsiJ) over its own hemisphere, while the next
  !> iteration needs the full-sphere average.  Both twins of a pair run the same
  !> iteration concurrently, and each overwrites its hemisphere as it goes, so a
  !> twin cannot read its partner's arrays during an execute: it would see the
  !> previous iteration or the current one, or a partial accumulation, depending
  !> on the schedule.  Instead, the twin that completes an iteration last forms
  !> the average and hands it to both.  This is safe without a lock on the data:
  !> the partner has finished the iteration, and its column cannot start the
  !> next one before this slab's column has finished the current one (the twin
  !> columns need each other, see rt_column_t%nbor_relations).
  !> Every iteration therefore reads the previous iteration's full-sphere values
  !> only, whatever the thread count or schedule.
  !>
  !> With rt_to_file each twin saves its hemisphere before announcing it, and
  !> the last one averages the two files and writes the result back to both.
  !> Jali is per patch and needs no exchange; it is saved alongside.
  !===============================================================================
  SUBROUTINE combine_twin_hemispheres (self)
    class(rt_patch_t) :: self
    class(rt_patch_t), pointer :: twin
    real(4), dimension(:,:,:,:,:), allocatable :: other
    character(8) :: Jnu_label
    logical :: last, with_PsiJ
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%combine_twin_hemispheres', itimer=itimer)
    twin => self%dn_twin
    if (.not. associated(twin)) then
      call shared%abort ('rt_patch_t%combine_twin_hemispheres: no twin')
    endif
    Jnu_label = merge('Jnu_out', 'Jnu    ', self%do_out_sp)
    ! add2meta=.false.: PsiJ and Jali are scratch, deleted again in finish(), and
    ! the readers pass mem= (never alloc=), so nothing consults the meta file.
    if (rt_to_file) then
      call self%save5(self%Jnu, trim(Jnu_label), key=self%key(1:3), add2meta=.false.)
      if (self%do_scat_precond) then
        call self%save5(self%PsiJ, 'PsiJ', key=self%key(1:3), add2meta=.false.)
        call self%save5(self%Jali, 'Jali', key=self%key(1:3), add2meta=.false.)
      endif
    endif

    !$omp critical (twin_hemisphere_combine)
    self%PsiJ_combined = .false.
    self%hemisphere_iter = self%iter
    last = twin%hemisphere_iter == self%iter
    !$omp end critical (twin_hemisphere_combine)
    if (last) then
      with_PsiJ = self%do_scat_precond .and. twin%do_scat_precond
      if (rt_to_file) then
        allocate(other, mold=self%Jnu)
        call twin%load5(trim(Jnu_label), key=twin%key(1:3), mem=other)
        call average(self%Jnu, other)
        call self%save5(self%Jnu, trim(Jnu_label), key=self%key(1:3), add2meta=save_Jnu)
        call twin%save5(self%Jnu, trim(Jnu_label), key=twin%key(1:3), add2meta=save_Jnu)
        if (with_PsiJ) then
          call twin%load5('PsiJ', key=twin%key(1:3), mem=other)
          call average(self%PsiJ, other)
          call self%save5(self%PsiJ, 'PsiJ', key=self%key(1:3), add2meta=.false.)
          call twin%save5(self%PsiJ, 'PsiJ', key=twin%key(1:3), add2meta=.false.)
        endif
        deallocate(other)
      else
        call average(self%Jnu, twin%Jnu)
        twin%Jnu = self%Jnu
        if (with_PsiJ) then
          call average(self%PsiJ, twin%PsiJ)
          twin%PsiJ = self%PsiJ
        endif
      endif
      self%combined_iter = self%iter
      twin%combined_iter = self%iter
      self%PsiJ_combined = with_PsiJ
      twin%PsiJ_combined = with_PsiJ
    endif
    call trace%end(itimer)
  contains
    subroutine average (a, b)
      real(4), intent(inout) :: a(:,:,:,:,:)
      real(4), intent(in) :: b(:,:,:,:,:)
      integer :: f, x, y, z, s
      do s=1, size(a,5)
        do z=1, size(a,4)
          do y=1, size(a,3)
            do x=1, size(a,2)
              do f=1, size(a,1)
                a(f,x,y,z,s) = 0.5*(a(f,x,y,z,s) + b(f,x,y,z,s))
              end do
            end do
          end do
        end do
      end do
    end subroutine average
  END SUBROUTINE combine_twin_hemispheres

  !===============================================================================
  !> Before a twin's next iteration, check that the pair's combine has happened
  !> (see combine_twin_hemispheres) and, with rt_to_file, load the combined
  !> state for all sets at once.  Each set then reads its previous value before
  !> overwriting it, whether it lives in memory or came from the file.
  !===============================================================================
  SUBROUTINE restore_combined_state (self)
    class(rt_patch_t) :: self
    character(8) :: Jnu_label
    character(len=256) :: msg
    !-----------------------------------------------------------------------------
    if (self%combined_iter /= self%iter-1) then
      write(msg,'(A,I8,A,I8,A,2I6)') 'twin hemispheres not combined: id', self%id, &
        ' twin', self%dn_twin%id, ' iter, combined_iter', self%iter, self%combined_iter
      call console_error (msg)
      call shared%abort ('rt_patch_t%restore_combined_state: twin exchange out of step')
    endif
    if (rt_to_file) then
      Jnu_label = merge('Jnu_out', 'Jnu    ', self%do_out_sp)
      call self%load5(trim(Jnu_label), key=self%key(1:3), mem=self%Jnu)
      if (self%do_scat_precond) then
        if (self%PsiJ_combined) call self%load5('PsiJ', key=self%key(1:3), mem=self%PsiJ)
        if (any(self%Jali_valid)) call self%load5('Jali', key=self%key(1:3), mem=self%Jali)
      endif
    endif
  END SUBROUTINE restore_combined_state

  !===============================================================================
  !> Read set i_set of this patch's saved Jnu (label Jnu or Jnu_out), written by
  !> save5 as one stream of shape (nf, nx, ny, nz, n_set).  The sets are the
  !> slowest index, so one set is a contiguous block at a fixed offset.
  !===============================================================================
  SUBROUTINE load_Jnu_set (self, label, i_set, mem)
    class(rt_patch_t) :: self
    character(len=*), intent(in) :: label !< file label, Jnu or Jnu_out
    integer, intent(in) :: i_set          !< set to read
    real(4), intent(out) :: mem(:,:,:,:)  !< Jnu of the set
    character(len=256) :: filename
    character(len=4096) :: folder
    integer :: lunit
    integer(8) :: pos
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%load_Jnu_set', itimer=itimer)
    call self%get_meta_name(filename, folder, label, self%key(1:3))
    pos = 1_8 + int(i_set-1,8)*size(mem,kind=8)*int(storage_size(mem)/8,8)
    call self%save_lock%set_read('load')
    open(newunit=lunit, file=trim(folder)//filename, access='stream', status='old', action='read')
    read(lunit, pos=pos) mem
    close(lunit)
    call self%save_lock%unset_read('load')
    call trace%end(itimer)
  END SUBROUTINE load_Jnu_set

  !===============================================================================
  !> Saves chi, Snu and Inu of the main set, and Inu of the continuum set, for
  !> direction i_mu as save_chi, save_Snu and save_Inu ask.
  !===============================================================================
  SUBROUTINE save_debug_arrays(self, i_set, i_mu)
    class(rt_patch_t) :: self
    integer, intent(in) :: i_set !< set
    integer, intent(in) :: i_mu  !< direction
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
  !> The LTE continuum extinction chi_cont and its scattering part sigma, with
  !> beta and Bnu.  They come from the sampled table (do_sample_opac), from the
  !> composition directly or, in a 1D atmosphere, from the store that the first
  !> slab of the window fills.  On the make_eos path the first execution saves the
  !> table instead.  The sampled line-list table is filled here on first use.
  !===============================================================================
  SUBROUTINE lte_background_opacities (self, comp, include_absmet, do_sample_opac, do_make_eos)
    class(rt_patch_t) :: self
    class(composition_type) :: comp !< composition of the slab
    logical :: include_absmet       !< include the absmet background
    logical :: do_sample_opac       !< interpolate the sampled table
    logical :: do_make_eos          !< build the opacity table and save it (make_eos)
    real(4) :: h_kT, hv_kT
    type(spectrum_type) :: tab_sp
    type(rt_params_t), pointer :: rtp
    type(bg_opac_t), pointer :: store
    type(lock_t), pointer :: lock
    logical :: use_store
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
            ! every use of the Planck function in transfer units reads this one
            ! float, so that J = Bnu (initial_Jnu, fixed_Jnu) leaves exactly zero
            ! in J - Bnu and in the deviation integrals of bound_bound_rates
            self%Bnu(f,x,y,z) = self%beta(f,x,y,z) / (1 - self%beta(f,x,y,z))
          end do
        end do
      end do
    end do

    ! Separate LTE line-list opacity table (non-departure species only), keyed by
    ! i_nu like opac and shared across i_mu twins. Built velocity-free here; the
    ! Doppler shift is applied when it is interpolated into the foreground.  It
    ! precedes the 1D background store, whose early return would skip it.
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

    ! The background of a 1D atmosphere depends on i_nu alone (it takes no Doppler
    ! shift), so the first rt_patch of an i_nu stores chi_cont and sigma while its
    ! twins wait on the lock; they and every later iteration and pass copy them back.
    use_store = associated(self%m3d_parent%bg_point_tables) .and. .not. do_make_eos
    if (use_store) then
      lock  => self%m3d_parent%bg_point_locks (self%i_nu)
      store => self%m3d_parent%bg_point_tables(self%i_nu)
      call lock%set('bg_point_table')
      if (store%initialised) then
        self%chi_cont = store%chi
        self%sigma    = store%sigma
        call lock%unset('bg_point_table')
        call trace%end(itimer)
        return
      endif
    endif

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
            ! The background table is continuum-only, so it needs only the nodes of
            ! the background continua.  With sample_linelist the line opacity lives
            ! in the separate ll_opac table (built below), not baked into the
            ! background, so it can be interpolated with a per-cell Doppler shift
            ! into the foreground rather than the velocity-blind background.  The
            ! absmet opacity is added after the interpolation, from the osmet table.
            call comp%get_table_sp(self%sub_sp, tab_sp)
            call self%generate_opacity_table(tab_sp, include_absmet=.false., &
              include_llist=.false., build_deriv=.true., disjoint=.true.)
          endif
        endif
        if (associated(self%opac_lock)) call self%opac_lock%unset('opac_table')
      endif
    endif

    ! getting background LTE absorption
    if (do_make_eos .and. self%iter==0) then
      call self%save_opac_table()
    elseif (do_sample_opac) then
      call self%opac%interpolate(comp%temp, comp%ne, comp%chi500, &
        self%chi_cont, self%sigma, self%sp%lam, self%fs, self%fe)
      if (include_absmet) call comp%add_absmet_extinction(self%sub_sp, self%chi_cont, self%sigma)
    else
      call comp%opacity(self%sub_sp, self%chi_cont, self%sigma, &
        include_absmet=include_absmet, include_llist=.false., &
        beta=self%beta, sorted=self%sp%sorted)
    endif

    if (use_store) then
      allocate(store%chi, source=self%chi_cont)
      allocate(store%sigma, source=self%sigma)
      call mem_register(store%chi, label='bg_point_table')
      call mem_register(store%sigma, label='bg_point_table')
      store%initialised = .true.
      call lock%unset('bg_point_table')
    endif
    call trace%end(itimer)
  END SUBROUTINE lte_background_opacities

  !===============================================================================
  !> The background emissivity eta = chi_cont*Bnu + sigma*(Jnu - Bnu) of set
  !> i_set, leaving Jnu - Bnu in Jnu.  Jnu is the previous iteration's
  !> full-sphere mean intensity (accelerated by the scattering ALI where it is
  !> active) or, in the first iteration, the one initial_Jnu chooses.
  !===============================================================================
  SUBROUTINE lte_background_emission (self, comp, i_set, vec, do_rotate_atmos)
    class(rt_patch_t), target :: self
    class(composition_type) :: comp        !< composition of the slab
    integer, intent(in) :: i_set           !< set
    real(8), intent(in) :: vec(3)          !< direction of the ray, read with do_rotate_atmos
    logical, intent(in) :: do_rotate_atmos !< the composition is ray-aligned
    real(4), dimension(:,:,:,:), allocatable :: Jali_old
    real(8) :: op, Jfs, Jali_new, op_raw, op_max, amp_max
    integer :: n_neg
    real(4) :: wscat
    character(8) :: Jnu_label
    logical :: Jnu_to_Bnu, do_scat_now
    integer :: fs, fe
    integer, save :: itimer
    integer :: x, y, z, f
    character(len=256) :: msg
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%lte_background_emission', itimer=itimer)

    fs = self%fs
    fe = self%fe

    select case (uppercase(initial_Jnu))
    case ('B', 'BNU', 'PLANCK')
      Jnu_to_Bnu = .true.
    case ('0', 'ZERO')
      Jnu_to_Bnu = .false.
    case default
      call console_warn (ch_m3d, 'unknown initial_Jnu, use "Planck" or "zero"; using "Planck"')
      Jnu_to_Bnu = .true.
    end select

    if (self%iter > 0 .or. (self%do_final .and. rt_out%maxiter > 0)) then
      Jnu_to_Bnu = .false.
      !-------------------------------------------------------------------------
      ! A twin already holds the previous iteration's full-sphere Jnu in its own
      ! array (combine_twin_hemispheres, restore_combined_state).  The do_final
      ! ray reads the output hemisphere pair, which has finished iterating:
      ! the last combine wrote the same full-sphere Jnu to both halves, so one
      ! twin's copy of this ray's set is all it needs.  With rt_to_file that is
      ! a 1/n_set slice of one file, read for each final ray and set.
      !-------------------------------------------------------------------------
      if (self%do_final) then
        if (rt_to_file) then
          Jnu_label = merge('Jnu_out', 'Jnu    ', self%do_out_sp)
          call self%up_twin%load_Jnu_set(trim(Jnu_label), i_set, self%Jnu(:,:,:,:,i_set))
        else
          self%Jnu(:,:,:,:,i_set) = self%up_twin%Jnu(:,:,:,:,i_set)
        endif
      endif
      !-------------------------------------------------------------------------
      ! Jnu enters eta via comp%scattering, so a corrupt Jnu (from the twins'
      ! exchange, whichever way it is passed) shows up as an infinite source
      ! function.  Separate that from a corrupt boundary (see
      ! set_upstream_boundary) and from a corrupt opacity.
      !-------------------------------------------------------------------------
      if (check_NaN) then
        do z=1, self%nz
          do y=1, self%ny
            do x=1, self%nx
              do f=fs, fe
                if (.not. ieee_is_finite(self%Jnu(f,x,y,z,i_set))) then
                  write(msg,'(A,I8,A,2(I8),A,4(I4))') &
                    'RACE: bad Jnu in', self%id, ' twins', self%up_twin%id, self%dn_twin%id, &
                    ' at f,x,y,z=', f, x, y, z
                  call console_error (msg)
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
      ! J_FS, the average of the twins' hemispheres.  The source it feeds is
      !   eta = chi_tot*Bnu + sigma_scat*(J - Bnu),
      ! so J appears on both sides and plain Lambda-iteration converges only at
      ! the scattering albedo per iteration.  Replace J_FS with
      !   Jali = (J_FS - op*Jali_old)/(1 - op),   op = sigma_scat*PsiJ_full,
      ! the standard (1 - Lambda* a)^-1 acceleration.  Psi is stored as
      ! Lambda_diag/chi (rt_solvers_mod), so the dimensionless self-coupling is
      ! sigma_scat*PsiJ: multiply by the scattering *opacity*, mirroring the
      ! bound-free `phi = PsiJ * cont%alpha`.  PsiJ_full is the twins' average,
      ! because each twin integrates only its own hemisphere; unlike bound-free,
      ! which consumes PsiJ linearly, this update divides by (1-op) and so must
      ! see the full sphere.  combine_twin_hemispheres leaves it in self%PsiJ.
      !
      ! op <= 1 identically: sigma is the scattering part of chi_cont
      ! by construction on both the exact and the sampled path, and
      ! Lambda_diag <= 1.  So scat_op_max is not there to rescue an out-of-range
      ! ratio; it bounds the amplification 1/(1-op) where op legitimately
      ! approaches 1 (see its declaration in m3d_patch_mod).  The difference
      ! in the numerator cancels heavily as op -> 1, so it is formed in real(8).
      ! Where there is no scattering op = 0 and Jali = J_FS exactly.
      !-------------------------------------------------------------------------
      ! The full-sphere PsiJ exists only where both twins accumulated one in the
      ! previous iteration.  do_scat_precond follows any(get_rates), which drops
      ! as soon as the atom converges, so the pair may not have (see
      ! PsiJ_combined).  Skip the correction for this iteration then.  This is a
      ! local decision only: self%do_scat_precond stays up so PsiJ keeps
      ! accumulating and the recurrence can resume in the next iteration.
      do_scat_now = self%do_scat_precond .and. self%PsiJ_combined
      if (do_scat_now) then
        allocate(Jali_old(fs:fe,self%nx,self%ny,self%nz))
        if (self%Jali_valid(i_set)) then
          Jali_old = self%Jali(:,:,:,:,i_set)
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
                op = self%sigma(f,x,y,z) * dble(self%PsiJ(f,x,y,z,i_set))
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
    ! term.  Because sigma is built together with chi_cont (summed after
    ! interpolation on the sampled path, accumulated side by side on the exact
    ! one), the effective single-scattering albedo sigma/chi_cont is <= 1 by
    ! construction, which is what makes this plain Lambda-iteration a
    ! contraction rather than an amplifier.
    !
    ! Jnu is left holding (Jnu - Bnu) unconditionally.  The foreground's
    ! line-scattering term needs the same quantity, but Jnu is zeroed for a
    ! fresh angle accumulation before it runs, so keep a copy in dJnu.
    !
    ! diag_no_scat_loop severs the scattering Lambda-loop while holding every
    ! optical depth fixed, by zeroing the weight rather than branching inside
    ! the loop (see m3d_debug_mod).
    !---------------------------------------------------------------------------
    wscat = 1.0
    if (diag_no_scat_loop) wscat = 0.0
    do z=1, self%nz
      do y=1, self%ny
        do x=1, self%nx
          do f=fs, fe
            associate(Jnu=>self%Jnu(f,x,y,z,i_set), chi=>self%chi_cont(f,x,y,z), &
                      eta=>self%ehat(f,x,y,z), Bnu=>self%Bnu(f,x,y,z), &
                      sigma=>self%sigma(f,x,y,z))
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
  !> Fills the opacity table opac on the nodes of tab_comp, normalised by chi500.
  !> Both components are always filled.  disjoint (the sampling path) keeps
  !> val = true absorption and scat = scattering, each normalised, so that
  !> sample_interpolate can sum them after interpolation and the sampled albedo
  !> stays <= 1.  Without it (the make_eos path) val is the total and scat stays
  !> in 1/cm, which is what save_opac_table writes out.
  !===============================================================================
  SUBROUTINE generate_opacity_table (self, sp, include_absmet, include_llist, build_deriv, &
                                     disjoint)
    class(rt_patch_t), target :: self
    type(spectrum_type) :: sp                 !< frequencies of the table
    logical, intent(in) :: include_absmet     !< include the absmet background
    logical, intent(in) :: include_llist      !< include the line list
    logical, intent(in) :: build_deriv        !< prepare the PCA axes and node derivatives for interpolation
    logical, intent(in), optional :: disjoint !< keep absorption and scattering apart
    class(composition_type), pointer :: comp
    real(4), dimension(:,:,:,:), pointer :: all_chi, all_scat
    logical :: do_disjoint
    real(4) :: over_chi500
    integer :: j,k, ntemp, nne, nnu
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
  !> Build the separate sampled LTE line-list table (sample_linelist). tab_comp is
  !> masked to the non-departure species (see mask_departure_species), so the line
  !> opacity added here is pure LTE. In LTE the two working forms coincide
  !> (xhat_line = ehat_line = chi_line/(1-beta)); the table stores that single
  !> quantity, normalised by chi500, and interpolation adds it to both work slots.
  !> The pure line delta is isolated by subtracting the continuum working xhat that
  !> add_linelist needs in place for its opacity-threshold reference.
  !===============================================================================
  SUBROUTINE generate_ll_opacity_table (self, sp, include_absmet)
    class(rt_patch_t), target :: self
    type(spectrum_type) :: sp             !< frequencies of the window
    logical, intent(in) :: include_absmet !< include the absmet background in the reference continuum
    class(composition_type), pointer :: comp
    type(spectrum_type) :: psp
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
    ! hands back the working form directly (continuum_opacities returns
    ! xhat = ehat = chi/(1-beta)), so there is nothing to convert here; only the
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
      call comp%add_linelist_opacities(psp, all_chi, all_eta, dlam, &
        one_m_eps=one_m_eps, chi_scat=chi_scat, beta=all_beta)
      allocate(self%ll_opac%scat(nnu, ntemp, nne), source=0.)
      call mem_register(self%ll_opac%scat, label='ll_opac_table')
    else
      call comp%add_linelist_opacities(psp, all_chi, all_eta, dlam)
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
  !> Writes the make_eos opacity table per unit mass [cm^2/g]: the absorption to
  !> 'opa' and the scattering to 'scat'.
  !===============================================================================
  SUBROUTINE save_opac_table (self)
    class(rt_patch_t), target :: self
    class(composition_type), pointer :: comp
    real(4) :: over_rho
    integer :: j,k, ntemp, nrho
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
  !> Zeroes the cells 1:nx, 1:ny of array over its z and last axes.
  !===============================================================================
  SUBROUTINE reset_array (self, array)
    class(rt_patch_t), target :: self
    real(4), dimension(:,:,:,:,:), intent(inout) :: array !< array per (f, x, y, z, k)
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
  !> Adds the bound-free extinction and emissivity of the atom's continua to xhat
  !> and ehat, with the populations of set i_set.
  !===============================================================================
  SUBROUTINE add_bound_free_opacities (self, comp, atom, trans, adata, i_set)
    class(rt_patch_t), target :: self
    class(composition_type) :: comp     !< composition of the slab
    class(atom_type) :: atom            !< model atom
    class(atom_trans_t) :: trans        !< the atom's transitions in the window
    class(atom_data_t), target :: adata !< the atom's data on the mesh
    integer, intent(in) :: i_set        !< set
    class(cont_type), pointer :: cont
    real(4) :: ni, tmp, alpha
    real(8) :: dev
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
        dev = atom%ev(cont%j) - atom%ev(cont%i)
        do z=1, self%nz
          do y=1, self%ny
            do x=1, self%nx
              ni  = adata%nstar(x,y,z,cont%i) * adata%b(x,y,z,cont%i) * pop_ratio(x,y,z)
              tmp = adata%nstar(x,y,z,cont%i) * adata%b(x,y,z,cont%j) * pop_ratio(x,y,z)
              ! ehat holds eta/beta, which for a cut continuum (see radiative_cut) is
              ! the overflowing n_j n*_i/n*_j.  Its source function is below
              ! exp(-cut_x) b_j/b_i of 2h nu^3/c^2, so the continuum only absorbs.
              if (radiative_cut(dev, comp%e_kT(x,y,z))) tmp = 0
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
  !> The Voigt damping parameter of each line other than hydrogen, in units of
  !> the Doppler width (trans%adamp); hprofl_mod handles the hydrogen lines.
  !===============================================================================
  SUBROUTINE bound_bound_damping (self, comp, atom, trans, adata)
    class(rt_patch_t), target :: self
    class(composition_type) :: comp !< composition of the slab
    class(atom_type) :: atom        !< model atom
    class(atom_trans_t) :: trans    !< the atom's transitions in the window
    class(atom_data_t) :: adata     !< the atom's data on the mesh
    class(line_type), pointer :: line
    real(4), dimension(:,:,:), allocatable :: tothi, vbar2, ne, temp
    real(4), dimension(:,:,:), allocatable :: stark_pe, unsold
    integer :: i
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin (trim(self%kind)//'%bound_bound_damping', itimer=itimer)
    allocate(temp, ne, tothi, vbar2, stark_pe, unsold, mold=adata%inv_df)
    ne = real(comp%ne)
    temp = real(comp%temp)
    tothi = real(comp%tothi) ! total neutral hydrogen density
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
  !> Adds the line opacity and emissivity of the atom in the current direction to
  !> the work arrays and, with split_rates, to the tracked contributions.  In the
  !> first iteration of the internal spectrum it also sums each line's profile
  !> integral over the directions (adata%wphi).
  !===============================================================================
  SUBROUTINE add_bound_bound_opacities (self, comp, atom, trans, adata, wmu, vproj, LTE, &
                                        i_set, iatom)
    class(rt_patch_t), target :: self
    class(composition_type) :: comp                !< composition of the slab
    class(atom_type) :: atom                       !< model atom
    class(atom_trans_t) :: trans                   !< the atom's transitions in the window
    class(atom_data_t), target :: adata            !< the atom's data on the mesh
    real(4) :: wmu                                 !< quadrature weight of the direction
    real(4), contiguous, dimension(:,:,:) :: vproj !< velocity along the ray [c]
    logical :: LTE                                 !< take b = 1 (the LTE companion set)
    integer, intent(in) :: i_set                   !< set
    integer, intent(in) :: iatom                   !< index of the atom
    class(line_type), pointer :: line
    type(hline_t) :: hline
    type(hcell_t) :: hcell
    real(4), dimension(:,:,:), allocatable :: wphi, temp, ne, fo, vbar2h, tothi, t3nhe
    real(4), dimension(:), allocatable :: vphi
    real(4) :: idf, vpr, a, ni
    real(8) :: q
    real(4) :: lower_off, ekt, wj, chi_fac, eta_fac, bi, bj, Rij
    integer :: i, fs, fe
    logical :: do_LTE, get_wphi, sum_wphi
    real(8), allocatable :: pop_ratio(:,:,:)
    integer, save :: itimer=0
    integer :: x, y, z, f, index, k
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
      if (line%hydrogen) hline = hydrogen_line(line%i, line%j, line%nu0)
      fs = line%iblue
      fe = line%ired
      sum_wphi = get_wphi .and. sums_over_windows(self, line)

      ! A line whose lower level lies lower_off above the atom's level i (a component
      ! of a combined level) absorbs from that level's share exp(-lower_off/kT).  Both
      ! terms of the line opacity take that share: the stimulated emission factor
      ! beta(nu) already holds the line's own photon energy, so the line source
      ! function is the Planck function in LTE, whatever the upper level's offset.
      lower_off = line%lower_ev - atom%ev(line%i)

      ! Adding bound-bound opacities
      do z=1, self%nz
        do y=1, self%ny
          do x=1, self%nx
            idf = adata%inv_df(x,y,z)
            ekt = comp%e_kT(x,y,z)
            ni = adata%nstar(x,y,z,line%i) * pop_ratio(x,y,z)

            if (do_LTE) then
              bi = 1.0
              bj = 1.0
            else
              bi = adata%b(x,y,z,line%i)
              bj = adata%b(x,y,z,line%j)
            endif

            if (lower_off /= 0) then
              bi = bi * exp(-lower_off*ekt)
              bj = bj * exp(-lower_off*ekt)
            endif

            associate(phi => vphi(fs:fe))
              if (line%hydrogen) then
                k = self%hvcs(iatom)%slot(i)
                if (k > 0) hcell%vcs = self%hvcs(iatom)%p(x,y,z,k)
                call hcell_prepare(hline, temp(x,y,z), ne(x,y,z), fo(x,y,z), vbar2h(x,y,z), &
                                   tothi(x,y,z), t3nhe(x,y,z), idf, hcell)
                if (k > 0) self%hvcs(iatom)%p(x,y,z,k) = hcell%vcs
                do f=fs, fe
                  q = line%q(f) - vproj(x,y,z)
                  vphi(f) = hprofile(hline, hcell, q)
                end do
                if (H_occupation_probability) then
                  wj = comp%elements(1)%ions(1)%levels(line%j)%w(x,y,z)
                  bi = bi * wj
                  bj = bj * wj
                endif
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
              if (self%split_rates) then
                index = self%split_offset(iatom) + i
                do f=fs, fe
                  call add_largest(vphi(f)*eta_fac, index, &
                    self%split_opac(f,x,y,z,1), self%split_opac(f,x,y,z,2), &
                    self%split_id(f,x,y,z))
                  call add_largest(vphi(f)*chi_fac, index, &
                    self%split_chi(f,x,y,z,1), self%split_chi(f,x,y,z,2), &
                    self%split_chi_id(f,x,y,z))
                end do
              endif

              if (sum_wphi) then
                wphi(x,y,z) = sum(phi * line%wq(fs:fe))
              endif
            end associate
          end do
        end do
      end do

      if (sum_wphi) then
        call lock_mem(self%m3d_parent, 'wphi')
        associate(total => adata%wphi(line%ntrans))
          if (.not. allocated(total%w)) &
            call os%allocate(total%w, shape(wphi), label='wphi', zero=.true.)
          total%w = total%w + wmu * wphi
        end associate
        call unlock_mem(self%m3d_parent, 'wphi')
      endif
    end do
    call trace%end(itimer)
  END SUBROUTINE add_bound_bound_opacities

  !===============================================================================
  !> Does the profile integral of this line have to be summed over wavelength
  !> windows?  It does when a window edge cuts the line (it then reaches the edge
  !> of this window).  Hydrogen lines keep the summed integral as well, since
  !> their rates are not normalised per direction.
  !===============================================================================
  LOGICAL FUNCTION sums_over_windows (self, line)
    class(rt_patch_t), intent(in) :: self
    class(line_type), intent(in) :: line !< line of the window
    !-----------------------------------------------------------------------------
    sums_over_windows = line%iblue <= self%fs .or. line%ired >= self%fe &
                        .or. line%hydrogen
  END FUNCTION sums_over_windows

  !===============================================================================
  !> Adds the current direction's share of the atom's bound-bound rates to
  !> trans%gamma_ij and gamma_ji, which hold a thermal part and the deviation from
  !> it per line and cell (see below); line_rate_coefficients forms the rates.
  !===============================================================================
  SUBROUTINE bound_bound_rates (self, comp, atom, trans, adata, wmu, vproj, iatom)
    class(rt_patch_t), target :: self
    class(composition_type) :: comp                !< composition of the slab
    class(atom_type) :: atom                       !< model atom
    class(atom_trans_t) :: trans                   !< the atom's transitions in the window
    class(atom_data_t), target :: adata            !< the atom's data on the mesh
    real(4) :: wmu                                 !< quadrature weight of the direction
    real(4), contiguous, dimension(:,:,:) :: vproj !< velocity along the ray [c]
    integer, intent(in) :: iatom                   !< index of the atom
    class(line_type), pointer :: line
    type(hline_t) :: hline
    type(hcell_t) :: hcell
    real(4), dimension(:,:,:), allocatable :: temp, ne, fo, vbar2h, tothi, t3nhe
    real(4), dimension(:), allocatable :: wmunu, vph
    real(4) :: idf, vpr, a, ni, nj, phi
    real(4) :: vvv, sum_I, sum_P, sum_A, sum_T, sum_U, nb, inv_nref
    real(4) :: g_th, g_dev
    real(8) :: beta_ref
    real(8) :: q
    real(4) :: lower_off, ekt, eta_coef, g_ratio
    real(4) :: split_bi, split_bj, split_Rij, abs_coef
    real(4) :: psi, eta_u, eta_t, chi_t, eta_other, other, chi_other, up, spont
    real(8), allocatable :: pop_ratio(:,:,:)
    integer :: i, fs, fe, index
    logical :: split_rates, per_direction, own_norm, do_scat
    real(4), dimension(:), allocatable :: wqf
    real(4) :: sum_w, wtot
    integer, save :: itimer=0
    integer :: x, y, z, f, k
    !-----------------------------------------------------------------------------
    call trace%begin (trim(self%kind)//'%bound_bound_rates', itimer=itimer)
    ! Each line stores two integrals per cell instead of the rates themselves:
    !
    !   th  = A int phi * spont * Bnu/n_ref          (thermal part)
    !   dev = A int phi * (up - Bnu*spont)           (deviation from it)
    !
    ! with up and spont the (preconditioned) absorption and spontaneous integrands,
    ! Bnu = beta/(1-beta) the Planck occupation number at each frequency, and
    ! n_ref = beta_ref/(1-beta_ref) the Planck occupation number at the level
    ! energy difference, beta_ref = n*_j g_i/(n*_i g_j).  line_rate_coefficients
    ! (atom_mod) forms the rates from them in real(8):
    !
    !   R_ij = n_ref*th + dev,   R_ji = g_i/g_j * ((n_ref+1)*th + dev)
    !
    ! The thermal parts are then in detailed balance to real(8) rounding, n*_i R_ij
    ! = n*_j R_ji for J = Bnu, whatever the profile, the frequency quadrature, the
    ! line's own nu0 and the real(4) rounding of th.  (n_ref+1)*Bnu/n_ref equals
    ! (Bnu+1)*beta/beta_ref, so the downward thermal part carries the
    ! exp(-h(nu-nu0)/kT) of the emissivity add_bound_bound_opacities uses; the
    ! deviation part is taken with the absorption profile in both directions.
    ! For J = Bnu without preconditioning (fixed_Jnu), dev is exactly zero.
    split_rates = self%split_rates
    do_scat = self%do_line_scat
    ! A fixed-Jnu spectrum takes the rates unnormalised (set_fixed_Jnu)
    per_direction = .not. (self%do_full_sp .and. fixed_Jnu > 0)
    wtot = sum(rt_iter%wmu)
    allocate(wmunu(self%fs:self%fe), vph(self%fs:self%fe), wqf(self%fs:self%fe))
    ! The preconditioning correction has to subtract the same local term the formal
    ! solver added, so it is built from the emissivity coefficient that
    ! add_bound_bound_opacities used, pop_ratio included.
    allocate(pop_ratio, mold=adata%totn)
    call comp%population_ratio(atom%id, self%i_set_main, pop_ratio)
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
      if (line%hydrogen) hline = hydrogen_line(line%i, line%j, line%nu0)
      fs = line%iblue
      fe = line%ired
      wmunu(fs:fe) = wmu * line%wnu(fs:fe) * line%Aij
      ! A line inside one window is normalised by its own profile integral in this
      ! direction, which differs from the direction-summed wphi only by the
      ! variation of the quadrature error between directions.
      own_norm = per_direction .and. .not. sums_over_windows(self, line)
      if (own_norm) wqf(fs:fe) = line%wq(fs:fe)
      ! the lower-level share of add_bound_bound_opacities, for both levels
      lower_off = line%lower_ev - atom%ev(line%i)

      do z=1,self%nz
        do y=1,self%ny
          do x=1,self%nx
            ni = adata%b(x,y,z,line%i) * adata%nstar(x,y,z,line%i)
            nj = adata%b(x,y,z,line%j) * adata%nstar(x,y,z,line%j)
            idf = adata%inv_df(x,y,z)
            vpr = vproj(x,y,z)
            a = trans%adamp(x,y,z,i)
            ekt = comp%e_kT(x,y,z)
            ! 1/n_ref, held below the real(4) range: it only reaches the bound in
            ! cells far beyond radiative_cut, whose rates the matrix does not use
            inv_nref = 0.0
            if (adata%nstar(x,y,z,line%i) > 0 .and. adata%nstar(x,y,z,line%j) > 0) then
              beta_ref = adata%nstar(x,y,z,line%j) * atom%g(line%i) &
                       / (adata%nstar(x,y,z,line%i) * atom%g(line%j))
              inv_nref = real(min((1d0 - beta_ref)/beta_ref, 1d35), 4)
            endif
            if (lower_off /= 0) then
              ni = ni * exp(-lower_off*ekt)
              nj = nj * exp(-lower_off*ekt)
            endif
            if (preconditioning) then
              ! Same arithmetic, in the same order, as add_bound_bound_opacities:
              ! the split tracking can only take this line's emissivity back out
              ! of the local source term if it is the identical float.
              split_bj = adata%b(x,y,z,line%j)
              if (lower_off /= 0) split_bj = split_bj * exp(-lower_off*ekt)
              if (line%hydrogen .and. H_occupation_probability) &
                split_bj = split_bj * comp%elements(1)%ions(1)%levels(line%j)%w(x,y,z)
              split_Rij = adata%nstar(x,y,z,line%i) * pop_ratio(x,y,z)
              split_Rij = split_Rij * line%Bij
              eta_coef = split_Rij * split_bj
              split_bi = adata%b(x,y,z,line%i)
              if (lower_off /= 0) split_bi = split_bi * exp(-lower_off*ekt)
              if (line%hydrogen .and. H_occupation_probability) &
                split_bi = split_bi * comp%elements(1)%ions(1)%levels(line%j)%w(x,y,z)
              abs_coef = split_Rij * split_bi
              ! line_rate_coefficients scales the downward integral by g_i/g_j, so this
              ! ratio is what makes the two corrections cancel in the net rate at the
              ! old populations.
              g_ratio = 0.0
              if (nj > 0.0) g_ratio = (ni/nj) * (atom%g(line%j)/atom%g(line%i))
            endif

            sum_A = 0.0
            sum_I = 0.0
            sum_P = 0.0
            sum_T = 0.0
            sum_U = 0.0
            sum_w = 0.0
            ! One f-loop per rate formulation and profile kind, so that none
            ! carries the rate_formulation branch.
            if (line%hydrogen) then
              k = self%hvcs(iatom)%slot(i)
              if (k > 0) hcell%vcs = self%hvcs(iatom)%p(x,y,z,k)
              call hcell_prepare(hline, temp(x,y,z), ne(x,y,z), fo(x,y,z), vbar2h(x,y,z), &
                                 tothi(x,y,z), t3nhe(x,y,z), idf, hcell)
              if (k > 0) self%hvcs(iatom)%p(x,y,z,k) = hcell%vcs
              if (split_rates) then
                do f=fs, fe
                  q = line%q(f) - vpr
                  vph(f) = hprofile(hline, hcell, q)
                end do
              else
                do f=fs, fe
                  associate(Inu=>self%work(f,x,y,z,3), Psi=>self%Psi(f,x,y,z))
                    q = line%q(f) - vpr
                    phi = hprofile(hline, hcell, q)
                    vvv = phi * wmunu(f)
                    sum_A = sum_A + vvv * (self%Bnu(f,x,y,z)*inv_nref)
                    sum_I = sum_I + vvv * Inu
                    sum_P = sum_P + vvv * phi * Psi * self%beta(f,x,y,z)
                  end associate
                end do
              endif
            else
              if (split_rates) then
                vph(fs:fe) = qlg_bot_profile(a, line%q(fs:fe), vpr, idf)
              else
                do f=fs, fe
                  associate(Inu=>self%work(f,x,y,z,3), Psi=>self%Psi(f,x,y,z))
                    phi = qlg_bot_profile(a, line%q(f), vpr, idf)
                    vvv = phi * wmunu(f)
                    if (own_norm) sum_w = sum_w + phi * wqf(f)
                    sum_A = sum_A + vvv * (self%Bnu(f,x,y,z)*inv_nref)
                    sum_I = sum_I + vvv * Inu
                    sum_P = sum_P + vvv * phi * Psi * self%beta(f,x,y,z)
                  end associate
                end do
              endif
            endif

            ! The split rate integrands, from the profile the branches above
            ! stored: Rybicki & Hummer (1992), Eqs. (2.24) and (2.25), which
            ! precondition within the same transition.  Irem is the formal
            ! solution without its own local term, so adding back the other
            ! contributions' share of it never forms Inu - Psi*eta_t as a
            ! difference of two nearly equal numbers.  The spontaneous factor
            ! 1 - Psi*chi_t is built the same way, as (chi_other + escape*chi_t)/chi
            ! with escape = 1 - Psi*chi from the solver: in the core of a very
            ! thick line both 1 - Psi*chi and 1 - chi_t/chi are near the fp32
            ! rounding of one, and so is their true value, 1/dtau.  Where Psi is
            ! zero (fixed_Jnu) the unpreconditioned integrands are taken.  The
            ! tracked contributions carry no beta, and the line-scattering
            ! emissivity is not among them: it emits into the source term but
            ! adds no extinction of its own.
            if (split_rates .and. preconditioning) then
              index = self%split_offset(iatom) + i
              do f=fs, fe
                psi = self%Psi(f,x,y,z)
                eta_u = vph(f)*eta_coef
                eta_t = eta_u*self%beta(f,x,y,z)
                chi_t = vph(f)*abs_coef - eta_t
                eta_other = self%beta(f,x,y,z)*other_contribution(eta_u, index, &
                  self%split_opac(f,x,y,z,1), self%split_opac(f,x,y,z,2), self%split_id(f,x,y,z))
                chi_other = other_contribution(vph(f)*abs_coef, index, self%split_chi(f,x,y,z,1), &
                  self%split_chi(f,x,y,z,2), self%split_chi_id(f,x,y,z)) - eta_other
                other = eta_other
                if (do_scat) other = other + self%chi_scat(f,x,y,z)*self%dJnu(f,x,y,z)
                up = merge(self%work(f,x,y,z,3), self%Irem(f,x,y,z) + psi*other, psi == 0.)
                spont = merge(1.0, (chi_other + self%escape(f,x,y,z)*chi_t) &
                                   / self%work(f,x,y,z,1), psi == 0.)
                vvv = vph(f) * wmunu(f)
                nb = self%Bnu(f,x,y,z)
                sum_T = sum_T + vvv*spont*(nb*inv_nref)
                sum_U = sum_U + vvv*(up - nb*spont)
              end do
            else if (split_rates) then
              do f=fs, fe
                vvv = vph(f) * wmunu(f)
                nb = self%Bnu(f,x,y,z)
                sum_T = sum_T + vvv*(nb*inv_nref)
                sum_U = sum_U + vvv*(self%work(f,x,y,z,3) - nb)
              end do
            endif

            if (split_rates) then
              g_th = sum_T
              g_dev = sum_U
            else
              ! The 'original' formulation keeps its own sums, with the
              ! preconditioning corrections that cancel in the net rate at the old
              ! populations; its upward and downward integrals, sum_I - sum_P*eta_coef
              ! and sum_I + sum_A - sum_P*eta_coef*g_ratio, are split into th and dev
              ! after the fact, so its rates keep the real(4) rounding of a
              ! difference.
              g_th = sum_A
              g_dev = sum_I
              if (preconditioning) then
                g_th = g_th - sum_P * eta_coef * (g_ratio - 1.0)
                g_dev = g_dev - sum_P * eta_coef
              endif
              if (inv_nref > 0.0) g_dev = g_dev - g_th/inv_nref
            endif
            if (own_norm) then
              if (split_rates) sum_w = sum(vph(fs:fe) * wqf(fs:fe))
              ! lines below the strength cutoff have no profile weight and no rate
              if (sum_w > 0.0) then
                g_th = g_th / (sum_w * wtot)
                g_dev = g_dev / (sum_w * wtot)
              endif
            endif

            trans%gamma_ij(x,y,z,i) = trans%gamma_ij(x,y,z,i) + g_th
            trans%gamma_ji(x,y,z,i) = trans%gamma_ji(x,y,z,i) + g_dev
          end do
        end do
      end do
    end do
    call trace%end(itimer)
  END SUBROUTINE bound_bound_rates

  !===============================================================================
  !> Keep the largest-magnitude term of a running sum apart from the rest, so
  !> that one term can later be removed from the total without a cancelling
  !> subtraction.  The remainder never contains the largest term, even when a
  !> later term replaces it.
  !===============================================================================
  pure subroutine add_largest (value, index, largest, rest, largest_index)
    real(4), intent(in) :: value            !< the term
    integer, intent(in) :: index            !< its contribution index
    real(4), intent(inout) :: largest       !< largest term so far
    real(4), intent(inout) :: rest          !< sum of the other terms
    integer, intent(inout) :: largest_index !< contribution index of the largest term
    !-----------------------------------------------------------------------------
    ! This and other_contribution live here, next to their only callers, rather
    ! than in a module of their own: they sit in the innermost frequency loops of
    ! add_bound_bound_opacities and bound_bound_rates, and the build has no
    ! link-time optimisation, so across a module boundary they would be a call
    ! per frequency point per line per cell.
    if (abs(value) > abs(largest)) then
      rest = rest + largest
      largest = value
      largest_index = index
    else
      rest = rest + value
    endif
  end subroutine add_largest

  !===============================================================================
  !> The total less one named contribution, formed by addition only.
  !===============================================================================
  pure function other_contribution (value, index, largest, rest, largest_index) result(other)
    real(4), intent(in) :: value         !< the contribution to remove
    integer, intent(in) :: index         !< its contribution index
    real(4), intent(in) :: largest       !< largest term
    real(4), intent(in) :: rest          !< sum of the other terms
    integer, intent(in) :: largest_index !< contribution index of the largest term
    real(4) :: other
    !-----------------------------------------------------------------------------
    if (index == largest_index) then
      other = rest
    else
      other = largest + (rest-value)
    endif
  end function other_contribution

  !===============================================================================
  !> Adds the atom's bound-free rates from the angle-averaged Jnu of the main set,
  !> preconditioned with PsiJ when preconditioning is set.
  !===============================================================================
  SUBROUTINE bound_free_rates (self, atom, trans, adata)
    class(rt_patch_t), target :: self
    class(atom_type) :: atom            !< model atom
    class(atom_trans_t) :: trans        !< the atom's transitions in the window
    class(atom_data_t), target :: adata !< the atom's data on the mesh
    class(cont_type), pointer :: cont
    real(4), dimension(:), allocatable :: w_alpha, wmunu
    real(4) :: phi
    real(4) :: g_ij, g_ji, vij, vji, sum_abs, sum_ems, sum_pre
    real(8) :: dev
    integer :: i, j, fs, fe
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin (trim(self%kind)//'%bound_free_rates', itimer=itimer)
    allocate(wmunu(self%fs:self%fe))
    allocate(w_alpha(self%fs:self%fe))
    ! Computing bound-free radiative rates based on half integrated mean intensity.
    ! Cells where radiative_cut holds get no radiative rates: the emission sum would
    ! be built from clamped beta and its LTE ratio overflows real(4).  The atom's
    ! rate matrix puts a detailed-balance pair in their place.
    do j=1, size(trans%cont)
      cont => trans%cont(j)
      if (.not. cont%active) cycle
      i = j + size(trans%line)
      fs = cont%iblue
      fe = cont%ired
      wmunu(fs:fe) = 0.5 * cont%wnu(fs:fe)
      w_alpha(fs:fe) = self%twohnu3_c2(fs:fe) * cont%alpha(fs:fe) * wmunu(fs:fe)
      dev = atom%ev(cont%j) - atom%ev(cont%i)

      do z=1,self%nz
        do y=1,self%ny
          do x=1,self%nx
            if (radiative_cut(dev, self%m3d_parent%composition%e_kT(x,y,z))) cycle
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
  !> Fills the upstream ghost plane of work: no incoming radiation at the top of
  !> the domain, a plane extrapolated from the first cells at the bottom (with
  !> Inu = Snu), and the exit plane of the upstream slab (mem) everywhere else.
  !===============================================================================
  SUBROUTINE set_upstream_boundary (self, i_mu, i_set, ng)
    class(rt_patch_t) :: self
    integer :: i_mu  !< direction
    integer :: i_set !< set
    integer :: ng    !< number of ghost planes
    integer :: i
    real(4) :: xxold(3), xxnew(ng)
    integer, save :: itimer=0
    integer :: x, y, z, f
    character(len=256) :: msg
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
      ! The incoming boundary comes straight from the upstream slab, via the
      ! buffer the column handed over.  If it is already non-finite here, the
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
                  write(msg,'(A,I8,A,I8,A,3(I4),A,I3,A,2(I4),A,ES12.4)') &
                    'RACE: bad upstream boundary in', self%id, ' from', self%upstream%id, &
                    ' at f,x,y=', f, x, y, ' w=', i, ' i_mu,i_set=', i_mu, i_set, &
                    ' value=', self%work(f,x,y,z,i)
                  call console_error (msg)
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
  !> Adds the intensity of the current direction to Jnu with weight 2*wmu, since
  !> each twin covers one hemisphere, and with save_PsiJ Psi to PsiJ.
  !===============================================================================
  SUBROUTINE integrate_intensity (self, wmu, i_set, save_PsiJ)
    class(rt_patch_t) :: self
    real(4) :: wmu       !< quadrature weight of the direction
    integer :: i_set     !< set
    logical :: save_PsiJ !< accumulate PsiJ too
    real(4) :: wmu2
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
  !> Adds the slab's rates to the atom's rates on the composition patch, under its
  !> lock.
  !===============================================================================
  SUBROUTINE integrate_rates (self, atom, trans, adata)
    class(rt_patch_t) :: self
    class(atom_type) :: atom     !< model atom
    class(atom_trans_t) :: trans !< the atom's transitions in the window
    class(atom_data_t) :: adata  !< the atom's data on the mesh
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
  !> Copies the exit plane of work to mem for the next slab, allocating mem on
  !> first use, and samples the debug point for print_post_update_info.
  !===============================================================================
  SUBROUTINE set_downstream_boundary (self, i_mu, i_set)
    class(rt_patch_t) :: self
    integer :: i_mu  !< direction
    integer :: i_set !< set
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
          call console_warn (ch_m3d, 'NaN in rt_patch at pcheck')
        endif
      endif
      if (i_mu==self%mus .and. i_set==self%i_set_main) then   ! sample the spectrum set once
        memcheck(1:3) = self%work(self%fcheck, p(1), p(2), p(3),:)
        associate(beta=>self%beta(self%fcheck, p(1), p(2), p(3)))
          memcheck(2:3) = memcheck(2:3) * (1 - beta) / beta
        end associate
      endif
    endif

    if (.not. associated(self%mem)) then
      call self%allocate_boundary()
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
  !> ehat enters holding the physical emissivity eta (from lte_background_emission)
  !> and leaves holding the working form eta/beta.  xhat must therefore be formed
  !> first, while ehat is still eta: the two statements do not commute.
  !===============================================================================
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
  !> Seeds the extinction and source slots of work with the direction-independent
  !> xhat and ehat.
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
  !> Close the accumulation: turn the working pair (xhat, ehat) into the (chi, Snu)
  !> pair the formal solvers consume, and apply the line-scattering source on the
  !> way.  chi = xhat - ehat*beta is the physical extinction, scattering opacity
  !> included (the lines went in as absorbers, so chi is whole).  eta = ehat*beta +
  !> chi_scat*(Jnu - Bnu) redirects the (1-eps) fraction of the line source from Bnu
  !> to Jnu, and Snu = eta/chi.  |chi| is held at least maser_chi_floor*xhat, so
  !> that an inverted population keeps a bounded source function.
  !===============================================================================
  SUBROUTINE work_to_transfer_units (self)
    class(rt_patch_t) :: self
    logical :: do_scat, split_rates
    integer, save :: itimer=0
    integer :: x, y, z, f, i, n_floored, first(4)
    real(4) :: xhat, eta, eta_rest, chi, floor, first_val(3)
    real(4), allocatable :: chiv(:)
    character(len=256) :: msg
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%work_to_transfer_units', itimer=itimer)
    allocate(chiv(self%fs:self%fe))
    n_floored = 0
    do_scat = self%do_line_scat
    split_rates = self%split_rates
    do z=1, self%nz
      do y=1, self%ny
        do x=1, self%nx
          ! chi (to chiv) and eta (to the ehat slot) come first, one f-loop per
          ! rate formulation, so that neither carries the other's branch; the
          ! maser floor and Snu = eta/chi follow in a loop of their own.
          if (split_rates) then
            do f=self%fs, self%fe
              ! Take exactly the total the split tracking holds, so that removing
              ! one transition's contribution reproduces the solver's local term.
              ! Each part is scaled by beta on its own and the line-scattering
              ! source joins the rest, in the order the rates add them up.
              eta_rest = self%beta(f,x,y,z)*self%split_opac(f,x,y,z,2)
              eta = self%beta(f,x,y,z)*self%split_opac(f,x,y,z,1) + eta_rest
              chiv(f) = self%work(f,x,y,z,1) - eta
              if (do_scat) eta = self%beta(f,x,y,z)*self%split_opac(f,x,y,z,1) &
                               + (eta_rest + self%chi_scat(f,x,y,z)*self%dJnu(f,x,y,z))
              self%work(f,x,y,z,2) = eta
            end do
          else
            do f=self%fs, self%fe
              eta = self%work(f,x,y,z,2) * self%beta(f,x,y,z)
              chiv(f) = self%work(f,x,y,z,1) - eta
              if (do_scat) eta = eta + self%chi_scat(f,x,y,z)*self%dJnu(f,x,y,z)
              self%work(f,x,y,z,2) = eta
            end do
          endif
          do f=self%fs, self%fe
            xhat = self%work(f,x,y,z,1)
            eta  = self%work(f,x,y,z,2)
            chi  = chiv(f)
            ! Each line adds phi*Rij*(bi - bj*beta) to chi, which goes negative for
            ! an inverted population.  chi may then pass through zero, and eta/chi
            ! diverges.  Keep the sign (the inversion is real) but hold |chi| away
            ! from zero, so that Snu stays bounded by 1/maser_chi_floor.
            floor = max(real(maser_chi_floor) * abs(xhat), tiny(1.0))
            if (suppress_masers) then
              chi = max(chi, floor)
            elseif (abs(chi) < floor) then
              chi = sign(floor, chi)
              if (n_floored == 0) then
                first = [f, x, y, z]
                first_val = [chi, xhat, eta]
              endif
              n_floored = n_floored + 1
            endif
            self%work(f,x,y,z,1) = chi
            self%work(f,x,y,z,2) = eta / chi
          end do
        end do
      end do
    end do
    ! one warning per execution, for the first 20 that floor chi
    if (check_NaN .and. n_floored > 0) then
      !$omp atomic capture
      n_negative_chi = n_negative_chi + 1
      i = n_negative_chi
      !$omp end atomic
      if (i <= 20) then
        write(msg,'(A,A,I0,A,4(I5),4(1x,ES13.5))') trim(self%m3d_patch_info()), &
          '  CHI floored in ', n_floored, ' points, first f,x,y,z=', first, first_val, &
          self%beta(first(1),first(2),first(3),first(4))
        call console_warn (ch_m3d, msg)
      endif
    endif
    call trace%end(itimer)
  END SUBROUTINE work_to_transfer_units

  !===============================================================================
  !> Shifts each z-plane of array horizontally (periodically, with bilinear
  !> interpolation) so that the rays of direction vec become vertical, with the
  !> plane at pivot_z left in place.
  !===============================================================================
  SUBROUTINE rotate_array (self, vec, array)
    class(rt_patch_t) :: self
    real(8), intent(in) :: vec(3)                         !< direction of the ray
    real(4), intent(inout), dimension(:,:,:,:,:) :: array !< array per (f, x, y, z, k)
    real(4), dimension(2) :: delta, deltak
    real(4), dimension(:,:,:), allocatable :: shifted
    real(4) :: bottom, bot2pivot, height, f1, f2, a, b, c, d
    integer :: up1, up2, x1, x2, y1, y2, i
    integer :: x, y, z
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin('rt_patch%rotate_array', itimer=itimer)
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
  !> Sets the intensity of the current direction to fixed_Jnu*Bnu in place of a
  !> formal solution, and the profile integrals wphi of the atoms to 1.
  !===============================================================================
  SUBROUTINE set_fixed_Jnu (self)
    class(rt_patch_t) :: self
    integer :: i, j
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%set_fixed_Jnu', itimer=itimer)
    do z=1, self%nz
      do y=1, self%ny
        do x=1, self%nx
          do f=self%fs, self%fe
            associate(Bnu=>self%Bnu(f,x,y,z), Inu=>self%work(f,x,y,z,3))
              Inu = fixed_Jnu * Bnu
            end associate
          end do
        end do
      end do
    end do
    do i=1, size(self%atoms, 1)
      associate(adata => self%atoms(i)%data)
        if (allocated(adata%wphi)) then
          do j=1, size(adata%wphi)
            if (allocated(adata%wphi(j)%w)) adata%wphi(j)%w = 1
          end do
        endif
      end associate
    end do
    call trace%end(itimer)
  END SUBROUTINE set_fixed_Jnu

  !===============================================================================
  !> The velocity of each cell along direction vec.
  !===============================================================================
  SUBROUTINE get_vproj (self, vproj, vec, rotate, i_mu)
    class(rt_patch_t) :: self
    real(4), dimension(:,:,:) :: vproj !< velocity along the ray [c]
    real(8), dimension(:) :: vec       !< direction of the ray
    logical :: rotate                  !< take the ray-aligned composition of direction i_mu
    integer :: i_mu                    !< direction
    class(composition_type), pointer :: comp
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
  !> The contribution function of the emergent intensity per direction: the
  !> continuum set stores its chi, eta and Inu, and the main set combines them
  !> with its own into the contribution to the line depression (depression_cntrbf)
  !> or the emissivity.  The main set's last direction saves 'cntrbf'.
  !===============================================================================
  SUBROUTINE save_contribution_function (self, i_set, i_mu)
    class(rt_patch_t) :: self
    integer :: i_set !< set
    integer :: i_mu  !< direction
    real(4) :: net_chi, net_src, net_eta, lin_chi, lin_eta
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    ! The continuum set comes first, so it stashes and the main set combines.
    ! The stash keeps all three continuum quantities (chi, eta, Inu) rather than
    ! pre-combining them, so that the main set forms the differences chi-chi_cnt
    ! and eta-eta_cnt directly: deep in the atmosphere those are catastrophically
    ! cancelling, and any algebraically equivalent regrouping that multiplies out
    ! first loses most of the significant digits there.
    !
    ! The stash has a direction axis because execute loops over the sets outside
    ! the directions (unless rotate_atmos), and with velocities chi and eta
    ! depend on the direction.
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
  !> Prints the debug sample of set_downstream_boundary at console level 2.
  !===============================================================================
  SUBROUTINE print_post_update_info (self)
    class(rt_patch_t) :: self
    !-----------------------------------------------------------------------------
    if (console_on(ch_m3d, 2, all_ranks=.true.)) then
      call self%print_to_stdout(arr=memcheck(1:3))
    endif
  END SUBROUTINE print_post_update_info

  !===============================================================================
  !> Slab states are swept by their rt_column_t and never enter the task list.
  !> They still get a link of their own, since patch_t%init takes the level and
  !> hash key (the file key) from self%link%task.
  !===============================================================================
  SUBROUTINE append_to_task_list (self)
    class(rt_patch_t), target :: self
    class(m3d_halo_t), pointer :: link
    !-----------------------------------------------------------------------------
    allocate (link)
    call link%init()
    link%task => self
    self%link => link
  END SUBROUTINE append_to_task_list

  !===============================================================================
  !> Finishes the slab.  Once both twins have finished, frees their Jnu or, with
  !> rt_to_file, removes its files (unless save_Jnu) and those of the scattering
  !> ALI.
  !===============================================================================
  SUBROUTINE finish(self)
    class(rt_patch_t) :: self
    character(8) :: Jnu_label
    logical :: cleanup
    ! ----------------------------------------------------------------------------
    call self%m3d_patch_t%finish()
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
        if (rt_to_file) then
          if (.not. save_Jnu) then
            Jnu_label = merge('Jnu_out', 'Jnu    ', self%do_out_sp)
            call self%up_twin%remove(Jnu_label, key=self%up_twin%key(1:3))
            call self%dn_twin%remove(Jnu_label, key=self%dn_twin%key(1:3))
          endif
          !-----------------------------------------------------------------
          ! The scattering ALI persists PsiJ and Jali alongside Jnu (see
          ! combine_twin_hemispheres), but unlike Jnu they are pure iteration
          ! scratch (there is no save_PsiJ), so they go regardless of save_Jnu.
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
  !> patch_t%init allocates mem here.  A slab's mem is only its boundary buffer,
  !> which allocate_boundary creates when a sweep first needs it.
  !===============================================================================
  SUBROUTINE allocate_mem(self, m_in)
    class(rt_patch_t) :: self
    integer, intent(in), optional :: m_in(6) !< unused; required by patch_t%allocate_mem
  END SUBROUTINE allocate_mem

  !===============================================================================
  !> Keep the exit plane of a frozen slab's last sweep (rt_patch_pre_update) for
  !> the active slab after it: in mem, or with rt_to_file in a file, freeing mem.
  !===============================================================================
  SUBROUTINE retain_exit_plane (self)
    class(rt_patch_t) :: self
    !-----------------------------------------------------------------------------
    if (.not. rt_to_file) return
    call self%save6(self%mem, 'exit_plane', key=self%key(1:3), add2meta=.false.)
    call self%deallocate_mem()
  END SUBROUTINE retain_exit_plane

  !===============================================================================
  !> Give the next slab of the column, nxt, a copy of the retained exit plane.
  !===============================================================================
  SUBROUTINE restore_exit_plane (self, nxt)
    class(rt_patch_t) :: self
    class(rt_patch_t) :: nxt !< the next slab of the column
    !-----------------------------------------------------------------------------
    call nxt%allocate_boundary()
    if (rt_to_file) then
      call self%load6('exit_plane', key=self%key(1:3), mem=nxt%mem)
    else
      nxt%mem = self%mem
    endif
  END SUBROUTINE restore_exit_plane

  !===============================================================================
  !> Free the retained exit plane once no slab reads it.  Repeated calls do nothing.
  !===============================================================================
  SUBROUTINE release_exit_plane (self)
    class(rt_patch_t) :: self
    !-----------------------------------------------------------------------------
    if (rt_to_file) then
      call self%remove('exit_plane', key=self%key(1:3))
    elseif (associated(self%mem)) then
      call self%deallocate_mem()
    endif
  END SUBROUTINE release_exit_plane

  !===============================================================================
  !> The boundary buffer: the exit plane (chi, S, I) of every ray and set.
  !===============================================================================
  SUBROUTINE allocate_boundary(self)
    class(rt_patch_t) :: self
    real(4), dimension(:,:,:,:,:,:), pointer, contiguous :: tmp
    integer :: m(6)
    ! ----------------------------------------------------------------------------
    m = [self%nf, self%n(1), self%n(2), self%n_mu*self%n_set, 1, 3]
    call self%patch_t%allocate_mem(m)
    tmp => self%mem
    self%mem(self%fs:self%fe,1:self%n(1),1:self%n(2),self%mus:self%mue,1:self%n_set,1:3) => tmp
  END SUBROUTINE allocate_boundary

  !===============================================================================
  !> Saves the emergent intensity of the final pass from mem: 'ie', averaged over
  !> the surface [erg/s/cm^2/A/sr], 'i3' per surface cell with save_resolved, and
  !> the flux 'flux' [erg/s/cm^2/A].
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
    ! nphi_copies), matching tumult/flux.py's integrate_flux.
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
  !> Tracks the relative change between iterations of this window's emergent flux
  !> at the top of the domain, and folds its largest value into the global value
  !> of the round through the RMA counters that the population convergence check
  !> of nlte_patch_t%execute uses as well.  It runs in every iteration, whichever
  !> criterion decides (flux_convergence); nlte_patch_t prints the report (worst
  !> change, wavelength, i_nu) from the shared counters.
  !>
  !> rt_iter%wmu are the isotropic weights of the Jnu angle sum, normalised to 1
  !> over the full sphere, so the flux takes the extra factor mu = rt_iter%vec(3,i):
  !> F = 4*pi * sum_i rt_iter%wmu(i) * rt_iter%vec(3,i) * I_i over the upward
  !> directions mus:mue.  The normalisation cancels in the relative change.
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

      ! expected = number of top-of-domain window reporters per round (bb and bf
      ! windows with decouple_continuum).  round_pos is this report's 1-based
      ! position within its round, from the monotonic report counter: round_pos==1
      ! starts a fresh round (resets the running max to this report), and the
      ! decision is evaluated only at round_pos==expected, once every reporter has
      ! folded its value into the round's running max.  The population check
      ! instead resets when a reporter (a z column) is seen again, which is safe
      ! with many reporters but with as few as two (n_windows=1 with
      ! decouple_continuum) would evaluate convergence after a single report.
      expected = merge(2*n_windows, n_windows, decouple_continuum)

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
  !> Registers the memory of arr with os.
  !===============================================================================
  SUBROUTINE mem_register(arr, label)
    real(4), dimension(..) :: arr       !< the array
    character(len=*), optional :: label !< label of the allocation
    !-----------------------------------------------------------------------------
    call os%mem_register(storage_size(arr), product(shape(arr)), label)
  END SUBROUTINE mem_register

  !===============================================================================
  !> Unregisters the memory of arr with os.
  !===============================================================================
  SUBROUTINE mem_unregister(arr, label)
    real(4), dimension(..) :: arr       !< the array
    character(len=*), optional :: label !< label of the allocation
    !-----------------------------------------------------------------------------
    call os%mem_register(storage_size(arr), -product(shape(arr)), label)
  END SUBROUTINE mem_unregister

  !===============================================================================
  !> Verify that the output patches of a segmented (window_file) grid tile it
  !> exactly: in i_nu order each patch is non-empty, starts where the previous
  !> one ended (the windows are adjacent in grid index) and lies inside one window.
  !===============================================================================
  SUBROUTINE verify_output_patch_coverage(output_sp, offset)
    class(atomic_spectrum_type), intent(in) :: output_sp !< output spectrum
    integer, intent(in) :: offset(:)                     !< patches before each window, and the total last
    integer :: i_nu, i_seg, fs, fe, prev_fe, n_out
    character(len=256) :: msg
    !---------------------------------------------------------------------------
    if (output_sp%n_segments == 0) return

    if (size(offset) /= output_sp%n_segments + 1) then
      write(msg,'(A,I0,A,I0)') &
        'output patch offsets for ', size(offset)-1, ' windows on a grid of ', output_sp%n_segments
      call console_error (msg)
      error stop 'output patch offsets and n_segments mismatch'
    endif

    n_out = offset(size(offset))
    prev_fe = 0
    do i_nu = -1, -n_out, -1
      call output_patch_range(output_sp, offset, i_nu, fs, fe)
      i_seg = output_patch_window(offset, i_nu)
      if (fs /= prev_fe + 1 .or. fe < fs .or. &
          fs < output_sp%segment_first(i_seg) .or. fe > output_sp%segment_last(i_seg)) then
        write(msg,'(A,I0,A,I0,A,I0,A,I0)') 'output patch i_nu=', i_nu, &
          ' covers [', fs, ',', fe, '] in window ', i_seg
        call console_error (msg)
        error stop 'output patches do not tile the windows'
      endif
      prev_fe = fe
    enddo
    if (prev_fe /= output_sp%nnu) then
      write(msg,'(A,I0,A,I0)') 'output patches end at ', prev_fe, ' of ', output_sp%nnu
      call console_error (msg)
      error stop 'incomplete segment coverage'
    endif

    if (console_on(ch_m3d, 1)) write(console_unit,'(A,I0,A,I0,A)') 'Verified: ', n_out, &
      ' output patches tile ', output_sp%n_segments, ' windows'
  END SUBROUTINE verify_output_patch_coverage

END MODULE rt_patch_mod
