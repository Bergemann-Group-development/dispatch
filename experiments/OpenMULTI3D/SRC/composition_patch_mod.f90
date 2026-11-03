#include "../HOP/read_namelist.inc"
MODULE composition_patch_mod
  USE io_mod
  USE namelist_check_mod, only: namelist_error
  USE os_mod
  USE omp_mod
  USE task_mod, only: task_t
  USE link_mod, only: nbor_t
  USE gpatch_mod
  USE handyfuncs_mod
  USE optab_mod, only: optab, pca_opac_grid
  USE osmet_mod, only: osmet, eosmet
  USE ion_mod, only: particle_type, ion_type, debye_shield
  USE molecule_mod, only: molec_type, irwin_dis
  USE element_mod, only: element_type
  USE m3d_patch_mod,   only: m3d_patch_t, composition, dlam, &
  use_rho, use_ne, overlapping_bb, overlapping_bf, binned_opac, initial_Jnu, &
  rotate_atmos, decouple_continuum, n_windows, n_windows_out, out_patch_offset, use_vel, &
  make_eos, overlapping_bb_in_output, atmos_format, atmos_file, &
  n_nlte_patches, nlte_finished
  USE rt_params_mod, only: rt_iter, rt_out, pivot_z
  USE m3d_debug_mod
  USE read_atmos_mod, only: ia, interpolate_atmos, pchip_interpolation
  USE interpol_mod, only: pchip_slopes
  USE composition_mod, only: composition_type, exclude_elements, only_elements, &
  ref_lam, debug_opac, abundances, relative_abundances, mhd_eos, helium_scattering, &
  departure_pointer_t, read_departures, match_nlte_levels, &
  sample_linelist, opac_table_vmic, &
  line_scattering
  USE atom_mod, only: level_partf
  USE nlte_patch_mod,  only: nlte_patch_t, output_grid_and_list
  USE rt_patch_mod, only: rt_patch_t, rt_patch_ptr_t, output_set_da, output_set_element, verify_output_patch_coverage
  USE bftrans_mod, only: H_occupation_probability
  USE bbtrans_mod, only: tolerance, disable_rad, disable_stk, disable_vdw, &
                         molecule_damping_estimates
  USE linelist_mod, only: ll_profile, voigt_kernel
  USE op_math_mod, only: interval_t, merge_intervals
  USE qlg, only: qlg_init
  USE hprofl_mod, only: hydrogen_BPO
  USE hvcs_mod, only: hydrogen_profiles, vcs_data_dir
  USE anstee_mod, only: anstee_file
  USE trace_mod
  USE console_mod
  USE scaling_mod
  USE units_mod
  USE ray_mod, only: set_carlson_quadratures, one_dimensional
  USE bits_mod
  USE atom_mod, only: atom_type, atomic_spectrum_type
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  
  !-----------------------------------------------------------------------------
  type, public, extends(m3d_patch_t) :: composition_patch_t
    class(nlte_patch_t), pointer :: nlte_patch => null()
    !> the RT slab states of this depth slab, which init_rt_columns groups into
    !> rt_column_t tasks
    type(rt_patch_ptr_t), allocatable :: rt_slabs(:)
    !> the spectra and atoms the RT slabs are made from, kept until create_rt_slabs;
    !> the same object as nlte_patch when that is linked
    class(nlte_patch_t), pointer :: spectra => null()
  contains
    procedure :: pre_init => composition_pre_init
    procedure :: execute
    procedure :: slab_groups
    procedure :: create_rt_slabs
    procedure :: create_sub_patches
    procedure :: print_post_update_info
    procedure :: save_abundances_to_file
    procedure :: init_tabulated_comp
    procedure :: estimate_eos_memory
    procedure :: init_eos_comp
    procedure :: save_eos_variables
    procedure :: save_number_densities
    procedure :: nbor_relations => composition_patch_nbor_relations
  end type

  !> The RT slabs of one spectrum: windows nu(1):nu(2) times directions mu(1):mu(2)
  type, public :: slab_group_t
    integer :: nu(2), mu(2)
    class(atomic_spectrum_type), pointer :: sp => null()
    type(atom_type), pointer :: atoms(:) => null()
  end type slab_group_t
  
  real(8), public :: FeH=-99.9, ldne=0.25, ldrho=0.25, ldtemp=0.05
  logical, public :: save_atmos=.true., save_eos=.false.
  logical, public :: pca_opac=.true.       ! .false. -> axis-aligned (lnT,lnNe) sampling grid
  real(4), public, dimension(:), allocatable :: external_abund
  real(4), public :: alpha_enhancement=0.0, r_process_enhancement=0.0
  character(len=256), save :: save_species=''
  character(len=4096), save :: absdat_file='./input/absdat/TS_absdat.txt'
  character(len=4096), save :: abund_file='./input/abund/abund_magg22.txt'
  character(len=4096), save :: isotope_file='./input/isotopes/isotopes_blp.txt'
  character(len=4096), save :: molec_file='./input/molec/molec_data_BC.txt'
  character(len=4096), save :: opac_table=''
  character(len=4096), save :: absmet_file=''
  character(len=4096), save :: linelist_folder=''
  character(len=4096), save, dimension(:), allocatable :: linelist_file

  public :: read_composition_params, max_n_windows, set_n_windows_out

contains

!===============================================================================
SUBROUTINE read_composition_params ()
  integer :: i
  character(len=256), allocatable :: names(:)
  real(4) :: tmolim=-99
  integer, save:: itimer=0, verbose=console_inherit
  namelist /composition_params/ verbose, absdat_file, abund_file, &
  isotope_file, molec_file, H_occupation_probability, &
  mhd_eos, debye_shield, &
  tmolim, irwin_dis, save_species, &
  ref_lam, debug_opac, abundances, relative_abundances, alpha_enhancement, r_process_enhancement, &
  disable_rad, disable_stk, disable_vdw, molecule_damping_estimates, helium_scattering, &
  hydrogen_BPO, hydrogen_profiles, vcs_data_dir, eosmet
  namelist /opac_table_params/ opac_table, opac_table_vmic, pca_opac, &
  ldne, ldrho, ldtemp, absmet_file
  namelist /linelist_params/ dlam, exclude_elements, linelist_file, linelist_folder, &
  tolerance, only_elements, line_scattering, ll_profile, voigt_kernel, anstee_file, &
  match_nlte_levels
  ! ----------------------------------------------------------------------------
  call trace%begin('composition%init', itimer=itimer)

  allocate(linelist_file(128))
  linelist_file = ''

  READ_NAMELIST(composition_params)
  write (io%namelist, composition_params)

  READ_NAMELIST(opac_table_params)
  write (io%namelist, opac_table_params)

  READ_NAMELIST(linelist_params)

  ! Build QLG-bot's lookup tables. Must precede any opacity accumulation -- see the
  ! kernel note in HOP/linelist_mod.f90.
  call qlg_init

  if (linelist_folder /= '') then
    call os%get_filenames(linelist_folder, names)
    deallocate(linelist_file)
    allocate(linelist_file(size(names)))
    do i=1, size(names)
      linelist_file(i) = trim(linelist_folder)//'/'//trim(names(i))
    end do
  else
    call truncate_array(linelist_file, count(linelist_file/=''))
  endif
  write (io%namelist, linelist_params)

  call console_set_level (verbose, ch_hop)

  if (absmet_file /= '') then
    call osmet%init(absmet_file)
  endif

  call composition%init(absdat_file, abund_file, FeH, &
  isotope_file=isotope_file, molec_file=molec_file, tmolim=tmolim, &
  external_abund=external_abund, alpha_enhancement=alpha_enhancement, &
  r_process_enhancement=r_process_enhancement)
  call trace%end(itimer)
END SUBROUTINE read_composition_params

!===============================================================================
SUBROUTINE composition_pre_init (self)
  class(composition_patch_t), target:: self
  class(nlte_patch_t), pointer :: nlte_patch
  class(atomic_spectrum_type), save, pointer :: output_sp
  real(8) :: ll_lam(2)
  integer :: i, nw
  integer, save:: itimer=0
  logical, save :: first_time=.true., external_sp
  ! ----------------------------------------------------------------------------
  call trace%begin('composition_patch_t%composition_pre_init', itimer=itimer)
  self%kind = 'composition'
  call self%m3d_patch_t%pre_init()
  nw = merge(1+size(rt_out%vec,2), 1, rotate_atmos)
  ! no mem: the atmosphere lives in the composition objects, the radiation in the
  ! rt patches, and the populations in atom%data
  call self%set_nv_and_nw(0,nw)
  ! this patch will interpolate from the input snapshot; execute releases it again
  call ia%acquire()

  if (first_time) then
    allocate(output_sp)
    binned_opac = opac_table /= ''
      
    if (binned_opac) then
      call optab%init(opac_table, output_sp)
    else
      call output_sp%init()
      output_sp%name = 'output'
    endif

    external_sp = output_sp%nnu > 0
    output_grid_and_list = external_sp .and. any(linelist_file /= '')
  endif

  allocate(nlte_patch)
  call nlte_patch%pre_init_atom()
  call nlte_patch%pre_init_spectrum(output_sp)

  if (associated(nlte_patch%atoms) .or. read_departures) then
    call nlte_patch%pre_init_and_link(self)
    self%nlte_patch => nlte_patch
    !$omp atomic update
    n_nlte_patches = n_nlte_patches + 1
    nlte_finished = .false. ! the output tasks now have something to wait for
  endif

  if (first_time) then
    first_time = .false.
    ! The line-list readers keep lines inside [lam(1)-dlam, lam(end)+dlam], so
    ! collect the combined wavelength range of every spectrum whose pass draws on
    ! the line list: the output spectrum when the output pass includes (or samples)
    ! it, or when make_eos bakes it into the tables; the internal iteration
    ! spectra when the iteration pass does.
    !
    ! Phase B: build intervals from spectrum segments instead of min/max range.
    ! For output spectrum with multiple segments, use those segment bounds.
    ! For non-segmented internal spectra, use full range as single interval.
    block
      type(interval_t), allocatable :: ll_intervals(:)
      integer :: n_intervals

      ll_lam = [huge(1d0), -huge(1d0)]
      call build_linelist_intervals(composition, output_sp, nlte_patch, ll_lam, dlam, &
                                    ll_intervals, n_intervals)

      if (ll_lam(1) <= ll_lam(2)) then
        do i=1, size(linelist_file)
          ! Pass the intervals to init_linelist so build_linelist selects on the
          ! union of windows, not a single min/max span.
          if (allocated(ll_intervals) .and. n_intervals > 0) then
            call composition%init_linelist(linelist_file(i), dlam, ll_lam, ll_intervals)
          else
            call composition%init_linelist(linelist_file(i), dlam, ll_lam)
          end if
        end do
      elseif (any(linelist_file /= '')) then
        ! Without a spectrum_params grid the output pass runs on the model atom's
        ! frequencies (or not at all), which carry no line-list window.
        if (rt_out%include_linelist .and. size(rt_out%vec,2) > 0 .and. .not. external_sp) then
          call console_warn (ch_hop, 'linelist_params: line list not read: the output pass has no ' // &
            'spectrum; set aa_blue, aa_red and R_sampling in spectrum_params')
        else
          call console_warn (ch_hop, 'linelist_params: line list not read: neither the output ' // &
            'nor the iteration pass includes it')
        endif
      endif
      if (allocated(ll_intervals)) deallocate(ll_intervals)
    end block
    call nlte_patch%exclude_double_lines()
  endif

  ! create_rt_slabs makes the RT slabs once n_windows is chosen
  self%spectra => nlte_patch

  if (save_species /= '') then
    call self%save_number_densities(composition, pre_init=.true.)
  endif

  call self%save_abundances_to_file()
  call self%post_pre_init()
  call trace%end(itimer)
END SUBROUTINE composition_pre_init

!===============================================================================
!> The largest useful n_windows: no window of the internal or output spectrum
!> may be empty.
!===============================================================================
INTEGER FUNCTION max_n_windows (spectra) RESULT (n)
  class(nlte_patch_t), intent(in) :: spectra
  !-----------------------------------------------------------------------------
  n = huge(1)
  if (spectra%output_sp%nnu > 0) n = min(n, spectra%output_sp%nnu)
  if (decouple_continuum) then
    if (associated(spectra%bb_sp)) then
      if (spectra%bb_sp%nnu > 0) n = min(n, spectra%bb_sp%nnu)
    endif
    if (associated(spectra%bf_sp)) then
      if (spectra%bf_sp%nnu > 0) n = min(n, spectra%bf_sp%nnu)
    endif
  elseif (associated(spectra%sp)) then
    if (spectra%sp%nnu > 0) n = min(n, spectra%sp%nnu)
  endif
END FUNCTION max_n_windows

!===============================================================================
!> Output frequency patches.  The output windows (the segments of a window_file
!> grid, or the whole grid) share n_windows patches in proportion to their number
!> of points, so the output pass has about n_windows tasks per depth for the
!> threads.  A window gets at least one patch, and only as many as leave each
!> patch min_out_patch points or more: a patch of a few points costs about as
!> much fixed overhead (opacity tables and their Doppler padding, the lines of
!> its window it visits) as a large one.
!===============================================================================
SUBROUTINE set_n_windows_out (output_sp)
  class(atomic_spectrum_type), intent(in) :: output_sp
  integer, parameter :: min_out_patch = 16
  integer, allocatable :: npts(:)
  integer :: i
  !-----------------------------------------------------------------------------
  if (allocated(out_patch_offset)) deallocate(out_patch_offset)
  if (output_sp%nnu == 0) then
    allocate(out_patch_offset(1), source=0)
    n_windows_out = 0
    return
  endif
  if (output_sp%n_segments > 0) then
    allocate(npts, source=output_sp%segment_last - output_sp%segment_first + 1)
  else
    allocate(npts(1), source=output_sp%nnu)
  endif
  allocate(out_patch_offset(size(npts)+1))
  out_patch_offset(1) = 0
  do i=1, size(npts)
    out_patch_offset(i+1) = out_patch_offset(i) + max(1, min(npts(i)/min_out_patch, &
      nint(real(n_windows,8)*npts(i)/sum(npts))))
  end do
  n_windows_out = out_patch_offset(size(npts)+1)
END SUBROUTINE set_n_windows_out

!===============================================================================
!> The RT slabs this patch makes for the current n_windows and n_windows_out:
!> both directions of every internal window when the NLTE iteration runs, both
!> directions of every output window for the output Lambda iteration, and the
!> final output ray.
!===============================================================================
FUNCTION slab_groups (self) RESULT (groups)
  class(composition_patch_t), intent(in) :: self
  type(slab_group_t), allocatable :: groups(:)
  class(nlte_patch_t), pointer :: spectra
  type(slab_group_t) :: g(4)
  integer :: n
  !-----------------------------------------------------------------------------
  spectra => self%spectra
  n = 0
  if (associated(spectra%atoms) .and. rt_iter%maxiter > 0) then
    if (decouple_continuum) then
      n = n + 1; g(n) = slab_group_t([1, n_windows], [1,2], spectra%bb_sp, spectra%bb_atoms)
      n = n + 1; g(n) = slab_group_t([1, n_windows]+n_windows, [1,2], spectra%bf_sp, spectra%bf_atoms)
    else
      n = n + 1; g(n) = slab_group_t([1, n_windows], [1,2], spectra%sp, spectra%atoms)
    endif
  endif
  ! rt_out%solver='none' disables the output pass entirely -- create no output
  ! RT patches (the final i_mu=-1 ray or the Lambda-iteration i_mu=[1,2] ones).
  if (spectra%output_sp%nnu > 0 .and. rt_out%solver /= 'none') then
    if (size(rt_out%vec,2) > 0 .or. (make_eos .and. self%key(1)==1)) then
      if (rt_out%maxiter > 0) then
        n = n + 1; g(n) = slab_group_t([-n_windows_out, -1], [1,2], spectra%output_sp, spectra%output_atoms)
      endif
      n = n + 1; g(n) = slab_group_t([-n_windows_out, -1], [-1,-1], spectra%output_sp, spectra%output_atoms)
    endif
  endif
  allocate(groups, source=g(1:n))
END FUNCTION slab_groups

!===============================================================================
!> Make the RT slabs of this depth slab, once n_windows and n_windows_out are
!> set (init_rt_columns), with the opacity tables they share.
!===============================================================================
SUBROUTINE create_rt_slabs (self)
  class(composition_patch_t) :: self
  class(nlte_patch_t), pointer :: nlte_patch
  type(slab_group_t), allocatable :: groups(:)
  integer :: i
  !-----------------------------------------------------------------------------
  nlte_patch => self%spectra
  self%n_windows = n_windows

  ! One sampled opacity table per i_nu, shared by the rt_patches created below,
  ! which differ only in i_mu and would otherwise build identical copies.  The
  ! internal spectrum uses i_nu in [1, n_windows], the bound-free half of a decoupled
  ! continuum [n_windows+1, 2*n_windows], and the output spectrum [-n_windows_out, -1]; slot 0 is
  ! never used.  make_eos is excluded: save_opac_table overwrites the table in
  ! place, so those rt_patches keep a private one.
  if ((rt_iter%sample_back_opac .or. rt_out%sample_back_opac) .and. .not. make_eos) then
    allocate(self%opac_tables(-n_windows_out:2*n_windows))
    allocate(self%opac_locks (-n_windows_out:2*n_windows))
    do i=-n_windows_out, 2*n_windows
      call self%opac_locks(i)%init('opac', i)
    end do
    ! Separate LTE line-list opacity table (sample_linelist): same i_nu keying.
    if (sample_linelist) then
      allocate(self%ll_opac_tables(-n_windows_out:2*n_windows))
      allocate(self%ll_opac_locks (-n_windows_out:2*n_windows))
      do i=-n_windows_out, 2*n_windows
        call self%ll_opac_locks(i)%init('ll_opac', i)
      end do
    endif
  endif
  ! In 1D the LTE background opacity is the same for every direction, set and
  ! iteration (it takes no Doppler shift, in 3D neither), so each i_nu of either
  ! pass computes it once (see lte_background_opacities in rt_patch_mod).
  ! make_eos writes tables instead.
  if (one_dimensional .and. .not. make_eos) then
    allocate(self%bg_point_tables(-n_windows_out:2*n_windows))
    allocate(self%bg_point_locks (-n_windows_out:2*n_windows))
    do i=-n_windows_out, 2*n_windows
      call self%bg_point_locks(i)%init('bg_point', i)
    end do
  endif
  ! In a static 1D atmosphere the detailed line opacity is isotropic, so the
  ! output pass computes it once per output i_nu and set and reuses it for every
  ! direction, scattering iteration and the final pass (see stored_linelist in
  ! rt_patch_mod); a column holds few enough points to keep it.  With a vertical
  ! velocity every direction sees its own shift, and the store stays empty.
  if (one_dimensional .and. rt_out%include_linelist .and. &
      .not. rt_out%sample_linelist) then
    allocate(self%ll_point_tables(-n_windows_out:-1))
    allocate(self%ll_point_locks (-n_windows_out:-1))
    do i=-n_windows_out, -1
      call self%ll_point_locks(i)%init('ll_point', i)
    end do
  endif

  allocate(groups, source=self%slab_groups())
  do i=1, size(groups)
    call self%create_sub_patches(groups(i)%nu, groups(i)%mu, groups(i)%sp, groups(i)%atoms)
  end do

  if (nlte_patch%output_sp%n_segments > 0) &
    call verify_output_patch_coverage(nlte_patch%output_sp, out_patch_offset)

  self%spectra => null()
  if (.not. associated(self%nlte_patch)) then
    deallocate(nlte_patch)
  endif
END SUBROUTINE create_rt_slabs

!===============================================================================
!> Phase B: Build line-list intervals from output and internal spectra.
!> When the output spectrum has multiple segments (multi-window), use those
!> segment bounds as intervals. For non-segmented internal spectra, use the
!> full range as a single interval. Expand each interval by dlam.
!===============================================================================
SUBROUTINE build_linelist_intervals(composition, output_sp, nlte_patch, ll_lam, dlam, &
                                    intervals_out, n_intervals_out)
  class(composition_type), intent(in) :: composition
  class(atomic_spectrum_type), intent(in), pointer :: output_sp
  type(nlte_patch_t), intent(in) :: nlte_patch
  real(dp), intent(out) :: ll_lam(2)
  real(dp), intent(in) :: dlam
  type(interval_t), intent(out), allocatable :: intervals_out(:)
  integer, intent(out) :: n_intervals_out
  real(dp), allocatable :: windows(:,:)
  real(dp) :: sp_min, sp_max
  integer :: i, n_windows
  logical :: output_uses_ll, internal_uses_ll
  !-----------------------------------------------------------------------------

  ! Determine which spectra use the line list
  output_uses_ll = rt_out%include_linelist .and. (size(rt_out%vec,2) > 0 .or. make_eos)
  internal_uses_ll = rt_iter%include_linelist .and. rt_iter%maxiter > 0

  ! Initialize: will be updated to encompass all used spectra
  ll_lam = [huge(1d0), -huge(1d0)]
  n_intervals_out = 0

  ! Build window array from contributing spectra
  n_windows = 0
  if (output_uses_ll .and. associated(output_sp)) then
    if (output_sp%nnu > 0) then
      if (output_sp%n_segments > 0) then
        ! Multi-segment output: use segment bounds as windows
        n_windows = output_sp%n_segments
        allocate(windows(n_windows, 2))
        do i = 1, n_windows
          windows(i, 1) = output_sp%segment_blue_vac(i)
          windows(i, 2) = output_sp%segment_red_vac(i)
        end do
      else
        ! Single-segment or legacy output: use full range
        if (allocated(windows)) deallocate(windows)
        allocate(windows(1, 2))
        windows(1, 1) = minval(output_sp%lam)
        windows(1, 2) = maxval(output_sp%lam)
        n_windows = 1
      end if
    end if
  end if

  ! Add internal spectrum ranges
  if (internal_uses_ll) then
    if (decouple_continuum) then
      if (associated(nlte_patch%bb_sp) .and. nlte_patch%bb_sp%nnu > 0) then
        if (.not. allocated(windows)) allocate(windows(1, 2))
        sp_min = minval(nlte_patch%bb_sp%lam)
        sp_max = maxval(nlte_patch%bb_sp%lam)
        n_windows = n_windows + 1
        if (size(windows, 1) < n_windows) then
          call expand_window_array(windows, n_windows)
        end if
        windows(n_windows, 1) = sp_min
        windows(n_windows, 2) = sp_max
      end if
      if (associated(nlte_patch%bf_sp) .and. nlte_patch%bf_sp%nnu > 0) then
        if (.not. allocated(windows)) allocate(windows(1, 2))
        sp_min = minval(nlte_patch%bf_sp%lam)
        sp_max = maxval(nlte_patch%bf_sp%lam)
        n_windows = n_windows + 1
        if (size(windows, 1) < n_windows) then
          call expand_window_array(windows, n_windows)
        end if
        windows(n_windows, 1) = sp_min
        windows(n_windows, 2) = sp_max
      end if
    else
      if (associated(nlte_patch%sp) .and. nlte_patch%sp%nnu > 0) then
        if (.not. allocated(windows)) allocate(windows(1, 2))
        sp_min = minval(nlte_patch%sp%lam)
        sp_max = maxval(nlte_patch%sp%lam)
        n_windows = n_windows + 1
        if (size(windows, 1) < n_windows) then
          call expand_window_array(windows, n_windows)
        end if
        windows(n_windows, 1) = sp_min
        windows(n_windows, 2) = sp_max
      end if
    end if
  end if

  ! Merge intervals and return them
  if (n_windows > 0 .and. allocated(windows)) then
    ! Call merge_intervals to build the union of windows with dlam expansion
    call merge_intervals(windows, n_windows, dlam, intervals_out, n_intervals_out)

    ! Also compute min/max for backward compatibility (passed as ll_lam)
    ll_lam(1) = minval(windows(:, 1))
    ll_lam(2) = maxval(windows(:, 2))
  end if

  if (allocated(windows)) deallocate(windows)

contains

  subroutine expand_window_array(arr, new_size)
    real(dp), allocatable, intent(inout) :: arr(:,:)
    integer, intent(in) :: new_size
    real(dp), allocatable :: tmp(:,:)
    integer :: old_size
    old_size = size(arr, 1)
    allocate(tmp(new_size, 2))
    tmp(1:old_size, :) = arr
    call move_alloc(tmp, arr)
  end subroutine expand_window_array

END SUBROUTINE build_linelist_intervals

!===============================================================================
SUBROUTINE create_sub_patches (self, nu, mu, sp, atoms)
  class(composition_patch_t) :: self
  integer, intent(in) :: nu(2), mu(2)
  class(atomic_spectrum_type), intent(in), pointer :: sp
  type(atom_type), intent(in), pointer, optional :: atoms(:)
  class(rt_patch_t), pointer :: rt_patch
  type(rt_patch_ptr_t), allocatable :: grown(:)
  integer :: n, m
  !-----------------------------------------------------------------------------
  if (.not. associated(sp)) return
  if (sp%nnu == 0) return
  if (.not. allocated(self%rt_slabs)) allocate(self%rt_slabs(0))
  do m=mu(1), mu(2)
    do n=nu(1), nu(2)
      allocate(rt_patch)
      call rt_patch%pre_init_and_link(self, i_mu=m, i_nu=n, sp=sp, atoms=atoms)
      allocate(grown, source=[self%rt_slabs, rt_patch_ptr_t(rt_patch)])
      call move_alloc(grown, self%rt_slabs)
    end do
  end do
END SUBROUTINE create_sub_patches

!===============================================================================
SUBROUTINE execute (self)
  class(composition_patch_t) :: self
  class(composition_type), pointer :: comp
  real(4), dimension(:,:,:), allocatable :: ne, rho, temp, vx, vy, vz, vmic
  real(4), dimension(:,:,:), allocatable :: gtx, gty, gtz, gpx, gpy, gpz
  real(4), dimension(:,:,:,:,:), allocatable :: atmos
  real(8), dimension(:,:,:), allocatable :: ne8, rho8, temp8
  integer :: i, n(3)
  logical :: need_gradients, compute_gradients
  integer, save:: itimer=0
  ! ----------------------------------------------------------------------------
  call trace%begin('composition_patch_t%execute', itimer=itimer)
  n = self%n

  ! the chain-rule 'hermite' RT solvers need the ln(T) and Pe gradients of the
  ! snapshot interpolant (the Pe gradients are of Pe itself, not its log); they
  ! are computed whenever they can be (so saved atmospheres always carry them),
  ! and stored on the compositions only when a hermite solver is actually
  ! selected
  need_gradients = (rt_out%solver == 'hermite' .or. rt_iter%solver == 'hermite' .or. &
                     rt_out%solver == 'hermite2' .or. rt_iter%solver == 'hermite2')
  compute_gradients = (need_gradients .or. save_atmos) &
                      .and. pchip_interpolation .and. .not. make_eos
  if (need_gradients .and. .not. pchip_interpolation) &
    error stop "the 'hermite'/'hermite2' solvers require pchip_interpolation=.true."

  if (self%iter == 0) then
    allocate(self%composition)
    if (rotate_atmos) allocate(self%rotated_composition(self%n_out))

    allocate(temp8(n(1),n(2),n(3)))
    allocate(rho8 (n(1),n(2),n(3)))
    allocate(ne8  (n(1),n(2),n(3)))
    allocate(temp(n(1),n(2),n(3)), ne(n(1),n(2),n(3)), rho(n(1),n(2),n(3)))
    allocate(vx  (n(1),n(2),n(3)), vy(n(1),n(2),n(3)), vz (n(1),n(2),n(3)))
    allocate(vmic(n(1),n(2),n(3)))
    allocate(gtx (n(1),n(2),n(3)), gty(n(1),n(2),n(3)), gtz(n(1),n(2),n(3)))
    allocate(gpx (n(1),n(2),n(3)), gpy(n(1),n(2),n(3)), gpz(n(1),n(2),n(3)))
    gtx = 0.; gty = 0.; gtz = 0.
    gpx = 0.; gpy = 0.; gpz = 0.

    ! a previously computed atmosphere is read back in the same 14-slot layout
    ! it was written in.  NOTE: the Pe-gradient slots now hold dPe/d (linear),
    ! not dlnPe/d -- files written before that convention change load without
    ! error but carry wrong gradients, so they must be regenerated too.
    ! a run that rotates reads the unrotated copy from a save without rotated
    ! ones, and interpolates the ray-aligned copies from the original snapshot,
    ! which read_external has loaded into ia for that case
    if (atmos_format == 'save') then
      call self%load5('atmos', run_dir=atmos_file, alloc=atmos)
      if (size(atmos,4) < self%i_dpe_z) &
        error stop 'saved atmosphere predates the gradient slots -- regenerate it'
      if (size(atmos,5) /= self%nw .and. size(atmos,5) /= 1) &
        error stop 'saved atmosphere holds rotated copies for other output directions'
      if (size(atmos,5) < self%nw .and. .not. allocated(ia%ltemp)) &
        error stop 'rotated copies of a saved atmosphere need its original snapshot'
    elseif (save_atmos) then
      allocate(atmos(n(1),n(2),n(3),self%i_dpe_z,self%nw))
    endif

    do i=1, self%nw
      if (atmos_format == 'save' .and. i <= size(atmos,5)) then
        ne   = atmos(:,:,:,self%i_ne,   i)
        rho  = atmos(:,:,:,self%i_rho,  i)
        temp = atmos(:,:,:,self%i_temp, i)
        vx   = atmos(:,:,:,self%i_vx,   i)
        vy   = atmos(:,:,:,self%i_vy,   i)
        vz   = atmos(:,:,:,self%i_vz,   i)
        vmic = atmos(:,:,:,self%i_vmic, i)
        gtx  = atmos(:,:,:,self%i_dlnt_x,  i)
        gty  = atmos(:,:,:,self%i_dlnt_y,  i)
        gtz  = atmos(:,:,:,self%i_dlnt_z,  i)
        gpx  = atmos(:,:,:,self%i_dpe_x, i)
        gpy  = atmos(:,:,:,self%i_dpe_y, i)
        gpz  = atmos(:,:,:,self%i_dpe_z, i)
      elseif (i == 1) then
        if (compute_gradients) then
          call interpolate_atmos(n, self%size, self%position, temp, ne, rho, vx, vy, vz, vmic, &
                                 dlnt_dx=gtx, dlnt_dy=gty, dlnt_dz=gtz, &
                                 dpe_dx=gpx, dpe_dy=gpy, dpe_dz=gpz)
        else
          call interpolate_atmos(n, self%size, self%position, temp, ne, rho, vx, vy, vz, vmic)
        endif
      else
        if (compute_gradients) then
          ! rotated columns run along the ray; only the z pair exists
          gtx = 0.; gty = 0.; gpx = 0.; gpy = 0.
          call interpolate_atmos(n, self%size, self%position, temp, ne, rho, vx, vy, vz, vmic, &
                                 vec=rt_out%vec(:,i-1), bot2pivot=pivot_z, &
                                 dlnt_dz=gtz, dpe_dz=gpz)
        else
          call interpolate_atmos(n, self%size, self%position, temp, ne, rho, vx, vy, vz, vmic, &
                                 vec=rt_out%vec(:,i-1), bot2pivot=pivot_z)
        endif
      endif

      ! masked once here, rather than on every ray direction and iteration
      if (.not. use_vel(1)) vx = 0
      if (.not. use_vel(2)) vy = 0
      if (.not. use_vel(3)) vz = 0

      temp8 = dble(temp)
      rho8  = dble(rho)
      ne8   = dble(ne)

      ! copy composition with all the information from the abund
      ! and absdat file to make it threadprivate
      if (i==1) then
        self%composition = composition
      else
        self%rotated_composition(i-1) = composition
      endif
      comp => self%comp_ptr(i)

      ! this is needed, otherwise all pointers will still point
      ! to ions and molecules in original composition
      call comp%reassign_pointer()
      ! loading temp, ne and rho in to composition objekt
      ! ion and molecule equilibrium are solved immediately
      if (allocated(output_set_da)) then
        call comp%load_state_variables(temp8, ne8, rho8, use_rho=use_rho, use_ne=use_ne, &
                                       vmic=vmic, vx=vx, vy=vy, vz=vz, &
                                       set_da=output_set_da, set_element=output_set_element)
      else
        call comp%load_state_variables(temp8, ne8, rho8, use_rho=use_rho, use_ne=use_ne, &
                                       vmic=vmic, vx=vx, vy=vy, vz=vz)
      endif

      if (.not. use_ne) ne = real(comp%ne)
      if (.not. use_rho) rho = real(comp%rho)

      ! with use_ne=.false. there was no Pe field to interpolate: rebuild the
      ! Pe gradients from Pe = ne*kB*T on the patch grid instead
      if (compute_gradients .and. .not. use_ne .and. atmos_format /= 'save') &
        call pe_gradients_from_grid(comp, self%ds, gpx, gpy, gpz, horizontal=(i==1))

      if (need_gradients) then
        if (i == 1) then
          call comp%set_gradients(gtz, gpz, dlnt_dx=gtx, dlnt_dy=gty, &
                                  dpe_dx=gpx, dpe_dy=gpy)
        else
          call comp%set_gradients(gtz, gpz)
        endif
      endif

      if (save_atmos) then
        atmos(:,:,:,self%i_ne,   i) = ne
        atmos(:,:,:,self%i_rho,  i) = rho
        atmos(:,:,:,self%i_temp, i) = temp
        atmos(:,:,:,self%i_vx,   i) = vx
        atmos(:,:,:,self%i_vy,   i) = vy
        atmos(:,:,:,self%i_vz,   i) = vz
        atmos(:,:,:,self%i_chi500, i) = comp%chi500
        atmos(:,:,:,self%i_vmic,   i) = vmic
        atmos(:,:,:,self%i_dlnt_x,  i) = gtx
        atmos(:,:,:,self%i_dlnt_y,  i) = gty
        atmos(:,:,:,self%i_dlnt_z,  i) = gtz
        atmos(:,:,:,self%i_dpe_x, i) = gpx
        atmos(:,:,:,self%i_dpe_y, i) = gpy
        atmos(:,:,:,self%i_dpe_z, i) = gpz
      endif

      if (check_NaN) then
        call check_array(temp, 'temperature', check_zero=.true., abort=.true.)
        call check_array(rho, 'density', check_zero=.true., abort=.true.)
        call check_array(ne, 'electron density', check_zero=.true., abort=.true.)
      endif
    end do

    call ia%release()

    if (save_atmos) then
      call self%save5(atmos, 'atmos')
    endif

    if (save_eos) then
      call self%save_eos_variables()
    endif

    if (save_species /= '') then
      call self%save_number_densities(self%composition, pre_init=.false.)
    endif

    if (rt_iter%sample_back_opac .or. rt_out%sample_back_opac) then
      call self%init_tabulated_comp(ldtemp=ldtemp, ldne=ldne)
    elseif (make_eos) then
      call self%init_eos_comp(self%composition%temp, self%composition%rho, ldtemp=ldtemp, ldrho=ldrho)
    endif

    self%dtime = max(self%dtime, io%end_time - self%dtime)
    self%do_out_sp = .true. ! so that it is still executed after convergence
  else
    self%dtime = self%patch%dtime
    !print *, self%m3d_patch_info()
    !print *, 'updating composition more than once is not implemented yet. Stopping...'
    !stop
  endif
  call trace%end(itimer)
END SUBROUTINE execute

!===============================================================================
!> Rebuild the ln(Pe) gradients from the composition's own Pe = ne*kB*T on the
!> patch grid, for paths where no interpolation-time gradients exist
!> (use_ne=.false.).  Vertical slopes use the canonical monotone pchip node
!> slopes; horizontal partials (unrotated composition only) use centered FDs
!> with one-sided edges.  Grid-node slopes lose the sub-grid structure of the
!> snapshot interpolant, but keep the hermite solver working identically.
!===============================================================================
SUBROUTINE pe_gradients_from_grid (comp, ds, gpx, gpy, gpz, horizontal)
  class(composition_type) :: comp
  real(8), intent(in) :: ds(3)
  real(4), dimension(:,:,:), intent(out) :: gpx, gpy, gpz
  logical, intent(in) :: horizontal
  real(4), dimension(:,:,:), allocatable :: pe
  real(4), dimension(:), allocatable :: zz, dxw, dyw
  integer :: i, x, y, nx, ny, nz
  integer, save :: itimer=0
  ! ----------------------------------------------------------------------------
  call trace%begin('pe_gradients_from_grid', itimer=itimer)
  nx = size(gpz,1); ny = size(gpz,2); nz = size(gpz,3)
  allocate(pe(nx,ny,nz), zz(nz), dxw(nz-1), dyw(nz-1))
  pe = real(comp%pe, 4)
  zz = [ (real(i*ds(3),4), i=1, nz) ]
  do y=1, ny
    do x=1, nx
      call pchip_slopes(zz, pe(x,y,:), 1, nz, dxw, dyw, gpz(x,y,:))
    end do
  end do
  gpx = 0.
  gpy = 0.
  if (horizontal .and. nx > 2) then
    gpx(2:nx-1,:,:) = (pe(3:nx,:,:) - pe(1:nx-2,:,:)) / real(2.*ds(1), 4)
    gpx(1,   :,:)   = (pe(2, :,:)   - pe(1,     :,:)) / real(ds(1), 4)
    gpx(nx,  :,:)   = (pe(nx,:,:)   - pe(nx-1,  :,:)) / real(ds(1), 4)
  endif
  if (horizontal .and. ny > 2) then
    gpy(:,2:ny-1,:) = (pe(:,3:ny,:) - pe(:,1:ny-2,:)) / real(2.*ds(2), 4)
    gpy(:,1,   :)   = (pe(:,2, :)   - pe(:,1,     :)) / real(ds(2), 4)
    gpy(:,ny,  :)   = (pe(:,ny,:)   - pe(:,ny-1,  :)) / real(ds(2), 4)
  endif
  call trace%end(itimer)
END SUBROUTINE pe_gradients_from_grid

!===============================================================================
!> Sample the background opacity on a grid built in the orthogonal PCA basis of
!> (lnT, lnNe) over all atmosphere cells, rather than the axis-aligned box in
!> (T, lnNe). Since T and Ne are strongly correlated the atmosphere hugs a thin
!> locus; a PC-rotated grid tracks it, dropping the unphysical off-diagonal
!> corners that crash the EOS and cutting the node count. ldtemp/ldne are the
!> requested lnT/lnNe resolutions, projected onto the PC axes for spacing.
SUBROUTINE init_tabulated_comp (self, ldtemp, ldne)
    class(composition_patch_t), target :: self
    class(composition_type), pointer :: comp
    real(8), intent(in) :: ldtemp, ldne
    integer :: i, n(3), n1, n2
    real(8), dimension(:,:,:), allocatable :: tab_temp, tab_ne, tab_rho
    real(4), dimension(:,:,:), allocatable :: tab_vmic
    real(4), dimension(:,:,:,:), allocatable :: temp, ne
    !---------------------------------------------------------------------------
    ! the PCA is taken over every cell of every rotated copy, so gather them; each
    ! copy resamples the cube along its own direction and has its own (T, ne)
    n = self%n
    allocate(temp(n(1),n(2),n(3),self%nw), ne(n(1),n(2),n(3),self%nw))
    do i=1, self%nw
      comp => self%comp_ptr(i)
      temp(:,:,:,i) = real(comp%temp)
      ne  (:,:,:,i) = real(comp%ne)
    end do

    allocate(self%tab_comp)
    self%tab_comp = composition
    call self%tab_comp%reassign_pointer()

    ! With sample_linelist, the sampled table must hold only the LTE (non-departure)
    ! species; strip the active NLTE atoms' species out of tab_comp's line list so
    ! their lines are instead computed in detail per cell (with departures). Done
    ! here, single-threaded at setup, before any concurrent table build reads it.
    if (sample_linelist .and. associated(self%nlte_patch)) then
      if (associated(self%nlte_patch%atoms)) &
        call self%tab_comp%mask_departure_species(self%nlte_patch%atoms%id)
    endif

    ! ---- PCA-rotated (T, ne) node grid, one node past the data hull; the
    ! transform is written straight into tab_comp%pca so the table fill can
    ! capture the axis metadata -------------------------------------------------
    call pca_opac_grid(temp, ne, ldtemp, ldne, 1, tab_temp, tab_ne, self%tab_comp%pca, n1, n2, &
                       use_pca=pca_opac)

    allocate(tab_rho (n1, n2, 1))
    allocate(tab_vmic(n1, n2, 1))
    tab_vmic = opac_table_vmic

    call self%tab_comp%load_state_variables(tab_temp, tab_ne, tab_rho, vmic=tab_vmic, &
                                            use_ne=.true., use_rho=.false.)
END SUBROUTINE init_tabulated_comp

!===============================================================================
SUBROUTINE init_eos_comp (self, temp, rho, ldtemp, ldrho)
  class(composition_patch_t), target :: self
  real(8), dimension(:,:,:), intent(in) :: temp, rho
  real(8), intent(in) :: ldtemp, ldrho
  real(8) :: minTemp, maxTemp, minRho, maxRho
  integer :: i, j, ntemp, nrho
  real(8), dimension(:,:,:), allocatable :: tab_temp, tab_ne, tab_rho
  real(4), dimension(:,:,:), allocatable :: eos, tab_vmic
  real(8) :: minLnRho, maxLnRho, minLnT, maxLnT
  !---------------------------------------------------------------------------
  allocate(self%tab_comp)
  self%tab_comp = composition
  call self%tab_comp%reassign_pointer()

  minTemp = minval(temp)
  maxTemp = maxval(temp)
  minRho = minval(rho)
  maxRho = maxval(rho)

  minLnRho = log(minRho)
  maxLnRho = log(maxRho)
  minLnT = log(minTemp)
  maxLnT = log(maxTemp)

  ntemp = ceiling((maxLnT - minLnT) / ldtemp) 
  nrho  = ceiling((maxLnRho - minLnRho) / ldrho) 

  allocate(tab_temp(ntemp, nrho, 1))
  allocate(tab_rho (ntemp, nrho, 1))
  allocate(tab_ne  (ntemp, nrho, 1))
  allocate(tab_vmic (ntemp, nrho, 1))

  do j=1, nrho
  do i=1, ntemp
    tab_temp(i,j,1) = minLnT + (i-1)*ldtemp
    tab_rho (i,j,1) = minLnRho + (j-1)*ldrho
  end do
  end do

  tab_rho = exp(tab_rho)
  tab_temp = exp(tab_temp)
  tab_vmic = opac_table_vmic
  call self%tab_comp%load_state_variables(tab_temp, tab_ne, tab_rho, use_rho=.true., use_ne=.false., vmic=tab_vmic)

  allocate(eos(6, ntemp, nrho))
  eos(1, :, :) = safe_to_real4(self%tab_comp%temp(:, :, 1))
  eos(2, :, :) = safe_to_real4(self%tab_comp%rho(:, :, 1))
  eos(3, :, :) = safe_to_real4(self%tab_comp%ne(:, :, 1))
  eos(4, :, :) = safe_to_real4(self%tab_comp%pg(:, :, 1))
  eos(5, :, :) = safe_to_real4(self%tab_comp%E(:, :, 1))
  eos(6, :, :) = self%tab_comp%chi500(:, :, 1)
  call self%save3(eos, 'eos')

contains
  elemental function safe_to_real4(d_val) result(s_val)
    real(8), intent(in) :: d_val
    real(4) :: s_val
    s_val = real(max(min(d_val, dble(huge(1.0_4))), -dble(huge(1.0_4))), 4)
  end function safe_to_real4
END SUBROUTINE init_eos_comp

!===============================================================================
SUBROUTINE save_eos_variables (self)
  class(composition_patch_t), target :: self
  class(composition_type), pointer :: comp
  real(4), dimension(:,:,:,:), allocatable :: eos
  integer :: nx, ny, nz
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('composition_patch%save_eos_variables', itimer=itimer)
  comp => self%composition
  nx = size(comp%temp,1)
  ny = size(comp%temp,2)
  nz = size(comp%temp,3)
  allocate(eos(nx,ny,nz,2))
  eos(:,:,:,1) = comp%pg
  eos(:,:,:,2) = comp%E
  call self%save4(eos, 'eos')
  call trace%end(itimer)
END SUBROUTINE save_eos_variables

!===============================================================================
SUBROUTINE print_post_update_info (self)
  class(composition_patch_t):: self
  class(composition_type), pointer :: comp
  integer :: p(3)
  !-----------------------------------------------------------------------------
  ! post_update also runs before the first execute, where there is nothing yet
  if (.not. associated(self%composition)) return
  if (console_on(ch_m3d, 2, all_ranks=.true.)) then
    p = self%pcheck
    comp => self%composition
    memcheck(1) = real(comp%temp(p(1),p(2),p(3)))
    memcheck(2) = real(comp%rho (p(1),p(2),p(3)))
    memcheck(3) = real(comp%ne  (p(1),p(2),p(3)))
    call self%print_to_stdout(arr=memcheck(1:3))
  endif
END SUBROUTINE print_post_update_info

!===============================================================================
SUBROUTINE save_abundances_to_file (self)
  class(composition_patch_t), target :: self
  class(element_type), pointer :: element
  real(8) :: abnd
  integer :: lunit, i
  !-----------------------------------------------------------------------------
  if (self%ipos(3)==0) then
    open(newunit=lunit, file=trim(io%run_dir)//'abund.txt', form='FORMATTED')
    do i=1, composition%nel
      element => composition%elements(i)
      abnd = log10(element%abund) + 12.0
      abnd = max(abnd, -99.9)
      write(lunit, '(A3, F8.4)') element%id, abnd
    end do
    close(lunit)
  endif
END SUBROUTINE save_abundances_to_file

!===============================================================================
SUBROUTINE save_number_densities (self, comp, pre_init)
  class(composition_patch_t), target :: self
  class(composition_type) :: comp
  class(ion_type), pointer :: ion
  class(molec_type), pointer :: mol
  class(particle_type), pointer :: particle
  integer :: i, nspecies
  logical, intent(in) :: pre_init
  character(len=256) :: local_text
  character(len=6), dimension(:), allocatable :: species
  character(len=6) :: key
  !-----------------------------------------------------------------------------
  local_text = save_species
  call count_words(local_text, nspecies, sep=',')
  allocate(species(nspecies))
  call replace_char(local_text, ' ', '*')
  call replace_char(local_text, ',', ' ')
  read(local_text, *) species

  do i=1, nspecies
    key = species(i)
    call replace_char(key, '*', ' ')
    key = adjustl(key)
    call comp%particle_by_id(trim(key), particle, verbose=.true.)
    if (associated(particle)) then
      if (pre_init) then
        particle%in_opac = .true.
      else
        call self%save3(real(particle%n), trim(key))
      endif
    endif
  end do
END SUBROUTINE save_number_densities

!===============================================================================
SUBROUTINE composition_patch_nbor_relations(self, another, needed, needs_me, download)
  class(composition_patch_t) :: self
  class(task_t), pointer :: another     !< another task in the list
  logical :: needed, needs_me, download !< Is the other task: needed,
                                           !    does it need me, should I download?
  !-----------------------------------------------------------------------------
  call trace%begin ('composition_patch_t%nbor_relations')

  if (another%kind == 'rt_column') then
    ! the columns through this slab need its composition
    needed = .false.; needs_me = .true.; download = .false.
    call trace%end()
    return
  endif

  call self%m3d_patch_t%nbor_relations(another, needed, needs_me, download)

  if (needed .or. needs_me) then
    if (another%id /= self%patch%id) then
      select case (another%kind)
      case ('nlte_patch')
        ! is needed after convergence and needs me
      case default
        needed = .false.; needs_me = .false.
      end select
    endif
  endif
  call trace%end()
END SUBROUTINE composition_patch_nbor_relations

!===============================================================================
!> Registered EOS allocations of one composition patch, in bytes.  Each
!> load_state_variables holds f, n, saha and part of every ion and molecule of
!> the patch, plus f and w of each H I level while the MHD occupation
!> probabilities iterate (transient).  n and part stay allocated for the species
!> that the opacity or the line list use, and each H I level keeps n, and w with
!> H_occupation_probability (retained).  With rotate_atmos the patch loads one
!> composition per output direction as well, one after the other.
!===============================================================================
SUBROUTINE estimate_eos_memory (self, transient, retained)
  class(composition_patch_t), intent(in) :: self
  real(8), intent(out) :: transient, retained
  class(ion_type), pointer :: HI
  real(8) :: field
  integer :: i, n_all, n_kept, n_comp, n_level
  !-----------------------------------------------------------------------------
  field = real(storage_size(0d0)/8,8)*product(real(self%n,8))
  n_all = composition%nmolec
  n_kept = count(composition%molecules%in_opac .or. composition%molecules%in_linelist)
  do i=1, composition%nel
    n_all = n_all + size(composition%elements(i)%ions)
    n_kept = n_kept + count(composition%elements(i)%ions%in_opac .or. &
                            composition%elements(i)%ions%in_linelist)
  end do
  call composition%ion_by_id('H I', HI)
  n_level = HI%nk
  n_comp = 1
  if (rotate_atmos) n_comp = 1 + self%n_out
  transient = 4*field*n_all
  if (mhd_eos) transient = transient + 2*field*n_level
  retained = 2*field*n_kept
  retained = retained + merge(2, 1, H_occupation_probability)*field*n_level
  retained = retained*n_comp
END SUBROUTINE estimate_eos_memory

END MODULE composition_patch_mod
