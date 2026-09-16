MODULE composition_patch_mod
  USE io_mod
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
  USE m3d_patch_mod,   only: m3d_patch_t, composition, comoving_frame, dlam, &
  use_rho, use_ne, overlapping_bb, overlapping_bf, binned_opac, initial_Jnu, &
  rotate_atmos, decouple_continuum, n_nu, n_nu_out, use_vel, &
  make_eos, overlapping_bb_in_output, atmos_format, atmos_file, &
  n_nlte_patches, nlte_finished
  USE rt_params_mod, only: rt_iter, rt_out, pivot_z
  USE m3d_debug_mod
  USE read_atmos_mod, only: ia, interpolate_atmos, pchip_interpolation
  USE interpol_mod, only: pchip_slopes
  USE composition_mod, only: composition_type, exclude_elements, only_elements, switch, &
  ref_lam, debug_opac, abundances, relative_abundances, mhd_eos, helium_scattering, &
  departure_pointer_t, read_departures, departure_indices, &
  sample_linelist, opac_table_vmic, &
  line_scattering
  USE atom_mod, only: use_classic_col, level_partf
  USE nlte_patch_mod,  only: nlte_patch_t
  USE rt_patch_mod, only: rt_patch_t, output_set_da, output_set_element, verify_output_patch_coverage
  USE bftrans_mod, only: H_occupation_probability
  USE bbtrans_mod, only: tolerance, bbtrans_type, disable_ga, disable_stark, disable_vdw, ga_fallback
  USE linelist_mod, only: ll_profile, voigt_kernel, voigt_kernel_init, interval_t, merge_intervals
  USE hprofl_mod, only: hydrogen_BPO
  USE hvcs_mod, only: hydrogen_vcs_table, vcs_data_dir
  USE trace_mod
  USE scaling_mod
  USE units_mod
  USE ray_mod, only: set_carlson_quadratures, one_dimensional
  USE list_mod, only: list_t
  USE link_mod, only: nbor_t
  USE bits_mod
  USE atom_mod, only: atom_type, atomic_spectrum_type
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  
  !-----------------------------------------------------------------------------
  type, public, extends(m3d_patch_t) :: composition_patch_t
    class(gpatch_t), pointer :: mhd => null() ! pointer to MHD patch
    class(nlte_patch_t), pointer :: nlte_patch => null()
  contains
    procedure :: pre_init => composition_pre_init
    procedure :: execute
    procedure :: create_sub_patches
    procedure :: print_post_update_info
    procedure :: save_abundances_to_file
    procedure :: init_tabulated_comp
    procedure :: init_eos_comp
    procedure :: save_eos_variables
    procedure :: save_number_densities
    procedure :: nbor_relations => composition_patch_nbor_relations
  end type
  
  real(8), public :: FeH=-99.9, dtemp=250., ldne=0.25, ldrho=0.25, ldtemp=0.05
  logical, public :: save_atmos=.true., save_eos=.false.
  logical, public :: pca_opac=.true.       ! .false. -> axis-aligned (lnT,lnNe) sampling grid
  real(4), public, dimension(:), allocatable :: external_abund
  real(4), public :: alpha_enhancement=0.0, r_process_enhancement=0.0
  character(len=256), save :: save_species=''
  character(len=256), save :: absdat_file='./input/TS_absdat.dat'
  character(len=256), save :: abund_file='./input/abund/abund_magg22'
  character(len=256), save :: isotope_file='./input/isotopes_blp'
  character(len=256), save :: molec_file='./input/molec_data_BC.txt'
  character(len=256), save :: opac_table=''
  character(len=256), save :: absmet_file=''
  character(len=512), save :: linelist_folder=''
  character(len=256), save, dimension(:), allocatable :: linelist_file

  public :: read_composition_params

contains

!===============================================================================
SUBROUTINE read_composition_params ()
  integer :: i, iostat
  real(4) :: tmolim=-99
  integer, save:: itimer=0, verbose=0
  namelist /composition_params/ verbose, absdat_file, abund_file, &
  isotope_file, molec_file, H_occupation_probability, &
  switch, mhd_eos, debye_shield, &
  tmolim, irwin_dis, save_species, &
  ref_lam, debug_opac, abundances, relative_abundances, alpha_enhancement, r_process_enhancement, &
  disable_ga, disable_stark, disable_vdw, helium_scattering, hydrogen_BPO, &
  hydrogen_vcs_table, vcs_data_dir, &
  ga_fallback, eosmet
  namelist /opac_table_params/ opac_table, opac_table_vmic, pca_opac, &
  dtemp, ldne, ldrho, ldtemp, absmet_file
  namelist /linelist_params/ dlam, exclude_elements, linelist_file, linelist_folder, &
  tolerance, only_elements, line_scattering, ll_profile, voigt_kernel
  ! ----------------------------------------------------------------------------
  call trace%begin('composition%init', itimer=itimer)

  allocate(linelist_file(128))
  linelist_file = ''

  rewind (io%input)
  read (io%input, composition_params, iostat=iostat)
  if (iostat > 0) then
    rewind (io%input)
    read (io%input, composition_params)
  end if
  call io%namelist_warning ('composition_params', iostat)

  write (io%namelist, composition_params)

  rewind (io%input)
  read (io%input, opac_table_params, iostat=iostat)
  if (iostat > 0) then
    rewind (io%input)
    read (io%input, opac_table_params)
  end if
  call io%namelist_warning ('opac_table_params', iostat)

  write (io%namelist, opac_table_params)

  rewind (io%input)
  read (io%input, linelist_params, iostat=iostat)
  if (iostat > 0) then
    rewind (io%input)
    read (io%input, linelist_params)
  end if
  call io%namelist_warning ('linelist_params', iostat)

  ! Build QLG-bot's lookup tables. Must precede any opacity accumulation -- see the
  ! kernel note in HOP/linelist_mod.f90.
  call voigt_kernel_init

  if (linelist_folder /= '') then
    deallocate(linelist_file)
    call os%get_filenames(linelist_folder, linelist_file)
    do i=1, size(linelist_file)
      linelist_file(i) = trim(linelist_folder)//'/'//trim(linelist_file(i))
    end do
  else
    call truncate_array(linelist_file, count(linelist_file/=''))
  endif
  write (io%namelist, linelist_params)

  if (verbose > 0) write (*,*) ' - Initializing the NEW opacity package'
  composition%verbose = verbose

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
  integer :: i, mu, nw, n_mu
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

    if (n_nu < 0) then
      n_nu = omp%nthreads
    endif
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
          ! Pass the intervals to init_linelist2 so build_linelist selects on the
          ! union of windows, not a single min/max span.
          if (allocated(ll_intervals) .and. n_intervals > 0) then
            call composition%init_linelist2(linelist_file(i), dlam, ll_lam, ll_intervals)
          else
            call composition%init_linelist2(linelist_file(i), dlam, ll_lam)
          end if
        end do
      elseif (any(linelist_file /= '')) then
        print *, 'linelist_params: no pass includes or samples the line list -- skipping line-list read'
      endif
      if (allocated(ll_intervals)) deallocate(ll_intervals)
    end block
    call nlte_patch%exclude_double_lines()
  endif

  ! make sure there are no more f-domain splits than threads
  if (output_sp%nnu > 0) n_nu = min(n_nu, output_sp%nnu)
  if (decouple_continuum) then
    if (nlte_patch%bb_sp%nnu > 0) n_nu = min(n_nu,  nlte_patch%bb_sp%nnu)
    if (nlte_patch%bf_sp%nnu > 0) n_nu = min(n_nu,  nlte_patch%bf_sp%nnu)
  else
    if (nlte_patch%sp%nnu > 0) n_nu = min(n_nu,  nlte_patch%sp%nnu)
  endif

  self%n_nu = n_nu
  n_mu = size(rt_out%vec,2)

  ! Output frequency patches.  A contiguous output grid is split like the internal
  ! one, into n_nu parts.  A window_file grid is split per window, so that no patch
  ! spans the gap between two windows (see output_patch_range).  With fewer windows
  ! than n_nu, each window is split into n_nu parts, so the output pass still has
  ! at least n_nu tasks per depth for the threads; the smallest window caps that
  ! count at one point per part.
  if (output_sp%n_segments > 0) then
    n_nu_out = output_sp%n_segments
    if (output_sp%n_segments < n_nu) n_nu_out = output_sp%n_segments * &
      min(n_nu, minval(output_sp%segment_last - output_sp%segment_first + 1))
  else
    n_nu_out = n_nu
  endif

  ! One sampled opacity table per i_nu, shared by the rt_patches created below,
  ! which differ only in i_mu and would otherwise build identical copies.  The
  ! internal spectrum uses i_nu in [1, n_nu], the bound-free half of a decoupled
  ! continuum [n_nu+1, 2*n_nu], and the output spectrum [-n_nu_out, -1]; slot 0 is
  ! never used.  make_eos is excluded: save_opac_table overwrites the table in
  ! place, so those rt_patches keep a private one.
  if ((rt_iter%sample_back_opac .or. rt_out%sample_back_opac) .and. .not. make_eos) then
    allocate(self%opac_tables(-n_nu_out:2*n_nu))
    allocate(self%opac_locks (-n_nu_out:2*n_nu))
    do i=-n_nu_out, 2*n_nu
      call self%opac_locks(i)%init('opac', i)
    end do
    ! Separate LTE line-list opacity table (sample_linelist): same i_nu keying.
    if (sample_linelist) then
      allocate(self%ll_opac_tables(-n_nu_out:2*n_nu))
      allocate(self%ll_opac_locks (-n_nu_out:2*n_nu))
      do i=-n_nu_out, 2*n_nu
        call self%ll_opac_locks(i)%init('ll_opac', i)
      end do
    endif
  endif
  ! In 1D the detailed line opacity is isotropic, so the output pass computes it
  ! once per output i_nu and set and reuses it for every direction, scattering
  ! iteration and the final pass (see stored_linelist in rt_patch_mod); a column
  ! holds few enough points to keep it.
  if (one_dimensional .and. rt_out%include_linelist .and. &
      .not. rt_out%sample_linelist) then
    allocate(self%ll_point_tables(-n_nu_out:-1))
    allocate(self%ll_point_locks (-n_nu_out:-1))
    do i=-n_nu_out, -1
      call self%ll_point_locks(i)%init('ll_point', i)
    end do
  endif
  if (associated(nlte_patch%atoms)) then
    if (rt_iter%maxiter > 0) then
      if (decouple_continuum) then
        call self%create_sub_patches([1, n_nu],      [1,2], nlte_patch%bb_sp, nlte_patch%bb_atoms)
        call self%create_sub_patches([1, n_nu]+n_nu, [1,2], nlte_patch%bf_sp, nlte_patch%bf_atoms)
      else
        call self%create_sub_patches([1, n_nu], [1,2], nlte_patch%sp, nlte_patch%atoms)
      endif
    endif
  endif

  ! rt_out%solver='none' disables the output pass entirely -- create no output
  ! RT patches (the final i_mu=-1 ray or the Lambda-iteration i_mu=[1,2] ones).
  if (output_sp%nnu > 0 .and. rt_out%solver /= 'none') then
    if (n_mu > 0 .or. (make_eos .and. self%key(1)==1)) then
      if (rt_out%maxiter > 0) then
        call self%create_sub_patches([-n_nu_out, -1], [1,2], output_sp, nlte_patch%output_atoms)
      endif
      call self%create_sub_patches([-n_nu_out, -1], [-1,-1], output_sp, nlte_patch%output_atoms)
    endif
  endif

  if (output_sp%n_segments > 0) then
    call verify_output_patch_coverage(output_sp, n_nu_out)
  endif

  if (.not. associated(self%nlte_patch)) then
    deallocate(nlte_patch)
  endif

  if (save_species /= '') then
    call self%save_number_densities(composition, pre_init=.true.)
  endif

  call self%save_abundances_to_file()
  call self%post_pre_init()
  call trace%end(itimer)
END SUBROUTINE composition_pre_init

!===============================================================================
! Phase B: Build line-list intervals from output and internal spectra.
! When the output spectrum has multiple segments (multi-window), use those
! segment bounds as intervals. For non-segmented internal spectra, use the
! full range as a single interval. Expand each interval by dlam.
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
  type(interval_t), allocatable :: temp_intervals(:)
  real(dp) :: sp_min, sp_max
  integer :: i, n_windows, nout
  logical :: output_uses_ll, internal_uses_ll
  !-----------------------------------------------------------------------------

  ! Determine which spectra use the line list
  output_uses_ll = (rt_out%include_linelist .or. rt_out%sample_linelist) .and. &
                   (size(rt_out%vec,2) > 0 .or. make_eos)
  internal_uses_ll = (rt_iter%include_linelist .or. rt_iter%sample_linelist) .and. &
                     (rt_iter%maxiter > 0)

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
  integer :: n, m
  logical :: new_nlte
  !-----------------------------------------------------------------------------
  if (.not. associated(sp)) return
  if (sp%nnu == 0) return
  do m=mu(1), mu(2)
    do n=nu(1), nu(2)
      allocate(rt_patch)
      call rt_patch%pre_init_and_link(self, i_mu=m, i_nu=n, sp=sp, atoms=atoms)
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
    stop "the 'hermite'/'hermite2' solvers require pchip_interpolation=.true."

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
    if (atmos_format == 'save') then
      call self%load5('atmos', run_dir=atmos_file, alloc=atmos)
      if (size(atmos,4) < self%i_dpe_z) &
        stop 'saved atmosphere predates the gradient slots -- regenerate it'
    elseif (save_atmos) then
      allocate(atmos(n(1),n(2),n(3),self%i_dpe_z,self%nw))
    endif

    do i=1, self%nw
      if (atmos_format == 'save') then
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
SUBROUTINE init_tabulated_comp (self, ldtemp, ldne)
    ! Sample the background opacity on a grid built in the orthogonal PCA basis of
    ! (lnT, lnNe) over all atmosphere cells, rather than the axis-aligned box in
    ! (T, lnNe). Since T and Ne are strongly correlated the atmosphere hugs a thin
    ! locus; a PC-rotated grid tracks it, dropping the unphysical off-diagonal
    ! corners that crash the EOS and cutting the node count. ldtemp/ldne are the
    ! requested lnT/lnNe resolutions, projected onto the PC axes for spacing.
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
  integer :: i, j, ntemp, nrho, nnu
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
  if (self%verbose > 1) then
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
  class(task_t), pointer :: another         ! another task in the list
  logical :: needed, needs_me, download    ! Is the other task: needed, 
                                           !    does it need me, should I download?
  !-----------------------------------------------------------------------------
  call trace%begin ('composition_patch_t%nbor_relations')

  call self%m3d_patch_t%nbor_relations(another, needed, needs_me, download)

  if (needed .or. needs_me) then
    if (another%id /= self%patch%id) then
      select case (another%kind)
      case ('nlte_patch')
        ! is needed after convergence and needs me
      case ('rt_patch')
        needed = .false.
      case default
        needed = .false.; needs_me = .false.
      end select
    endif
  endif
  call trace%end()
END SUBROUTINE composition_patch_nbor_relations

END MODULE composition_patch_mod
