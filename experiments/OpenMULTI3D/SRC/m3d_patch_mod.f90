MODULE m3d_patch_mod
    USE os_mod
    USE list_mod
    USE halo_mod
    USE gpatch_mod
    USE link_mod
    USE task_mod
    ! USE guards_mod
    USE snapshots_mod
    USE mesh_mod
    USE io_mod
    USE data_io_mod
    USE timer_mod
    USE trace_mod
    USE mpi_mod
    USE omp_mod
    USE omp_lock_mod
    USE shared_mod
    USE omp_timer_mod
    USE math_mod
    USE bits_mod
    USE const_mod
    USE handyfuncs_mod
    USE prec_mod, only: inf4
    USE patch_mod, only: patch_t, patch_t_init
    USE m3d_halo_mod, only: m3d_halo_t, print_finished, lock_mem, unlock_mem, print_lock, &
                            add_nbor_by_key
    USE spectrum_mod, only: lam_tolerance
    USE atom_mod, only: atom_type, atomic_spectrum_type, lse_solver, &
                        nq_fac, cswitch, ng_step, ng_start, ng_global, save_rate_matrix, &
                        max_update, min_pop_frac, conv_pop_floor, sor_w
    USE ng_mod, only: ng_rate_min, ng_rate_tol, ng_clamp_max, ng_min_iters_left, &
                      ng_cond_min, ng_max_change, ng_min_gap
    USE composition_mod, only: composition_type, departure_pointer_t
    USE ray_mod, only: one_dimensional
    USE rt_params_mod, only: rt_iter, rt_out, read_rt_params, echo_rt_params, &
                             setup_quadrature, ray_offset, pivot_z
    USE m3d_debug_mod
    USE bftrans_mod, only: h_occupation_probability
    USE element_mod, only: element_type
    USE ion_mod, only: ion_type, particle_type
    USE line_types_mod, only: linelist_type
    USE linelist_mod, only: print_ll_atom
    USE molecule_mod, only: molec_type
    USE optab_mod, only: sample_optab_type
    USE console_mod
    USE, intrinsic :: ieee_arithmetic
    implicit none
    private

    type, public :: key_link_t
      integer :: key(4)
      class(key_link_t), pointer :: next => null()
    end type key_link_t

    ! Detailed line-list increments to (xhat, ehat) at one output patch's own
    ! points for one set, as add_linelist_opacities returns them, plus the
    ! physical chi_scat when line_scattering is on.  Shape (fs:fe,nx,ny,nz).
    type, public :: point_opac_set_t
      real(4), dimension(:,:,:,:), allocatable :: chi, eta, scat
      logical :: initialised=.false.
    end type point_opac_set_t

    type, public :: point_opac_t
      type(point_opac_set_t), allocatable :: set(:)   ! (n_set)
    end type point_opac_t

    ! LTE background opacity (chi_cont, sigma) of one 1D rt_patch slice, as
    ! lte_background_opacities returns it.  Shape (fs:fe,nx,ny,nz).
    type, public :: bg_opac_t
      real(4), dimension(:,:,:,:), allocatable :: chi, sigma
      logical :: initialised=.false.
    end type bg_opac_t
  
    !-----------------------------------------------------------------------------
  
    type, public, extends(gpatch_t) :: m3d_patch_t
      class(gpatch_t),    pointer :: patch => null()         ! pointer to mhd patch
      class(m3d_patch_t), pointer :: m3d_parent => null()    ! pointer to parent patch
      type(composition_type), pointer :: composition => null()
      type(composition_type), pointer :: tab_comp => null()
      type(composition_type), pointer :: rotated_composition(:) => null()
      ! Sampled background-opacity tables, owned by the composition patch and
      ! shared by all its rt_patches.  Everything a table's content depends on is
      ! a function of i_nu alone (the spectrum, the fs:fe slice, include_absmet),
      ! so i_nu is the key; the rt_patches that differ only in i_mu would
      ! otherwise each fill a byte-identical copy.  One lock per table lets
      ! different i_nu build concurrently while twins wait for the first builder.
      type(sample_optab_type), pointer :: opac_tables(:) => null()   ! indexed by i_nu
      type(lock_t),            pointer :: opac_locks(:)  => null()   ! same bounds
      ! Separate sampled LTE line-list opacity table (sample_linelist), holding
      ! only non-departure species; same i_nu keying and lock pattern as opac.
      type(sample_optab_type), pointer :: ll_opac_tables(:) => null()   ! indexed by i_nu
      type(lock_t),            pointer :: ll_opac_locks(:)  => null()   ! same bounds
      ! Detailed line-list opacity of the 1D output pass, stored by the first
      ! foreground of each (i_nu, set) and copied back by the later scattering
      ! iterations and the final pass; indexed by the output i_nu in [-n_nu_out,-1].
      type(point_opac_t), pointer :: ll_point_tables(:) => null()
      type(lock_t),       pointer :: ll_point_locks(:)  => null()   ! same bounds
      ! LTE background opacity of a 1D atmosphere, stored by the first
      ! lte_background_opacities of each i_nu and copied back by its twins and
      ! later iterations; indexed like opac_tables, [-n_nu_out, 2*n_nu].
      type(bg_opac_t), pointer :: bg_point_tables(:) => null()
      type(lock_t),    pointer :: bg_point_locks(:)  => null()   ! same bounds
      class(atomic_spectrum_type), pointer :: sp => null()
      type(atom_type), dimension(:), pointer :: atoms => null()
      class(key_link_t), pointer :: key_list => null()
      type(departure_pointer_t), dimension(:,:), pointer :: departures
      real(4), dimension(:,:,:,:,:,:), pointer :: authority_mem => null()
      logical :: in_task_list = .false.
      logical :: initialized  = .false.
      real(8) :: t_offset=0
      integer :: iter=0, maxiter, master_maxiter, dir=0
      integer :: n_nu, n_out, zs, ze
      integer :: fcheck, pcheck(3)
      logical :: do_out_sp=.false., do_full_sp=.false., do_final=.false.
      logical :: wait_for_nlte=.false.
      logical :: compressed=.false., save_output=.false., saved=.false.
      logical :: print_update_info, checked_bounds=.false.
      logical :: has_atom=.false.
      integer :: i_nu=0, i_mu=0
      integer :: i_ne=1, i_rho=2, i_temp=3, i_vx=4, i_vy=5, i_vz=6
      ! remaining atmos-save slots: chi500, vmic, then the ln(T) and Pe
      ! gradients (the Pe slots hold dPe/d, not dlnPe/d)
      integer :: i_chi500=7, i_vmic=8
      integer :: i_dlnt_x=9, i_dlnt_y=10, i_dlnt_z=11
      integer :: i_dpe_x=12, i_dpe_y=13, i_dpe_z=14
      logical :: prepend=.false.
    contains
      procedure :: pre_init
      procedure :: post_pre_init
      procedure :: init
      procedure :: pre_update
      procedure :: update
      procedure :: execute
      procedure :: post_update
      procedure :: finish
      procedure :: is_finished
      procedure :: nbor_relations => m3d_nbor_relations
      procedure :: pre_init_and_link
      procedure :: add_to_key_list
      procedure :: set_nv_and_nw
      procedure :: set_initial_time
      procedure :: set_maxiter
      procedure :: init_bdries
      procedure :: init_level
      procedure :: task_info
      procedure :: append_to_task_list
      procedure :: extrapolate_guards
      procedure :: is_ahead_of
      procedure :: rotate
      procedure :: comp_ptr
      procedure :: allocate_mem
      procedure :: check_mem_NaN
      procedure :: m3d_patch_info
      procedure :: m3d_patch_info2
      procedure :: print_to_stdout
      procedure :: print_post_update_info
      procedure :: check_in_bounds
      procedure :: init_nbors_hash
      procedure :: check_nbor
      procedure :: hash_key
      procedure :: set_place
      procedure :: queue_before
      procedure :: sort_nbors
      procedure :: rotate_cube4
      procedure :: collect_line_diagnostics
    end type m3d_patch_t
  
    !-----------------------------------------------------------------------------
    integer, save, public:: verbose=console_inherit
    integer, public :: nz_total
    integer, public, save :: save_id=-1
    integer, public, save :: n_nu=-1
    integer, public, save :: n_nu_out=0   ! output frequency patch count, set in composition_pre_init
    real(8), public:: convlim=1d-3, convmax=tiny(1.0)
    logical, public, save :: converged=.false., use_vel(3)=.true., premature_finish=.true.
    ! Number of nlte_patches created on this rank, and whether all of them have
    ! converged.  Defaults to "nothing to wait for", so that a run without atoms,
    ! which creates no nlte_patch at all, is not held up (see is_ahead_of).
    integer, public, save :: n_nlte_patches=0
    logical, public, save :: nlte_finished=.true.
    logical, public, save :: flux_convergence=.false.
    logical, public, save :: preconditioning=.true., binned_opac=.false.
    ! Which preconditioning the bound-bound rates use.  'original' subtracts the
    ! local source term from the rate after the frequency integral, which for a
    ! thermalised thick line leaves no significant digits and can stall the
    ! iteration outright.  'largest' uses Rybicki & Hummer (1992) Sect. 2.4
    ! instead, taking the line's own share out of the local source term rather
    ! than out of the integrated rate.  Both have the same fixed point.  It needs
    ! the linear iteration solver, which is the only one that keeps the formal
    ! solution's non-local part; rt_patch falls back to 'original' otherwise.
    character(len=16), public, save :: rate_formulation='largest'
    ! Olson-Auer-Buchler local ALI on the continuum-scattering source term.
    ! Without it, eta = chi_tot*Bnu + sigma_scat*(J_old - Bnu) is a plain,
    ! unaccelerated Lambda-iteration that converges at the scattering albedo per
    ! iteration -- i.e. barely at all where the albedo approaches 1.  Off by
    ! default so existing runs are bit-identical; see fe607_convergence.md.
    logical, public, save :: scattering_preconditioning=.false.
    ! Cap on the dimensionless self-coupling op = sigma_scat*PsiJ, which sets the
    ! amplification 1/(1-op) the OAB update applies.  op <= 1 holds by
    ! construction now that sigma is carried as the scattering *part of*
    ! chi_cont on both the exact and the sampled path (see composition_mod
    ! opacity/disjoint and optab_mod sample_interpolate), and Lambda_diag <= 1.
    ! op does legitimately approach 1 in the bound-free half of the internal
    ! spectrum, where continuum scattering is nearly the whole opacity and
    ! Psi_diag -> 1, so the local operator is genuinely near-singular there and
    ! the cap bounds the amplification; it is not a roundoff guard.
    ! 0.5 caps the amplification at 2.  It is an empirical setting from a 1x1
    ! column sweep (converges at 0.2/0.5/0.7, diverges at 0.8/0.9), not a bound
    ! derived from the scheme, and it has not been re-measured in 3D.
    ! Irrelevant unless scattering_preconditioning=T.
    real(8), public, save :: scat_op_max=0.5d0
    ! Handling of population inversions (negative line opacity).  In NLTE a line
    ! contributes phi*Rij*(bi - bj*beta) to chi, which is negative for an inverted
    ! population.  chi can then pass through zero (Snu=eta/chi diverges) and dtau
    ! turns negative (exp(+|dtau|) amplification).  Both are physical -- Mg I does
    ! invert its Rydberg levels -- but unbounded in real(4).  These parameters keep
    ! the maser but bound it, so that the population feedback can quench it between
    ! iterations instead of the arithmetic overflowing first.
    real(8), public, save :: maser_chi_floor=1d-4  ! |chi| floor, as a fraction of xhat
    real(8), public, save :: maser_dtau_max=2d0    ! cap on the per-cell amplification
    real(8), public, save :: maser_Inu_max=1d30    ! absolute cap on |Inu|
    logical, public, save :: suppress_masers=.false. ! instead force chi>0 (kills the maser)
    logical, public, save :: save_resolved=.false.
    logical, public, save :: decouple_continuum=.true.
    logical, public, save :: use_rho=.true., use_ne=.true.
    logical, public, save :: overlapping_bb=.true., overlapping_bf=.true.
    logical, public, save :: overlapping_bb_in_output=.true., detailed_absmet=.false.
    logical, public, save :: column_by_column=.false.
    logical, public, save :: print_start=.false., print_send=.false.
    logical, public, save :: Jnu_to_file=.true.
    character, public, save :: amr='P'
    character(len=6), public, save :: initial_Jnu='Planck'
    logical, public, save :: rotate_atmos=.false.
    logical, public, save :: make_eos=.false.
    logical, public, save :: save_cntrbf=.false., depression_cntrbf=.true.
    logical, public, save :: save_Inu=.false., save_chi=.false.
    logical, public, save :: save_Snu=.false., save_Jnu=.false.
    real(8), public, save :: dlam=100.
    real(8), public, dimension(:), allocatable :: rank_lvlerr
    character(len=128), save :: save_patch_kind=''
    ! Optional input manifest for an extended output set axis (abundance ladder).
    ! It must live outside run_dir, whose text outputs are cleared at startup.
    character(len=512), public, save :: set_file=''
    real(4), public :: fixed_Jnu=-1.0
    ! Output-content selection formerly read via &atom_params
    logical, public, save :: save_populations=.true.
    logical, public, save :: bound_free_in_output=.false.
    integer, public, save :: save_rates=0
    character(len=128), public :: atmos_format='unknown'
    character(len=512), public :: atmos_file='./input/atmos/p5777_g+4.4_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod'
    type(mpi_counter_t), public :: atom_update_counter, global_maxerr, global_pos
    type(mpi_counter_t), public :: flux_update_counter, global_flux_err, global_flux_pos, global_flux_lambda
    ! flux_conv_signal: raised (>0) by rt_patch, cross-rank, once the emergent
    ! flux has converged. rt_patch computes the flux AFTER nlte_patch has already
    ! solved the populations against the previous radiation field, so rather than
    ! stopping immediately we signal nlte_patch to perform one final population
    ! solve against the converged field. flux_final_counter counts the nlte
    ! patches that have done that final solve; the last one (== dims(3)) sets
    ! the global `converged` flag, which then winds the run down as usual.
    type(mpi_counter_t), public :: flux_conv_signal, flux_final_counter
    type(composition_type), public, save, target :: composition
    ! the output atoms whose line lists write_linelist_out exports, and whether any
    ! final-pass patch has merged its diagnostics (set by collect_line_diagnostics)
    type(atom_type), dimension(:), pointer, save :: linelist_atoms => null()
    logical, save :: have_line_diagnostics = .false.
      
    public :: hash_string, pre_update, update, print_execution_history
    public :: release_output_tasks, write_linelist_out

    namelist /nlte_params/ preconditioning, rate_formulation, &
    scattering_preconditioning, &
    scat_op_max, decouple_continuum, &
    flux_convergence, initial_Jnu, fixed_Jnu, &
    ng_step, ng_start, ng_global, cswitch, lse_solver, &
    ng_rate_min, ng_rate_tol, ng_clamp_max, ng_min_iters_left, ng_cond_min, &
    ng_max_change, ng_min_gap, &
    convlim, convmax, premature_finish, &
    max_update, min_pop_frac, conv_pop_floor, sor_w, &
    maser_chi_floor, maser_dtau_max, maser_Inu_max, suppress_masers

    namelist /output_params/ save_id, save_patch_kind, save_resolved, &
    save_cntrbf, depression_cntrbf, save_Inu, save_chi, save_Snu, save_Jnu, &
    Jnu_to_file, overlapping_bb_in_output, &
    save_populations, save_rates, save_rate_matrix, bound_free_in_output, set_file

    namelist /m3d_debug_params/ verbose, fcheck, linecheck, lvlcheck, pcheck, &
    atomcheck, check_NaN, print_finished, print_lock, print_start, print_send, &
    diag_pops, diag_db, diag_no_scat_loop

    namelist /run_params/ n_nu, make_eos, rotate_atmos, &
    column_by_column, lam_tolerance

  contains
  
  !===============================================================================
  SUBROUTINE pre_init (self)
    class(m3d_patch_t), target:: self
    integer :: iostat, m(6), i
    logical, save :: first_time=.true.
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%pre_init', itimer=itimer)
    self%type = trim(self%kind)
    !-----------------------------------------------------------------------------
    ! Read parameters
    !-----------------------------------------------------------------------------
    call io%lock%set ('m3d_patch_t%pre_init')
    if (first_time) then
      first_time = .false.
      rewind (io%input)
      read (io%input, nlte_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, nlte_params)
      end if
      call io%namelist_warning ('nlte_params', iostat)
      select case (rate_formulation)
      case ('original', 'largest')
      case default
        error stop 'rate_formulation must be original or largest'
      end select
      rewind (io%input)
      read (io%input, output_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, output_params)
      end if
      call io%namelist_warning ('output_params', iostat)
      rewind (io%input)
      read (io%input, m3d_debug_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, m3d_debug_params)
      end if
      call io%namelist_warning ('m3d_debug_params', iostat)
      rewind (io%input)
      read (io%input, run_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, run_params)
      end if
      call io%namelist_warning ('run_params', iostat)
      call read_rt_params ()
      rotate_atmos = rotate_atmos .and. .not. (one_dimensional .or. column_by_column)
      !-------------------------------------------------------------------------
      ! Assign angles and integration weights
      !-------------------------------------------------------------------------
      call setup_quadrature (rt_iter, is_output=.false., ds=self%patch%ds)
      call setup_quadrature (rt_out,  is_output=.true.,  ds=self%patch%ds)

      call save_angle_quadrature()

      if (pivot_z < 0) then
        ! taking the top as the pivot plane
        pivot_z = self%patch%box(3)
      endif

      Jnu_to_file = (Jnu_to_file .and. .not. one_dimensional) .or. save_Jnu

      if (ng_step > 0) then
        if (ng_step < 3) then
          call console_warn (ch_m3d, 'NG acceleration needs 3 initial iterations; ng_step set to 3')
          ng_step = 3
        endif
      endif
      call console_set_level (verbose, ch_m3d)
    end if

    self%n_out = size(rt_out%vec, 2)
    self%do_final = self%i_mu < 0 ! this is the ray working on the final output spectrum after scattering
    self%do_out_sp = self%i_nu < 0 .or. self%i_mu < 0 ! this is any ray working on the output spectrum
    self%do_full_sp = .not. self%do_out_sp ! this is any ray working on the internal spectrum for the iteration
    ! The output-spectrum rays must not start when their own column converges, but
    ! only once every nlte_patch has; see is_ahead_of().  Taken here, rather than
    ! from do_out_sp directly, because composition_patch sets do_out_sp after its
    ! first execute, and those tasks are not to be held back.
    self%wait_for_nlte = self%do_out_sp
  
    if (self%pre_initialized .and. console_on(ch_m3d, 3)) &
      write(console_unit,'(a)') 'pre_init issued somewhere'
    call io%lock%unset ('m3d_patch_t%pre_init')
  
    call self%patch_t%pre_init()

    self%id = 0 ! this is assigned later
    self%in_task_list = .false.
  
    ! ----------------------------------------------------------------------------
    ! relevant synchronization items for cartesian setup
    ! ----------------------------------------------------------------------------
    self%n          = self%patch%n
    self%ng         = 0
    self%ipos       = self%patch%ipos
    self%rank       = self%patch%rank
    self%box        = self%patch%box
    self%size       = self%patch%size
    self%position   = self%patch%position
    self%ds         = self%patch%ds
    self%origin     = self%patch%origin
    self%llc_cart   = self%patch%llc_cart
    self%centre_nat = self%patch%centre_nat
    self%llc_nat    = self%patch%llc_nat
    self%status     = self%patch%status ! be aware that this copies all bits
    self%periodic   = self%patch%periodic
    self%mesh       => self%patch%mesh
    self%task_list  => self%patch%task_list

    self%ze = self%zs + self%n(3)

    select case (self%kind)
    case ('nlte_patch')
      self%save_output = .true.
    case default
      self%save_output = string_in_text(save_patch_kind, self%kind)
    end select

    self%save_output = self%save_output .or. self%id == save_id

    ! ----------------------------------------------------------------------------
    ! Do not check density or validate with it. This is rather important to be set
    ! ----------------------------------------------------------------------------
    call self%set (bits%no_density)
    call self%set (bits%no_filled)
    
    call self%set_place(self, -2)
    call self%init_level()

    self%nt = 1
    self%time  = self%patch%time
    self%dtime = self%patch%dtime
    self%dt_fixed = self%patch%dtime
    self%out_next = self%patch%out_next
    call self%task_t%init()
    self%lock%kind = 'm3d'
    self%lock%id = self%id
    self%key = self%hash_key() ! overwrite usual hash key

    if (.not. associated(self%mem_lock)) then
      allocate (self%mem_lock(self%nt))
      do i=1,self%nt
        call self%mem_lock(i)%init ('meml', self%id)
      end do
    endif
    call self%save_lock%init ('save', self%id)

    call self%set_maxiter()
    call trace%end (itimer)
  END SUBROUTINE pre_init

  !===============================================================================
  SUBROUTINE post_pre_init (self)
    class(m3d_patch_t) :: self
    logical, save :: first_time=.true.
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('m3d_patch_t%post_pre_init', itimer=itimer)
    if (first_time) then
      ! this is done here, so that overwritten parameters are saved accordingly
      write (io%namelist, nlte_params)
      write (io%namelist, output_params)
      write (io%namelist, m3d_debug_params)
      write (io%namelist, run_params)
      call echo_rt_params ()
      first_time = .false.
    endif
    !-----------------------------------------------------------------------------
    ! Preliminary region information, relative to the task position.  This may
    ! change in patch_t%init, but the settings below allow overlap() calls to work
    !-----------------------------------------------------------------------------
    self%authority%l = - self%size*0.5_8
    self%authority%u = + self%size*0.5_8
    self%interest%l  = self%authority%l - self%ng*self%size/self%n
    self%interest%u  = self%authority%u + self%ng*self%size/self%n

    self%gn = self%n + 2*self%ng
    if (self%ng(3)==0) then
      call self%set(bits%no_download)
    endif

    call self%init_bdries()
    call self%set_initial_time()
    call self%append_to_task_list()
    call trace%end(itimer)
  END SUBROUTINE post_pre_init

  !===============================================================================
  recursive SUBROUTINE pre_init_and_link (self, m3d_parent, i_mu, i_nu, sp, atoms)
    class(m3d_patch_t), target :: self, m3d_parent
    class(atomic_spectrum_type), pointer, optional :: sp
    type(atom_type), pointer, optional :: atoms(:)
    integer, optional :: i_mu, i_nu
    !-----------------------------------------------------------------------------
    if (present(i_mu)) self%i_mu = i_mu
    if (present(i_nu)) self%i_nu = i_nu
    if (present(sp)) self%sp => sp
    if (present(atoms)) self%atoms => atoms
    self%has_atom = associated(self%atoms)
    self%patch => m3d_parent%patch
    self%m3d_parent => m3d_parent
    self%zs = m3d_parent%zs
    call self%pre_init()
    call m3d_parent%add_to_key_list(self)
  END SUBROUTINE pre_init_and_link

  !===============================================================================
  SUBROUTINE set_nv_and_nw (self, nv, nw)
    class(m3d_patch_t) :: self
    integer :: nv, nw
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%set_nv_and_nw')
    self%nv = nv
    self%nw = nw
    self%mv = self%nv
    call trace%end()
  END SUBROUTINE set_nv_and_nw
  
  !===============================================================================
  SUBROUTINE init (self)
    class(m3d_patch_t) :: self
    integer :: li(3), ui(3), status
    integer, save:: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%init', itimer=itimer)
    ! internal is set in task%pre_init, which is called by patch%init
    status = self%status
    call patch_t_init(self) ! mem is allocated here
    call data_io%init (self)
    self%initialized=.true.
    self%status = status

    if (self%is_set(bits%no_download)) then
      ! ----------------------------------------------------------------------------
      ! Set unsigned and pervolume.
      ! Although quantities are generally per volume, but density is not stored
      ! locally, thus this flag must not be set.
      ! opacities can become negative in rare cases (masers and lasers)!
      ! ----------------------------------------------------------------------------
      self%pervolume = .false.
      self%unsigned  = .false.
    endif

    self%pcheck = pcheck
    self%fcheck = fcheck
    
    where (pcheck < 0)
      self%pcheck = pcheck + [self%n(1:2), nz_total] + 1
    endwhere
    self%pcheck(3) = self%pcheck(3) - self%zs + 1
    call trace%end(itimer)
  END SUBROUTINE init

  !===============================================================================
  SUBROUTINE set_initial_time (self)
    class(m3d_patch_t) :: self
    real(8) :: tnu
    integer :: thread
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%set_initial_time')
    if (self%do_out_sp) then
      if (self%do_final) then
        self%time = self%time + rt_iter%maxiter + rt_out%maxiter
      else
        self%time = self%time + rt_iter%maxiter
      endif
    endif

    if (self%dir/=0) then
      self%t_offset = 0.5d0
    else
      self%t_offset = 0.1d0
    endif

    self%time = self%time + self%t_offset * self%dtime

    self%t = self%time            ! in case already set
    self%iit = 1                   ! initial slots
    self%it = 1                    ! slot for current values
    self%new = 1                   ! slot for new values
    self%iit(self%nt) = self%new   ! slot for new values
    self%time = self%t(self%it)    ! initial time

    call trace%end()
  END SUBROUTINE set_initial_time

  !===============================================================================
  SUBROUTINE set_maxiter (self)
    class (m3d_patch_t) :: self
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%set_maxiter')
    self%master_maxiter = rt_iter%maxiter + rt_out%maxiter

    select case (self%kind)
    case ('composition')
      self%maxiter = 1
    case ('nlte_patch')
      self%maxiter = rt_iter%maxiter
    case ('rt_patch', 'rt_column')
      self%maxiter = rt_iter%maxiter - 1
    end select

    if (self%do_out_sp) then
      if (self%do_final) then
        self%maxiter = 0
      else
        self%maxiter = rt_out%maxiter - 1
      endif
    endif
    call trace%end()
  END SUBROUTINE set_maxiter
  
  !===============================================================================
  !> Initialize boundary bits
  !===============================================================================
  SUBROUTINE init_bdries (self)
    class(m3d_patch_t):: self
    !-----------------------------------------------------------------------------
    if (.not. associated (self%boundaries)) then
      allocate (self%boundaries)
      call self%boundaries%pre_init (self)
    endif
    if (self%do_final .and. rt_iter%maxiter==0) then
      call self%boundaries%clear(bits%zu)
    endif
  END SUBROUTINE init_bdries

  !===============================================================================
  SUBROUTINE pre_update (self)
    class(m3d_patch_t), target :: self
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%pre_update', itimer=itimer)
    if (print_start) then
      write(console_unit, *) 'Update on', self%id, 'has begun', self%finished, self%time
    endif
    if (.not. self%is_finished()) then
      if (converged .and. self%do_full_sp) then
        call self%finish()
        self%print_update_info = .false.
        call print_execution_history(self, 'SKP')
      else
        call print_execution_history(self, 'PRE')
        if (.not. self%initialized) then
          call self%init()
        endif
      endif
    endif
    if (.not. self%checked_bounds) then
      self%print_update_info = self%check_in_bounds()
      self%checked_bounds = .true.
    endif
    self%place = 0
    call trace%end(itimer)
  END SUBROUTINE pre_update

  !===============================================================================
  SUBROUTINE update (self)
    class(m3d_patch_t) :: self
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('m3d_patch_t%update', itimer=itimer)
    call lock_mem(self, 'update')
    if (.not. self%is_finished()) then
      call self%execute()
      call self%post_update()
    endif
    call unlock_mem(self, 'update')
    call trace%end(itimer)
  END SUBROUTINE update

  !===============================================================================
  SUBROUTINE execute (self)
    class(m3d_patch_t) :: self
    !-----------------------------------------------------------------------------
    ! fallback
  END SUBROUTINE execute
  
  !===============================================================================
  SUBROUTINE post_update (self)
    class(m3d_patch_t) :: self
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%post_update', itimer=itimer)
    if (self%kind == 'nlte_patch' .and. converged) then
      call print_execution_history(self, 'CON')
    else
      call print_execution_history(self, 'END')
    endif
    if (self%print_update_info) then
      call self%print_post_update_info()
    elseif (console_on(ch_m3d, 3)) then
      call self%print_to_stdout()
    endif
    call self%check_mem_NaN()
    call trace%end(itimer)
  END SUBROUTINE post_update

  !===============================================================================
  SUBROUTINE save_angle_quadrature ()
    integer :: lunit, iostat, i, nmu
    character(len=256) :: meta_file
    integer, save :: itimer=0
    logical, save :: first_time = .true.
    !-----------------------------------------------------------------------------
    if (first_time) then
      first_time = .false.
      call trace%begin('m3d_patch_t%save_angle_quadrature', itimer=itimer)
      meta_file = trim(io%run_dir) // 'angle_quadrature.txt'
      open(newunit=lunit, file=meta_file, form='formatted')
      nmu = size(rt_out%wmu)
      write(lunit, '(A3, 4(A12))') '#', 'weight', 'mux', 'muy', 'muz'
      do i=1, nmu
        write(lunit, '(I3, F12.8, SP, 3(F12.8))') i, rt_out%wmu(i), rt_out%vec(:,i)
      end do
      close(lunit)
      call trace%end(itimer)
    endif
  END SUBROUTINE save_angle_quadrature

  !===============================================================================
  SUBROUTINE rotate (self)
    class(m3d_patch_t):: self
    !-----------------------------------------------------------------------------
    self%iter = self%iter + 1
    if (self%iter > self%maxiter) then
      call self%finish()
    endif
    call self%task_t%rotate
  END SUBROUTINE rotate

  !===============================================================================
  FUNCTION m3d_patch_info (self) result(text)
    class(m3d_patch_t) :: self
    character(len=48) :: text
    integer :: i
    !-----------------------------------------------------------------------------
    i = omp%place
    if (console_on(ch_m3d, 3)) then
      write(text, '(I3, 1x, I3, 1x, F8.4, A12, I8, 3(I4))') i, self%place, self%time, &
                      trim(self%kind), self%id, self%ipos(3)+1, self%i_nu, self%i_mu
    else
      write(text, '(A12, I8, 3(I4))')  trim(self%kind), self%id, &
      self%ipos(3)+1, self%i_nu, self%i_mu
    endif
  END FUNCTION m3d_patch_info

  !===============================================================================
  FUNCTION m3d_patch_info2 (self) result(text)
    class(m3d_patch_t) :: self
    character(len=38) :: text
    !-----------------------------------------------------------------------------
    write(text, '(A12, I8)') trim(self%kind), self%id
  END FUNCTION m3d_patch_info2
  
  !===============================================================================
  SUBROUTINE finish (self)
    class(m3d_patch_t):: self
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%finish')
    self%time = io%end_time
    self%t = io%end_time
    call self%set_place(self, -2)
    call self%set(bits%no_download)
    call self%set(bits%ignore)
    call trace%end()
  END SUBROUTINE finish

  !===============================================================================
  LOGICAL FUNCTION is_finished (self)
    class(m3d_patch_t):: self
    !-----------------------------------------------------------------------------
    is_finished = self%time >= io%end_time
  END FUNCTION is_finished

  !===============================================================================
  SUBROUTINE m3d_nbor_relations(self, another, needed, needs_me, download)
    class(m3d_patch_t) :: self
    class(task_t), pointer :: another         ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                             !    does it need me, should I download?
    class(task_t), pointer :: other
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%m3d_nbor_relations')

    if (.not.associated(another)) then
      write (stderr,*) mpi%rank, 'another not associated'
      call io%abort ('task_t another not associated')
    end if

    needed=.true. ; needs_me=.true. ; download=.true.

    if (another%is_set (bits%frozen)) then
      needs_me = .false.
      needed = .false.
    endif

    if (self%rank /= mpi%rank) then
      !---------------------------------------------------------------------------
      ! virtual tasks don't need guard zones
      !---------------------------------------------------------------------------
      download = .false.
      if (another%rank /= mpi%rank) then
        !---------------------------------------------------------------------------
        ! task is virtual, nbor is virtual
        !---------------------------------------------------------------------------
        needed = .false.
        needs_me = .false.
      endif
    endif

    if (needed .or. needs_me) then
      select type(another)
      class is(m3d_patch_t)
        if (another%kind == self%kind) then
          if (another%patch%id == self%patch%id) then
            download=.false.
          elseif ((another%i_mu /= self%i_mu) .or. (another%i_nu /= self%i_nu)) then
            needed=.false. ; needs_me=.false. ; download=.false.
          endif
        else
          download=.false.
          if (another%patch%id /= self%patch%id) then
            needed=.false. ; needs_me=.false.
          else
            if (self%i_mu /= 0 .and. another%i_mu /= 0) then
              if (another%i_mu /= self%i_mu) then
                needed=.false. ; needs_me=.false. ; download=.false.
              endif
            endif

            if (self%i_nu /= 0 .and. another%i_nu /= 0) then
              if (another%i_nu /= self%i_nu) then
                needed=.false. ; needs_me=.false. ; download=.false.
              endif
            endif
          endif
        endif
      class default
        if (another%id == self%patch%id) then
          needs_me=.false. ; download=.false.
        else
          needed=.false. ; needs_me=.false. ; download=.false.
        endif
      end select
    endif
  
    call trace%end()
  END SUBROUTINE m3d_nbor_relations
  
  !===============================================================================
  SUBROUTINE append_to_task_list(self)
    class (m3d_patch_t), target :: self
    class (m3d_halo_t), pointer :: link
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%add_to_task_list', itimer=itimer)
    if (associated(self%task_list)) then
      if (self%maxiter >= 0) then
        allocate (link)
        call os%mem_register(storage_size(link), 1, label='link')
        call link%init()
        link%task => self
        self%link => link
        !call link%append_to(self%task_list)
        call self%task_list%append_link(self%link)
        call self%task_list%hash_table%add (self%key, self%link)
      endif
    endif
    call trace%end(itimer)
  END SUBROUTINE append_to_task_list
  
  !===============================================================================
  !> For diagnostic output of tasks and nbor lists to data/run/rank_RRRRR.log
  !===============================================================================
  SUBROUTINE task_info (self, unit)
    class(m3d_patch_t):: self
    integer, optional:: unit
    integer:: unit_L
    !-----------------------------------------------------------------------------
    unit_L = shared%option (unit, io_unit%mpi)
    write(unit_L,'(1x,a,2x,a10,i6,4i4,7x,2l1,2x,3f9.3)') &
      'm3d_patch_t%task_info: kind, id, rank, level, nu, mu, BV, pos', &
      self%kind, self%id, self%rank, self%level, self%i_nu, self%i_mu, &
      self%is_set(bits%boundary), self%is_set(bits%virtual), self%position
  END SUBROUTINE task_info
  
  !===============================================================================
  SUBROUTINE init_level (self)
    class(m3d_patch_t):: self
    !-----------------------------------------------------------------------------
    self%level = self%patch%level
    !self%level = 1
  END SUBROUTINE init_level
  
  !===============================================================================
  SUBROUTINE extrapolate_guards (self, target_mem)
    class(m3d_patch_t):: self
    real(4), pointer, contiguous, optional:: target_mem(:,:,:,:)
    !.............................................................................
    ! do nothing
  END SUBROUTINE
  
  !===============================================================================
  !> Check if self (which is an nbor task) is ahead of target (which is the one to
  !> possibly move to the ready queue), using self%dtime*target%grace as the grace
  !> period, since we want to limit the extrapolation in the nbor task to at most
  !> target%grace*self%dtime.
  !===============================================================================
  LOGICAL FUNCTION is_ahead_of (self, target)
    class (m3d_patch_t) :: self
    class (task_t), pointer :: target
    !.............................................................................
    real:: grace
    real(8):: nbtime, nbdtime, tgtime
    logical:: needed, needs_me, download
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%is_ahead_of')
    !---------------------------------------------------------------------------
    ! An nlte_patch that converges finishes, and thereby jumps to io%end_time,
    ! which would make it -- and every other task in its column -- ahead of the
    ! output-spectrum tasks there.  Those must not start on that account: the
    ! long characteristics would then run, and allocate their frequency-resolved
    ! arrays, while the rest of the domain is still iterating.  Withhold them
    ! until every nlte_patch has converged; release_output_tasks() then puts them
    ! back in play, since by that time no nbor of theirs updates any more.
    !---------------------------------------------------------------------------
    if (.not. nlte_finished) then
      select type (target)
      class is (m3d_patch_t)
        if (target%wait_for_nlte) then
          is_ahead_of = .false.
          call trace%end()
          return
        endif
      end select
    endif
    grace = self%grace + shared%extra_grace
                                            call self%lock%set ('is_ahead_of')
    !---------------------------------------------------------------------------
    ! dt(self%dt) is always consistent with t(:), while dtime changes in
    ! courant_condition(), before other attributes change in rotate()
    !---------------------------------------------------------------------------
    nbtime = self%time
    tgtime = target%time
    is_ahead_of = .false.
    if (self%is_set(bits%frozen)) then
      is_ahead_of = .true.
    elseif (nbtime == tgtime) then
      call self%nbor_relations(target, needed, needs_me, download)
      if (needed .eqv. needs_me) then
        select type(target)
        class is (m3d_patch_t)
          if (associated(target%m3d_parent)) then
            is_ahead_of = .not. target%m3d_parent%id == self%id
          else
            is_ahead_of = .true.
          endif
        end select
      endif
    elseif (nbtime > tgtime) then
      is_ahead_of = .true.
    end if
                                            call self%lock%unset ('is_ahead_of')
    call trace%end()
  END FUNCTION is_ahead_of

  !===============================================================================
  !> Put the output-spectrum tasks back in play, once the last nlte_patch has
  !> converged.  This has to be done explicitly: check_ready() is normally called
  !> on a task by its nbors, as they update, but the nbors of these tasks -- the
  !> nlte_patches and compositions of their slabs -- may all have finished by
  !> now, so no update that would check them is left.
  !===============================================================================
  SUBROUTINE release_output_tasks (task_list)
    class(list_t), pointer :: task_list
    class(link_t), pointer :: link
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%release_output_tasks')
    link => task_list%head
    do while (associated(link))
      select type (task => link%task)
      class is (m3d_patch_t)
        if (task%wait_for_nlte) then
          call task_list%check_ready (link, lock=.true.)
        endif
      end select
      link => link%next
    end do
    call trace%end()
  END SUBROUTINE release_output_tasks

  !===============================================================================
  !> The composition of the i-th rotated copy of the atmosphere.  i=1 is the
  !> unrotated one, held in self%composition; the rest resample the cube along
  !> their own ray direction and live in self%rotated_composition.
  !===============================================================================
  FUNCTION comp_ptr (self, i) RESULT (comp)
    class(m3d_patch_t), target :: self
    integer, intent(in) :: i
    class(composition_type), pointer :: comp
    !-----------------------------------------------------------------------------
    if (i == 1) then
      comp => self%composition
    else
      comp => self%rotated_composition(i-1)
    endif
  END FUNCTION comp_ptr

  !===============================================================================
  SUBROUTINE allocate_mem (self, m_in)
    class(m3d_patch_t):: self
    integer, intent(in), optional :: m_in(6)
    integer :: m(6), li(3), ui(3)
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    ! nv=0 marks a memory-less patch (composition, nlte): its state lives in the
    ! composition objects and in atom%data, and mem stays null
    if (self%nv == 0) return
    if (.not. associated(self%mem)) then
      call trace%begin ('m3d_patch_t%allocate_mem', itimer=itimer)
      if (present(m_in)) then
        m = m_in
      else
        m = [self%gn, self%nv, self%nt, self%nw]
        self%compressed = .false.
      endif
      call self%patch_t%allocate_mem(m)
      if (check_NaN) then
        self%mem = tiny(1.0)
      endif
      li = self%li
      ui = self%ui
      if (self%compressed) then
        li(3) = 1
        ui(3) = 1 
      endif
      self%authority_mem => self%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),:,:,:)
      call trace%end (itimer)
    endif
  END SUBROUTINE allocate_mem

  !===============================================================================
  ! Merge the line-list diagnostics of this patch's composition copies (unrotated,
  ! rotated and the sampled table's) into the global composition they were copied
  ! from, which write_linelist_out exports. Called by every rt_patch after its final
  ! pass when save_linelist is set, so every line strength and reach is the maximum
  ! over all patches.
  !-------------------------------------------------------------------------------
  SUBROUTINE collect_line_diagnostics (self)
    class(m3d_patch_t), target :: self
    integer :: i
    !-----------------------------------------------------------------------------
    !$omp critical (line_diagnostics)
    have_line_diagnostics = .true.
    if (associated(self%atoms)) linelist_atoms => self%atoms
    associate (p => self%m3d_parent)
      if (associated(p%composition)) call merge_copy(p%composition)
      if (associated(p%rotated_composition)) then
        do i=1, size(p%rotated_composition)
          call merge_copy(p%rotated_composition(i))
        end do
      endif
      if (associated(p%tab_comp)) call merge_copy(p%tab_comp)
    end associate
    !$omp end critical (line_diagnostics)
  contains
    subroutine merge_copy (copy)
      type(composition_type), target, intent(in) :: copy
      type(composition_type), pointer :: cp
      !---------------------------------------------------------------------------
      cp => copy
      if (associated(cp, composition)) return
      call composition%merge_line_diagnostics(copy)
    end subroutine merge_copy
  END SUBROUTINE collect_line_diagnostics

  !===============================================================================
  ! Write linelist_out.txt in the same block format as the input linelists, via the
  ! writers composition_type%write_linelists and linelist_mod%print_ll_atom, from the global
  ! composition and the model atoms that the final-pass patches merged their copies
  ! into. Called once by the main program after the task list has run, when every
  ! merge is complete; does nothing unless a patch collected diagnostics
  ! (save_linelist). This driver only walks the composition and the model atoms -- it
  ! stays here because linelist_mod sits below composition_mod/atom_mod in the module
  ! hierarchy and cannot see their types.
  !-------------------------------------------------------------------------------
  SUBROUTINE write_linelist_out ()
    integer :: lunit
    ! ----------------------------------------------------------------------------
    if (.not. have_line_diagnostics) return
    open(newunit=lunit, file=trim(io%run_dir)//'linelist_out.txt', form='formatted')
    call composition%write_linelists(lunit)
    call print_atoms()
    close(lunit)

  contains

    !-----------------------------------------------------------------------------
    ! model-atom SoA linelists; gather the per-stage species ids from the
    ! composition, then let linelist_mod write the (stage, isotope) blocks
    subroutine print_atoms ()
      class(atom_type), pointer :: atom
      class(element_type), pointer :: el
      class(molec_type), pointer :: ml
      character(len=20), allocatable :: ids(:)
      integer :: ia, s, s0, s1
      !---------------------------------------------------------------------------
      if (.not. associated(linelist_atoms)) return
      do ia=1, size(linelist_atoms)
        atom => linelist_atoms(ia)
        if (.not. associated(atom%trans)) cycle
        associate(sll => atom%trans%ll)
          if (sll%nline == 0) cycle
          nullify(el, ml)
          call composition%element_by_id(atom%id, el, verbose=.false.)
          if (.not. associated(el)) call composition%molecule_by_id(atom%id, ml, verbose=.false.)

          s0 = minval(sll%ipart(1:sll%nline))
          s1 = maxval(sll%ipart(1:sll%nline))
          if (allocated(ids)) deallocate(ids)
          allocate(ids(s0:s1))
          do s=s0, s1
            ids(s) = atom%id
            if (associated(el)) then
              if (s >= lbound(el%ions,1) .and. s <= ubound(el%ions,1)) ids(s) = el%ions(s)%id
            elseif (associated(ml)) then
              ids(s) = ml%id
            endif
          end do
          call print_ll_atom(lunit, sll, ids, s0)
        end associate
      end do
    end subroutine print_atoms
  END SUBROUTINE write_linelist_out

  !===============================================================================
  SUBROUTINE check_mem_NaN (self)
    class(m3d_patch_t):: self
    integer :: i,j,k,v,t,w
    integer, save :: itimer
    character(len=128) :: msg
    !-----------------------------------------------------------------------------
    if (.not. check_NaN) return
    call trace%begin ('m3d_patch_t%check_mem_NaN', itimer=itimer)
    if (associated(self%mem)) then
      do w=lbound(self%mem, 6), ubound(self%mem, 6)
      do t=lbound(self%mem, 5), ubound(self%mem, 5)
      do v=lbound(self%mem, 4), ubound(self%mem, 4)
      do k=lbound(self%mem, 3), ubound(self%mem, 3)
      do j=lbound(self%mem, 2), ubound(self%mem, 2)
      do i=lbound(self%mem, 1), ubound(self%mem, 1)
        if (self%mem(i,j,k,v,t,w) /= self%mem(i,j,k,v,t,w)) then
          write(msg, '(A,A,3(I4),I5,I3,I3)') trim(self%m3d_patch_info()), '  NaN', i,j,k,v,t,w
          call console_error (msg)
          error stop
        elseif (self%mem(i,j,k,v,t,w) > inf4) then
          write(msg, '(A,A,3(I4),I5,I3,I3)') trim(self%m3d_patch_info()), '  Infinity', i,j,k,v,t,w
          call console_error (msg)
          error stop
        !elseif (self%unsigned(w) .and. self%mem(i,j,k,v,t,w) < 0) then
        !  write(stdout, '(A,A,3(I4),I5,I3,I3)') self%m3d_patch_info(), '  negative', i,j,k,v,t,w
        !  stop
        endif
      end do
      end do
      end do
      end do
      end do
      end do
    endif
    call trace%end(itimer)
  END SUBROUTINE check_mem_NaN

  !===============================================================================
  SUBROUTINE set_place (self, other, place)
    class(m3d_patch_t) :: self
    class(task_t):: other
    integer :: place
    !-----------------------------------------------------------------------------
    other%place = place
    if (self%id /= other%id .and. console_on(ch_m3d, 3)) then
      write(console_unit,'(A, A, A, A, I4, A, I4)') trim(self%m3d_patch_info2()), ' set place of ', &
      trim(other%kind), ' ', other%id, ' to ', other%place
    endif
  END SUBROUTINE set_place

  !===============================================================================
  SUBROUTINE print_to_stdout (self, arr, text)
    class(m3d_patch_t) :: self
    real(8), dimension(:), optional :: arr
    character(len=*), optional :: text
    !-----------------------------------------------------------------------------
    !$omp critical (cr_print_to_stdout)
    if (present(arr)) then
      if (size(arr)==0) then
        write(console_unit,'(A, A)') trim(patch_info()), 'invalid index!'
      else
        write(console_unit,'(A, 6(ES15.7))') trim(patch_info()), arr
      endif
    elseif (present(text)) then
      write(console_unit,'(A, A, A)') trim(patch_info()), ' ', trim(text)
    else
      write(console_unit,'(A)') trim(patch_info())
    endif
    flush(console_unit)
    !$omp end critical (cr_print_to_stdout)
  contains
    function patch_info() result(text)
      character(len=32) :: text
      !-----------------------------------------------------------------------------
      write(text, '(A12, I4)')  trim(self%kind), self%key(1)
      if (self%i_nu /= 0) write(text, '(A, I4)') trim(text), self%i_nu
      if (self%i_mu /= 0) write(text, '(A, I4)') trim(text), self%i_mu
    end function
  END SUBROUTINE print_to_stdout

  !===============================================================================
  SUBROUTINE print_post_update_info (self)
    class(m3d_patch_t) :: self
    !-----------------------------------------------------------------------------
    if (console_on(ch_m3d, 2, all_ranks=.true.)) then
      call self%print_to_stdout()
    endif
  END SUBROUTINE print_post_update_info

  !===============================================================================
  FUNCTION check_in_bounds (self) result (cib)
    class(m3d_patch_t) :: self
    logical :: cib
    !-----------------------------------------------------------------------------
    cib = all(self%pcheck > 0) .and. all(self%pcheck <= self%n)
  END FUNCTION check_in_bounds

  !===============================================================================
  SUBROUTINE print_execution_history (self, status)
    class(task_t) :: self
    character(len=3) :: status
    !-----------------------------------------------------------------------------
    !$omp critical (cr_print_execution_history)
    write(io_unit%dispatcher, &
    '(F17.8, 2(I3, 1x), 1x, A3, 1x, A12, 1x, I8, 1x, 3(I3,1x), F17.8)') &
    wallclock(), mpi%rank, omp%thread, status, trim(self%kind), self%id, self%key(1:3), self%time
    flush(io_unit%dispatcher)
    !$omp end critical (cr_print_execution_history)
  END SUBROUTINE print_execution_history

  !===============================================================================
  !> Return a unique hash key for the patch
  !===============================================================================
  FUNCTION hash_key (self, position, size, level, box, origin) RESULT (key)
    class(m3d_patch_t):: self
    real(8), dimension(3), optional:: position, origin, size, box
    integer, optional:: level
    integer:: key(4)
    !-----------------------------------------------------------------------------
    key(1) = self%ipos(3) + 1 
    key(2) = self%i_nu
    key(3) = self%i_mu
    key(4) = hash_string(self%kind)
  END FUNCTION hash_key

  !===============================================================================
  !> Fast method to generate nbor lists in cases with very long task lists
  !===============================================================================
  SUBROUTINE init_nbors_hash (self, nbors)
    class(m3d_patch_t):: self
    class(nbor_t), pointer:: nbors
    !.............................................................................
    class(task_t), pointer:: task
    class(nbor_t), pointer:: nbor
    class(key_link_t), pointer :: klink
    integer, save:: itimer=0
    integer :: i, i_nu, key(4)
    !-----------------------------------------------------------------------------
    key(1) = self%key(1)
    key(2:3) = 0
    key(4) = hash_string(self%patch%kind)
    call self%check_nbor(key, nbors)

    if (associated(self%m3d_parent)) then
      key(4) = hash_string(self%m3d_parent%kind)
      call self%check_nbor(key, nbors)
    endif

    klink => self%key_list
    do while (associated(klink))
      call self%check_nbor(klink%key, nbors)
      klink => klink%next
    end do
  END SUBROUTINE init_nbors_hash

  !===============================================================================
  SUBROUTINE check_nbor (self, key, nbors)
    class(m3d_patch_t):: self
    class(nbor_t), pointer:: nbors
    integer:: key(4)
    !-----------------------------------------------------------------------------
    call add_nbor_by_key (self, key, nbors)
  END SUBROUTINE check_nbor

  !===============================================================================
  SUBROUTINE rotate_cube4 (self, vec, cube, do_log)
    class(m3d_patch_t) :: self
    real(8), intent(in) :: vec(3)
    real(4), intent(inout), dimension(:,:,:) :: cube
    logical, intent(in), optional :: do_log
    real(4), dimension(:,:), allocatable :: shifted
    real(4), dimension(2) :: delta, deltak
    real(4) :: bottom, bot2pivot, height, f1, f2, a, b, c, d
    integer :: x, y, z, up1, up2, x1, x2, y1, y2, i, nx, ny, nz
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin('m3d_patch_t%rotate_cube4', itimer=itimer)
    delta = vec(:2) / (vec(3) * self%ds(:2)) ! shift in units of grid cells
  
    ! bot2pivot is the distance between the bottom of the patch and the pivot plane
    bottom = self%position(3) - 0.5*self%size(3)
    bot2pivot = pivot_z - bottom
  
    allocate(shifted, mold=cube(:,:,1))
    nx = size(cube, 1)
    ny = size(cube, 2)
    nz = size(cube, 3)

    if (option(do_log)) then
      cube = log(cube)
    endif
  
    do z=1, nz
      associate(plane => cube(:,:,z))
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
        do y=1, ny
          y1 = modulo(y+up2-1, ny)+1
          y2 = modulo(y+up2,   ny)+1
      
          do x=1, nx
            x1 = modulo(x+up1-1, nx)+1
            x2 = modulo(x+up1,   nx)+1
    
            shifted(x,y) = a*plane(x1,y1)
            shifted(x,y) = b*plane(x2,y1) + shifted(x,y)
            shifted(x,y) = c*plane(x1,y2) + shifted(x,y)
            shifted(x,y) = d*plane(x2,y2) + shifted(x,y)
          end do
        end do
        plane = shifted
      end associate
    end do
    deallocate(shifted)
    if (option(do_log)) then
      cube = exp(cube)
    endif
    call trace%end(itimer)
  END SUBROUTINE rotate_cube4

  !===============================================================================
  SUBROUTINE add_to_key_list (self, other)
    class(m3d_patch_t), intent(inout) :: self, other
    class(key_link_t), pointer :: next, link, prev
    ! ----------------------------------------------------------------------------
    nullify(prev)
    next => self%key_list
    do while(associated(next))
      prev => next
      next => next%next
    end do
    allocate(link)
    link%key = other%key
    if (associated(prev)) then
      prev%next => link
    else
      self%key_list => link
    endif
  END SUBROUTINE add_to_key_list

  !===============================================================================
  FUNCTION queue_before (self, other) result(queue)
    class(m3d_patch_t):: self
    class(task_t):: other
    logical :: queue
    !-----------------------------------------------------------------------------
    if (self%prepend) then
      queue = other%time >= self%time
      !self%prepend = .false.
    else
      queue = self%task_t%queue_before(other)
    endif
  END FUNCTION queue_before

  !===============================================================================
  SUBROUTINE sort_nbors (self)
    class(m3d_patch_t) :: self
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: nbtask
    integer:: n_nbors
    !-----------------------------------------------------------------------------
    nbor => self%link%nbor
    do while(associated(nbor))
      nbtask => nbor%link%task
      nbtask%istep = nbtask%level
      nbtask%level = -count(nbtask%key == self%key)
      nbor => nbor%next
    end do  

    nbor => self%link%nbor
    call sort_nbors_by_level(self%link, nbor, nbor)

    do while(associated(nbor))
      nbtask => nbor%link%task
      nbtask%level = nbtask%istep
      nbtask%istep = 0
      nbor => nbor%next
    end do  
  END SUBROUTINE sort_nbors

END MODULE m3d_patch_mod
