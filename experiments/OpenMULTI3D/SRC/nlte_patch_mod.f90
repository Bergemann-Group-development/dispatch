MODULE nlte_patch_mod
  USE os_mod
  USE gpatch_mod
  USE io_mod
  USE scaling_mod
  USE trace_mod
  USE bits_mod
  USE task_mod, only: task_t
  USE link_mod, only: nbor_t
  USE handyfuncs_mod
  USE atom_mod, only: atom_type, atomic_spectrum_type, use_classic_col, kramer_bf, &
    red_cut, blue_cut, logf_cut, atom_data_t, atom_trans_t, &
    nq_fac, nq_min, nq_max, disable_bf, disable_bb, level_partf, &
    fmin, cbbesc, cbfesc, cbbhsc, cbfhsc, colscale, ch0_hmip, &
    classic_ce, classic_ci, classic_ch, classic_chi, upsilon_floor, &
    ng_step, ng_start, ng_global, save_rate_matrix, legacy_freql, sor_w, line_grid_file
  USE trans_mod, only: line_type, cont_type, transition_type
  USE const_mod, only: CC
  USE ion_mod, only: ion_type, particle_type, debye_shield, massive_type
  USE element_mod, only: element_type
  USE molecule_mod, only: molec_type
  USE line_types_mod, only: linelist_type
  USE composition_mod, only: composition_type, departure_pointer_t, read_departures
  USE m3d_patch_mod, only: m3d_patch_t, composition, key_link_t, &
    convlim, convmax, &
    decouple_continuum, overlapping_bb, overlapping_bf, converged, premature_finish, &
    rotate_atmos, rank_lvlerr, overlapping_bb_in_output, n_windows, n_windows_out, detailed_absmet, &
    atmos_file, flux_convergence, flux_conv_signal, global_flux_err, global_flux_pos, global_flux_lambda, &
    flux_final_counter, nlte_finished, nlte_last_change, nlte_last_iter, n_convmax_finished, &
    n_nlte_patches, release_output_tasks, &
    save_populations, bound_free_in_output, save_rates, save_every, stop_requested, stop_exit, &
    save_net_rates
  USE rt_params_mod, only: rt_out
  USE bftrans_mod, only: bftrans_type
  USE osmet_mod, only: osmet_type, osmet
  USE m3d_debug_mod
  USE shared_mod
  USE console_mod
  USE, intrinsic :: ieee_arithmetic
  implicit none

  !-----------------------------------------------------------------------------

  type, public, extends(m3d_patch_t) :: nlte_patch_t
    type(atom_type), dimension(:), pointer :: output_atoms=>null(), bb_atoms=>null(), bf_atoms=>null()
    class(atomic_spectrum_type), pointer :: output_sp=>null()
    class(atomic_spectrum_type), pointer :: bf_sp=>null()
    class(atomic_spectrum_type), pointer :: bb_sp=>null()
    real(8) :: maxerr=0, global=0
    real(8), dimension(:), allocatable :: lvlerr
    real(8) :: local
    integer :: natom=0
    logical :: update_skipped=.false., atom_allocated=.false., molec
    logical :: pop_meta_written=.false.   ! this patch has its rows in pop_<atom>_meta.txt
  contains
    procedure :: pre_init => nlte_patch_pre_init
    procedure :: execute
    procedure :: estimate_atom_memory
    procedure :: finish
    procedure :: project_populations
    procedure :: nbor_relations => nlte_patch_nbor_relations
    procedure :: pre_init_atom
    procedure :: pre_init_spectrum
    procedure :: init_departures
    procedure :: write_populations
    procedure :: check_stop_request
    procedure :: get_osmet_nu
    procedure :: exclude_double_continua
    procedure :: exclude_double_lines
    procedure :: print_post_update_info
    procedure :: check_in_bounds
    procedure :: init_nbors_hash
    procedure :: save_abundance_to_file
    procedure :: save_convergence_report
    !procedure :: load_mem
    procedure :: check_mem_NaN
  end type nlte_patch_t

  integer, save :: epos(3), elvl, eatom=1
  !-----------------------------------------------------------------------------
  ! Population convergence is judged per iteration (see execute).  sweep_iter is
  ! the iteration now under way, sweep_active the number of nlte_patches that
  ! update in it, sweep_updates how many have, and sweep_maxerr the largest
  ! change among them.  last_maxerr is the largest change of the previous,
  ! complete iteration.  n_finished counts the finished nlte_patches.
  !-----------------------------------------------------------------------------
  integer, save :: sweep_iter=-1, sweep_active=0, sweep_updates=0, n_finished=0
  real(8), save :: sweep_maxerr=0d0, last_maxerr=0d0
  logical, save :: exclude_from_line_list=.true.
  logical, save :: stale_stop_checked=.false.
  logical, public, save :: external_sp=.false.
  ! The output pass takes the lines of the atom's species from the line list, with the
  ! atom's departure coefficients, instead of from the atom's own transitions:
  ! output_from_line_list (the default) asks for it, output_grid_and_list says that the
  ! run has both an explicit output grid and a line list (set by composition_pre_init),
  ! and line_list_precedence is the decision.  Without an explicit grid the output
  ! spectrum is built from the atom's own lines and these have to stay.
  logical, save :: output_from_line_list=.true.
  logical, public, save :: output_grid_and_list=.false., line_list_precedence=.false.
  !-----------------------------------------------------------------------------

contains

  !===============================================================================
  SUBROUTINE nlte_patch_pre_init (self)
    class(nlte_patch_t), target:: self
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin('nlte_patch_t%nlte_patch_pre_init', itimer=itimer)
    self%kind = 'nlte_patch'
    call self%m3d_patch_t%pre_init()

    if (self%atom_allocated) then
      ! we do not need a memory array for this patch,
      ! as the atom data is allocated separately for each atom
      call self%set_nv_and_nw(0, 0)
      self%static_mem=.true.
    elseif (.not. read_departures) then
      self%maxiter = -1 ! prevents addition to task list
      call trace%end(itimer)
      return
    endif

    call self%post_pre_init()
    call trace%end(itimer)
  END SUBROUTINE nlte_patch_pre_init

  !===============================================================================
  SUBROUTINE pre_init_atom (self)
    class(nlte_patch_t) :: self
    type(atom_type), save, dimension(:), allocatable :: atoms
    type(atom_trans_t), save, dimension(:), pointer :: dc_trans, out_trans
    class(element_type), pointer :: element
    class(molec_type), pointer :: molecule
    class(ion_type), pointer :: ion
    class(massive_type), pointer :: species
    character(len=256), save, dimension(:), allocatable :: atom_file
    character(len=512), save :: atom_root=''
    real(8), dimension(:), allocatable :: abundance
    real(8), save :: qnorm=10.0
    integer, save :: itimer=0, natom=0, n_in_lines=0
    integer, save, allocatable :: in_lines(:)
    logical, save :: first_time=.true., exclude_trace_cont=.true., use_atom_abnd=.false.
    logical, save :: relative_abundance = .false.
    logical :: pandas_format
    integer :: i, j, iostat, nlvl
    namelist /atom_params/ atom_root, atom_file, qnorm, exclude_trace_cont, level_partf, &
      use_atom_abnd, logf_cut, red_cut, blue_cut, n_in_lines, abundance, relative_abundance, &
      exclude_from_line_list, output_from_line_list, use_classic_col, kramer_bf, overlapping_bf, &
      classic_ce, classic_ci, classic_ch, classic_chi, upsilon_floor, &
      overlapping_bb, detailed_absmet, &
      nq_fac, nq_min, nq_max, cbbesc, cbfesc, cbbhsc, cbfhsc, colscale, &
      disable_bf, disable_bb, fmin, &
      ch0_hmip, legacy_freql, line_grid_file
    namelist /line_mask/ in_lines
    !-----------------------------------------------------------------------------
    call trace%begin('nlte_patch_t%pre_init_atom', itimer=itimer)
    if (first_time) then
      first_time = .false.
      allocate(atom_file(64))
      allocate(abundance(64), source=-99.9d0)
      atom_file = ''

      rewind (io%input)
      read (io%input, atom_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, atom_params)
      end if
      call io%namelist_warning ('atom_params', iostat)

      line_list_precedence = output_from_line_list .and. output_grid_and_list

      natom = count(atom_file /= '')
      ! atom_root is prepended to every atom_file
      if (atom_root /= '') then
        do i = 1, natom
          atom_file(i) = trim(atom_root)//'/'//trim(atom_file(i))
        end do
      endif
      if (natom == 0) then
        if (console_on(ch_m3d, 2)) write(console_unit,'(a)') 'no atom specified'
        call trace%end(itimer)
        return
      else
        ! this is the atom to be read
        allocate(atoms(natom))

        allocate(out_trans(natom))
        if (decouple_continuum) then
          allocate(self%bb_atoms(natom))
          allocate(self%bf_atoms(natom))
          allocate(dc_trans(natom))
        endif

        atoms(atomcheck)%check = .true.
        call truncate_array(atom_file, natom)
        call truncate_array(abundance, natom)
      endif

      if (n_in_lines > 0) then
        allocate(in_lines(n_in_lines))
        rewind (io%input)
        read (io%input, line_mask, iostat=iostat)
        if (iostat > 0) then
          rewind (io%input)
          read (io%input, line_mask)
        end if
        write (io%namelist, line_mask)
      endif

      do i=1, natom
        associate(atom_file=>atom_file(i), atom=>atoms(i))
          atom%name = 'input'
          j = rindex(atom_file, '/')+1
          pandas_format = .not. startswith(atom_file(j:len_trim(atom_file)), 'atom.')
          if (console_on(ch_m3d, 1)) write(console_unit,'(A, A)') 'Reading atom file: ', trim(atom_file)
          if (pandas_format) then
            level_partf = .false.
            call atom%read_pandas(atom_file, qnorm, in_lines=in_lines)
            call atom%set_colliders(composition)
          else
            call atom%read_classic(atom_file, qnorm, in_lines=in_lines, read_coll=.true.)
          endif

          call composition%element_by_id(atom%id(1:3), element, verbose=.false.)
          if (associated(element)) then
            atom%species => element
            element%ions%in_opac = .true.
            atom%id = element%id
            if (trim(element%id) == 'H') then
              element%ions(1)%nk = max(element%ions(1)%nk, atom%nlevel)
            endif
          else
            call composition%molecule_by_id(atom%id(1:3), molecule, verbose=.false.)
            if (associated(molecule)) then
              atom%species => molecule
              molecule%in_opac = .true.
              atom%id = molecule%id
              do j=1, molecule%nions
                ion => molecule%ions(j)%p
                ion%in_opac = .true.
              end do
            endif
          endif

          call atom%prepare_damping(atom%species)

          if (abundance(i) > -99) atom%abnd = abundance(i)
          if (use_atom_abnd .or. abundance(i) > -99) then
            call composition%overwrite_abundance(atom%species, atom%abnd, atom%awgt, relative_abundance)
          else
            atom%abnd = atom%abnd + composition%FeH
          endif

          if (exclude_trace_cont) then
            call self%exclude_double_continua(atom)
          endif

          ! here we create copies of the line and cont arrays,
          ! which will have altered ired and iblue for example
          ! they will be pointed at later
          call out_trans(i)%allocate(bb=atom%trans%line, bf=atom%trans%cont)
          if (line_list_precedence) then
            ! the output pass leaves the atom's lines to the line list; the atom's
            ! own transitions stay active in the iterations
            do j=1, size(out_trans(i)%line)
              call out_trans(i)%line(j)%deactivate()
            end do
          endif
          if (decouple_continuum) then
            call dc_trans(i)%allocate(bb=atom%trans%line, bf=atom%trans%cont)
          endif

          ! global SoA copy of the bb transitions, exported in linelist_out.txt;
          ! build on the output transitions (out_trans), since those are what the
          ! output patch's output_atoms point at and what write_linelist_out walks
          call atom%build_ll(out_trans(i))

          if (pandas_format) then
            call os%copy_file(trim(atom_file), trim(io%run_dir)//'atoms/'//trim(atom%id)//'.txt')
          else
            call os%copy_file(trim(atom_file), trim(io%run_dir)//'atoms/'//'atom.'//trim(atom%id))
          endif

          if (ng_global .and. ng_step >= 0) then
            ! One shared Ng state per species: copying the atom copies the
            ! pointer, not the target, so every patch accumulates into this one
            ! object and they all extrapolate with the same coefficients.
            allocate(atom%ng)
          endif
        end associate
      end do
      if (natom > 0) then
        nlvl = sum(atoms%nlevel, 1)
        allocate(rank_lvlerr(nlvl), source=0.d0)
      endif
      write (io%namelist, atom_params)
    endif

    if (natom == 0) then
      call trace%end(itimer)
      return
    endif

    ! this is the patch atom, which contains local data
    allocate(self%atoms(natom))

    ! this is the output atom, which contains local data
    ! shared with self%atoms, but with different transitions
    allocate(self%output_atoms(natom))

    ! here we create a copy of atom, but the transitions will
    ! still be linked to the original atom
    self%atoms = atoms
    self%atoms%name = 'main'

    if (decouple_continuum) then
      allocate(self%bb_atoms(natom))
      allocate(self%bf_atoms(natom))
    endif

    do i=1, natom
      ! this pointer will be referenced by other atoms with same local data
      allocate(self%atoms(i)%data)
      if (rotate_atmos) then
        allocate(self%atoms(i)%rotated_data(size(rt_out%vec, 2)))
      endif

      ! here we create a copy of the patch atom, which will share the same
      ! atom%data array, but we link it to the output transitions
      self%output_atoms(i) = self%atoms(i)
      self%output_atoms(i)%trans => out_trans(i)
      self%output_atoms(i)%name = 'output'

      if (decouple_continuum) then
        self%bb_atoms(i) = self%atoms(i)
        self%bf_atoms(i) = self%atoms(i)
        allocate(self%bb_atoms(i)%trans) ! this will detach the original pointer and
        allocate(self%bf_atoms(i)%trans) ! create a new one, to be pointed elsewhere
        self%bb_atoms(i)%trans%line => dc_trans(i)%line
        self%bf_atoms(i)%trans%cont => dc_trans(i)%cont
        self%bb_atoms(i)%name = 'bound-bound'
        self%bf_atoms(i)%name = 'bound-free'
      endif
    end do
    self%atom_allocated = .true.
    self%natom = natom
    call trace%end(itimer)
  END SUBROUTINE pre_init_atom

  !===============================================================================
  SUBROUTINE pre_init_spectrum (self, output_sp)
    class(nlte_patch_t) :: self
    type(atomic_spectrum_type), target :: output_sp
    type(atomic_spectrum_type), save, target :: sp, bb_sp, bf_sp
    logical, save :: first_time=.true.
    logical :: save_sp, saved_already
    integer :: i
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('nlte_patch_t%pre_init_spectrum', itimer=itimer)
    if (first_time) then
      first_time = .false.
      saved_already = .false.
      external_sp = output_sp%nnu > 0
      save_sp = external_sp .or. (associated(self%atoms) .and. size(rt_out%vec,2)>0)
      if (associated(self%atoms)) then
        if (decouple_continuum) then
          if (detailed_absmet .and. overlapping_bf) call self%get_osmet_nu(bf_sp)
          call setup_spectrum(bf_sp, 'bf', self%bf_atoms, bb=.false., bf=.true.,  sort=overlapping_bf, set_active=.true.)
          call setup_spectrum(bb_sp, 'bb', self%bb_atoms, bb=.true.,  bf=.false., sort=overlapping_bb, set_active=.true.)
        else
          if (detailed_absmet) call self%get_osmet_nu(sp)
          call setup_spectrum(sp, 'all', self%atoms, bb=.true., bf=.true., sort=.true., set_active=.true.)
        endif

        if (.not. external_sp) then
          call setup_spectrum(output_sp, 'out', self%output_atoms, bb=.true., bf=bound_free_in_output, sort=overlapping_bb_in_output, set_active=.false., crop=.true.)
          saved_already = .true.
        else
          do i=1, size(self%output_atoms)
            call self%output_atoms(i)%trans%set_ired_iblue(output_sp, bb=.true., bf=.false., set_active=.false.)
          end do
        endif
      endif
      if (save_sp .and. .not. saved_already) then
        call output_sp%save('out')
      endif
    endif

    self%output_sp => output_sp

    if (decouple_continuum) then
      self%bf_sp => bf_sp
      self%bb_sp => bb_sp
    else
      self%sp => sp
    endif

    call trace%end(itimer)
  contains
    !===============================================================================
    subroutine setup_spectrum(sp, name, atoms, bb, bf, sort, set_active, crop)
      class(atomic_spectrum_type), target :: sp
      character(len=*), intent(in) :: name
      class(atom_type), intent(in), dimension(:) :: atoms
      integer, dimension(:,:), pointer, contiguous :: nac, nal
      logical, intent(in) :: bb, bf, sort, set_active
      logical, intent(in), optional :: crop
      integer :: i
      character(len=160) :: text, detail
      character(len=24) :: what
      !-----------------------------------------------------------------------------
      sp%name = name
      select case (name)
      case ('bf');  what = 'bound-free'
      case ('bb');  what = 'bound-bound'
      case default; what = 'bound-bound and bound-free'
      end select
      ! allocate list of active continua and lines
      do i=1, size(atoms)
        call atoms(i)%trans%append_transitions(sp, bb=bb, bf=bf)
      end do
      ! a file grid is shared by all lines, so it has to be located by wavelength
      if (sort .or. (bb .and. line_grid_file /= '')) call sp%sort()
      ! restrict the output grid to a bare aa_red/aa_blue window (crop_output); a no-op
      ! otherwise. Must precede set_ired_iblue, which is recomputed against this grid.
      if (option(crop)) call sp%crop()
      call os%allocate(nac, [sp%nnu, size(atoms)], zero=.true.)
      call os%allocate(nal, [sp%nnu, size(atoms)], zero=.true.)
      call os%allocate(sp%na, sp%nnu, zero=.true.)
      do i=1, size(atoms)
        sp%nac => nac(:,i)
        sp%nal => nal(:,i)
        call atoms(i)%trans%set_ired_iblue(sp, bb=bb, bf=bf, set_active=set_active)
      end do
      sp%na = sum(min(nal+nac, 1), dim=2)
      if (save_sp) then
        call sp%save(name)
      endif
      if (set_active .and. console_on(ch_m3d, 1)) then
        write(text,'(a,a,a,I0,a)') 'NLTE frequency grid (', trim(what), '): ', sp%nnu, ' frequencies'
        ! how many transition-frequency pairs the rate sums run over, and the largest
        ! number of transitions sharing a frequency
        if (console_on(ch_m3d, 2)) then
          write(detail,'(a,I0,a,I0,a)') ', ', sum(nac+nal), ' active transition-frequency pairs, at most ', &
            maxval(sum(nac+nal, dim=2)), ' at one frequency'
          text = trim(text) // trim(detail)
        endif
        write(console_unit,'(a)') trim(text)
      endif
      call os%deallocate(nac)
      call os%deallocate(nal)
    end subroutine setup_spectrum
  END SUBROUTINE pre_init_spectrum

  !===============================================================================
  SUBROUTINE execute (self)
    class(nlte_patch_t) :: self
    class(atom_type), pointer :: atom
    class(atom_data_t), pointer :: adata
    class(composition_type), pointer :: comp
    class(line_type), pointer :: line
    class(cont_type), pointer :: cont
    class(element_type), pointer :: element
    class(molec_type), pointer :: molec
    real(4), dimension(:,:,:), allocatable :: cij, cji
    real(8), dimension(:,:), allocatable :: w, c
    integer :: i, j, k, p(3), m(3), ng_epos(4), count, lunit
    character(len=32) :: rates_tag
    real(8) :: ng_emax
    logical :: do_ng, locally_converged, found, last
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin (trim(self%kind)//'%update', itimer=itimer)
    comp => self%m3d_parent%composition
    self%maxerr = 0
    locally_converged = premature_finish
    m = self%n
    ! ng_step < 0 disables Ng; 0 selects the automatic trigger; > 0 keeps the
    ! fixed ng_start/ng_step schedule.  The scheduling arithmetic now lives in
    ! ng_acceleration, which needs the raw iteration number either way.
    do_ng = ng_step >= 0

    if (save_rates>0) then
      allocate(cij(m(1),m(2),m(3)), source=0.0)
      allocate(cji(m(1),m(2),m(3)), source=0.0)
    endif

    do k=1, self%natom
      atom => self%atoms(k)
      adata => atom%data

      if (self%iter == 0) then
        call atom%init_atom_data(comp, adata, self%maxiter>0, found)
        if (.not. found) call self%save_abundance_to_file(atom)
      elseif (.not. adata%converged) then
        if (save_rate_matrix) then
          allocate(w(atom%nlevel, atom%nlevel), source=0.d0)
          allocate(c(atom%nlevel, atom%nlevel), source=0.d0)
        endif

        call atom%normalize_profiles(self%iter)
        diag_zs = self%zs
        call atom%update_populations(comp, self%local, self%pcheck, cij, cji, w, c)

        if (save_net_rates) then
          write(rates_tag,'(a,"_",i5.5,".bin")') trim(atom%id), self%zs
          call atom%write_rates(comp, trim(io%run_dir)//'save/rates_'//trim(rates_tag), &
                                trim(io%run_dir)//'save/colrates_'//trim(rates_tag), self%iter)
        endif

        if (save_rates==k) then
          call self%save3(real(adata%gamma_ij(:,:,:,linecheck)), label='rij')
          call self%save3(real(adata%gamma_ji(:,:,:,linecheck)), label='rji')
          call self%save3(cij, label='cij')
          call self%save3(cji, label='cji')
        endif

        if (save_rate_matrix) then
          if (self%iter==1 .and. w(1,1)/=0.d0) then
            call self%save2(real(w), label=trim(atom%id)//'_w')
            call self%save2(real(c), label=trim(atom%id)//'_c')
          endif
        endif

        if (adata%maxerr /= adata%maxerr) then
          call shared%abort('NaN in convergence error')
        endif

        adata%converged = abs(adata%maxerr) < convmax
        if (abs(self%maxerr) < abs(adata%maxerr)) then
          self%maxerr = adata%maxerr
          epos = adata%epos + [0,0,self%zs-1]
          elvl = adata%elvl
          eatom = k
        endif

        if (do_ng .and. .not. adata%converged) then
          ! self%global, not the patch's own maxerr: the extrapolation is one
          ! global step, so the headroom test that gates it has to see the run's
          ! convergence measure rather than whichever patch happened to decide.
          ! It lags by an iteration, which is immaterial for a coarse "is there
          ! anything left to gain" test.
          call atom%ng_acceleration(self%iter, self%global, convlim, ng_emax, ng_epos)
          if (ng_emax /= 0 .and. console_on(ch_m3d, 2)) then
            ng_epos = ng_epos + [0,0,0,self%zs-1]
            write(console_unit,'(A, ES14.7, 4(I6))') 'NG acceleration ', ng_emax, ng_epos
          endif
        endif
      endif
      
      locally_converged = locally_converged .and. adata%converged

      if (self%print_update_info .and. k==atomcheck) then
        p = self%pcheck
        if (allocated(adata%b)) then
          memcheck(1) = adata%b(p(1),p(2),p(3),lvlcheck)
        else
          memcheck(1) = 1
        endif
        if (allocated(adata%gamma_ij)) then
          if (linecheck<=atom%nrad) then
            memcheck(2) = adata%gamma_ij(p(1),p(2),p(3),linecheck)
            memcheck(3) = adata%gamma_ji(p(1),p(2),p(3),linecheck)
            if (linecheck<=atom%nline) then
              memcheck(4) = 1
              if (allocated(adata%wphi(linecheck)%w)) memcheck(4) = adata%wphi(linecheck)%w(p(1),p(2),p(3))
            endif
          endif
        else
          memcheck(2) = 10.0**(atom%abnd-12.0)
          memcheck(3) = comp%toth(p(1),p(2),p(3))
        endif
      endif

      call os%deallocate(adata%gamma_ij, label='atom%gamma_ij')
      call os%deallocate(adata%gamma_ji, label='atom%gamma_ji')
    end do

    !---------------------------------------------------------------------------
    ! Population convergence.  Every nlte_patch of an iteration updates before
    ! any RT task of that iteration runs, and those finish before the next
    ! iteration starts, so the first update of an iteration marks the previous
    ! one as complete.  Each patch compares the largest change of that complete
    ! iteration with convlim; they all see the same value, whatever the order
    ! they run in.  Once it is below convlim, every patch still performs this
    ! iteration's solve, against the radiation field of the converged
    ! populations, and the last one to do so sets `converged`.  The RT tasks of
    ! this iteration and the next iteration's nlte updates are then skipped
    ! (m3d_patch_t%pre_update), and each patch finishes with this solve.  The
    ! count is taken before a patch finishes early (below), which must not
    ! remove it from the iteration it has just solved.
    !---------------------------------------------------------------------------
    if (self%iter > 0) then
      !$omp critical (nlte_sweep)
      if (self%iter /= sweep_iter) then
        if (sweep_updates /= sweep_active) then
          call shared%abort ('nlte_patch_t%execute: an iteration started before the previous one had completed')
        endif
        sweep_iter = self%iter
        nlte_last_iter = self%iter
        sweep_active = n_nlte_patches - n_finished
        sweep_updates = 0
        last_maxerr = sweep_maxerr
        sweep_maxerr = 0d0
      endif
      if (abs(self%maxerr) > abs(sweep_maxerr)) sweep_maxerr = self%maxerr
      sweep_updates = sweep_updates + 1
      last = sweep_updates == sweep_active
      if (last) nlte_last_change = last_maxerr
      !$omp end critical (nlte_sweep)
      self%global = last_maxerr
      if (last .and. abs(self%global) > 0 .and. abs(self%global) < convlim &
          .and. .not. flux_convergence .and. self%iter < self%maxiter) then
        converged = .true.
      endif
      if (last .and. .not. converged .and. self%iter < self%maxiter) &
        call self%check_stop_request()
    endif

    ! locally_converged is population/convmax-based (premature_finish shortcut);
    ! it must not terminate the task ahead of the flux metric when
    ! flux_convergence is the authoritative one (see the matching guard in
    ! rt_patch_pre_update)
    if (locally_converged .and. .not. flux_convergence) then
      !$omp atomic
      n_convmax_finished = n_convmax_finished + 1
      call self%finish()
    endif

    if (self%iter == 0) then
      call self%init_departures()
      ! a stop request left over from an earlier run in this folder is void
      !$omp critical (stop_flag)
      if (.not. stale_stop_checked) then
        stale_stop_checked = .true.
        if (os%exists(trim(io%run_dir)//'stop_nlte.flag')) then
          open (newunit=lunit, file=trim(io%run_dir)//'stop_nlte.flag', status='old')
          close (lunit, status='delete')
          call console_warn (ch_m3d, 'removed a stop_nlte.flag left from an earlier run')
        endif
      endif
      !$omp end critical (stop_flag)
    endif

    if (save_populations .and. save_every > 0 .and. self%iter > 0 .and. &
        .not. self%is_finished()) then
      if (mod(self%iter, save_every) == 0) call self%write_populations(intermediate=.true.)
    endif

    ! Flux-convergence wind-down (one more population solve):
    ! rt_patch raised flux_conv_signal in the PREVIOUS iteration, after nlte_patch
    ! had already solved these populations against the old radiation field. This
    ! execute() has now re-solved them against the converged field, so it is the
    ! final solve for this z-patch -- finish it. `converged` (which freezes the
    ! remaining rt tasks) is only set once EVERY z-patch has performed this final
    ! solve (count == dims(3)); setting it earlier would let the first patch's
    ! solve make the others skip theirs via pre_update.
    if (flux_convergence) then
      !$omp critical (flux_final)
      if (flux_conv_signal%update(+0) > 0) then
        count = flux_final_counter%update(+1)
      else
        count = 0
      endif
      !$omp end critical (flux_final)
      if (count > 0) then
        if (count >= self%task_list%dims(3) .and. self%iter < self%maxiter) then
          converged = .true.
        endif
        call self%finish()
      endif
    endif

    if (self%is_finished() .or. self%iter==self%maxiter) then
      do i=1, self%natom
        if (.not. allocated(self%atoms(i)%data%wphi)) cycle
        do j=1, size(self%atoms(i)%data%wphi)
          if (allocated(self%atoms(i)%data%wphi(j)%w)) &
            call os%deallocate(self%atoms(i)%data%wphi(j)%w, label='wphi')
        end do
        deallocate(self%atoms(i)%data%wphi)
      end do
    endif
    call trace%end(itimer)
  END SUBROUTINE execute

  !===============================================================================
  !> Bytes of the per-atom arrays this patch holds while the NLTE iteration runs
  !> and while the output pass runs.  The populations, their history and the
  !> Doppler factor stay to the end; the radiative rate accumulators, one
  !> per radiative transition, exist only during the iteration.  Only the
  !> iterated atoms carry departure coefficients, and the profile integrals of
  !> the lines that window edges cut are left out (a few lines per window edge).
  !===============================================================================
  SUBROUTINE estimate_atom_memory (self, nlte, output)
    class(nlte_patch_t), intent(in) :: self
    real(8), intent(out) :: nlte, output
    real(8) :: cell, bytes
    integer :: i, nl
    !-----------------------------------------------------------------------------
    cell = product(real(self%n,8))
    nlte = 0d0
    output = 0d0
    do i=1, self%natom
      nl = self%atoms(i)%nlevel
      bytes = 8d0*nl + 8d0 + 4d0
      if (self%maxiter > 0) then
        bytes = bytes + 8d0*nl
        if (ng_step >= 0) bytes = bytes + 3*8d0*nl
      endif
      output = output + cell*bytes
      if (self%maxiter > 0) bytes = bytes + 2*4d0*self%atoms(i)%nrad
      nlte = nlte + cell*bytes
    end do
  END SUBROUTINE estimate_atom_memory

  !===============================================================================
  SUBROUTINE finish (self)
    class(nlte_patch_t):: self
    class(atom_type), pointer :: atom
    class(atom_data_t), pointer :: adata
    class(composition_type), pointer :: comp
    real(4), dimension(:,:,:,:,:,:), allocatable :: out_mem
    integer :: i, j, l, n
    !-----------------------------------------------------------------------------
    if (self%is_finished()) return

    ! First we move departures from the atom data to the departures structure
    l = 1
    do i=1, self%natom
      atom => self%atoms(i)
      adata => atom%data
      if (allocated(adata%b)) then
        self%departures(atom%idep, 1)%mem = adata%b
      else
        self%departures(atom%idep, 1)%mem = 1.0
      endif
      if (allocated(adata%lvlerr)) then
        do j=1, atom%nlevel
          rank_lvlerr(l) = merge(rank_lvlerr(l), adata%lvlerr(j), abs(rank_lvlerr(l)) > abs(adata%lvlerr(j)))
          l = l +1
        end do
      endif
    end do

    ! Then we fill the project the departures onto the ray aligned grid
    if (rotate_atmos) then
      call self%project_populations()
    endif

    ! finally we write the departures to file
    if (save_populations) then
      call self%write_populations()
    endif

    call self%m3d_patch_t%finish()
    !$omp atomic capture
    n_finished = n_finished + 1
    n = n_finished
    !$omp end atomic
    if (n == self%task_list%dims(3)) then
      call self%save_convergence_report()
    endif
    !---------------------------------------------------------------------------
    ! Once the last nlte_patch of this rank has converged, the populations there
    ! are final, and the output-spectrum tasks that have been held back all along
    ! (see m3d_patch_t%is_ahead_of) can go ahead.  Ranks are released one by one,
    ! which is safe: the long-characteristics sweep of a rank still waits for the
    ! rank below it through its upstream nbor, as it always does.
    !---------------------------------------------------------------------------
    if (n == n_nlte_patches) then
      !$omp atomic write
      nlte_finished = .true.
      call release_output_tasks (self%task_list)
    endif
  END SUBROUTINE finish

  !===============================================================================
  SUBROUTINE init_departures(self)
    class(nlte_patch_t), target :: self
    class(departure_pointer_t), pointer :: dep
    class(atom_type), pointer :: atom
    class(atom_data_t), pointer :: adata
    character(len=256) :: run_dir, save_folder
    character(len=6), allocatable :: keys(:)
    real(4), dimension(:,:,:,:,:,:), allocatable, target :: mem
    integer :: i, j, k, nprecomp, nmu, n(4), x, y, z, l, lunit, nzero
    character(len=160) :: msg
    integer :: itimer=0
    logical :: found
    !-----------------------------------------------------------------------------
    call trace%begin('composition_patch_t%init_departures', itimer=itimer)

    nprecomp = 0
    if (read_departures) then
      run_dir = atmos_file
      save_folder = trim(run_dir)//'/save/'
      if (.not. os%exists_dir(save_folder)) then
        call console_warn (ch_m3d, 'no departure folder found: '//trim(save_folder))
      endif
    endif

    ! Discover which of the requested atoms have precomputed departures by probing each
    ! one's saved metadata file directly. This replaces os%get_filenames, whose shared
    ! shell temp file ('.temp_files<thread>') races across the concurrent dispatch.x
    ! processes the synthesis driver launches — they share a working directory and clobber
    ! each other's temp file. inquire is pure Fortran: no shell, no temp file, no race.
    allocate(keys(self%natom))
    do i=1, self%natom
      atom => self%atoms(i)
      found = .false.
      if (read_departures) then
        inquire(file=trim(save_folder)//'pop_'//trim(atom%id)//'_meta.txt', exist=found)
      endif
      atom%precomputed = found
      if (found) then
        nprecomp = nprecomp + 1
        keys(nprecomp) = atom%id
        atom%idep = nprecomp
      endif
    end do
    j = nprecomp
    do i=1, self%natom
      atom => self%atoms(i)
      if (.not. atom%precomputed) then
        j = j + 1
        keys(j) = atom%id
        atom%idep = j
      endif
    end do

    nmu = merge(size(rt_out%vec, 2), 0, rotate_atmos)
    allocate(self%departures(j, 1+nmu))
    self%m3d_parent%departures => self%departures

    do i=1, nprecomp
      call self%load6('pop_'//trim(keys(i)), run_dir=run_dir, alloc=mem)
      dep => self%departures(i,1)
      allocate(dep%mem, source=mem(:,:,:,:,2,1))
      deallocate(mem)               ! load6 allocates it afresh for the next atom
      call os%mem_register(storage_size(dep%mem), size(dep%mem), label='departures')
      dep%nlevel = size(mem,4)
      allocate(dep%ev(dep%nlevel))
      open(newunit=lunit, file=trim(save_folder)//'pop_'//trim(keys(i))//'_ev.txt', form='formatted')
      read(lunit, '(A)') ! skip header
      read(lunit, *) dep%ev
      close(lunit)
    end do

    do i=1, self%natom
      atom => self%atoms(i)
      adata => atom%data
      ! every departure set carries the species key of its atom -- precomputed or not,
      ! and on every rotated copy; an unkeyed set would match no particle at all
      do k=1, 1+nmu
        call self%departures(atom%idep, k)%set_key(atom%id)
      end do
      dep => self%departures(atom%idep, 1)
      n = shape(adata%nstar)
      dep%nlevel = atom%nlevel
      if (atom%precomputed) then
        if (.not. allocated(adata%b)) then
          call os%allocate(adata%b, shape(adata%nstar), label='b')
        endif
        !-----------------------------------------------------------------------
        ! The departures were written as real(4), so any level that sat below
        ! tiny(1.0) came back as an exact zero.  Zero is a trap here: the update
        ! is multiplicative (b = b*e), so such a level can never be repopulated,
        ! and its column of the rate matrix is identically zero.  Restore it to
        ! the smallest value the file could have represented; the population
        ! floor in update_populations then lifts it to something sensible.
        !-----------------------------------------------------------------------
        nzero = 0
        do l=1, n(4)
          do z=1, n(3)
            do y=1, n(2)
              do x=1, n(1)
                if (dep%mem(x,y,z,l) > 0.0) then
                  adata%b(x,y,z,l) = dep%mem(x,y,z,l)
                else
                  adata%b(x,y,z,l) = tiny(1.0)
                  nzero = nzero + 1
                endif
              end do
            end do
          end do
        end do
        if (nzero > 0) then
          write(msg,'(I0,A,A,A)') nzero, &
            ' underflowed departure coefficients read for ', trim(atom%id), &
            ' -- the run they came from had not converged'
          call console_warn (ch_m3d, msg)
        endif
      else
        allocate(dep%ev, source=atom%ev_keys())
        allocate(dep%mem(n(1), n(2), n(3), n(4)))
        call os%mem_register(storage_size(dep%mem), size(dep%mem), label='departures')
      endif
    end do
    call trace%end(itimer)
  END SUBROUTINE init_departures

  !===============================================================================
  !> Write this patch's LTE populations and departure coefficients to
  !> save/pop_<atom>_*.  An intermediate write (save_every) stores the unrotated
  !> populations only; the ray-aligned slots keep 1 until the final write, which
  !> follows project_populations.  Each write replaces the patch's file, and the
  !> meta row is added once.
  SUBROUTINE write_populations (self, intermediate)
    class(nlte_patch_t), intent(inout), target :: self
    logical, intent(in), optional :: intermediate
    class(atom_type), pointer :: atom
    class(atom_data_t), pointer :: adata
    class(departure_pointer_t), pointer :: dep
    real(4), dimension(:,:,:,:,:,:), allocatable :: out_mem
    integer :: x, y, z, l, m, n(4), i, j, k, lunit, nmu
    real(8), allocatable :: ev_key(:)
    !---------------------------------------------------------------------------------
    call trace%begin('nlte_patch_t%write_populations')
    do i=1, self%natom
      atom => self%atoms(i)
      adata => atom%data
      n = shape(adata%nstar)
      nmu = size(self%departures, 2)
      allocate(out_mem(n(1), n(2), n(3), n(4), 2, nmu), source=1.0)
      do m=1, merge(1, nmu, shared%option(intermediate, .false.))
        if (m>1) adata => atom%rotated_data(m-1)
        do l=1, n(4)
          do z=1, n(3)
            do y=1, n(2)
              do x=1, n(1)
                out_mem(x,y,z,l,1,m) = real(adata%nstar(x,y,z,l))
                if (allocated(adata%b)) then
                  out_mem(x,y,z,l,2,m) = real(adata%b(x,y,z,l))
                end if
              end do
            end do
          end do
        end do
      end do

      call self%save6(out_mem, label='pop_'//trim(atom%id), add2meta=.not. self%pop_meta_written)
      !$omp critical (write_atom_ev)
      if (.not. os%exists(trim(io%run_dir)//'save/pop_'//trim(atom%id)//'_ev.txt')) then
        open(newunit=lunit, file=trim(io%run_dir)//'save/pop_'//trim(atom%id)//'_ev.txt', form='formatted')
        write(lunit, '(A)') '# Level Energies (eV) + 100*Ionization Stage'
        allocate(ev_key, source=atom%ev_keys())   ! explicit: the build assumes norealloc_lhs
        do l=1, atom%nlevel
          write(lunit, '(F16.6)') ev_key(l)
        end do
        deallocate(ev_key)
        close(lunit)
      endif
      !$omp end critical (write_atom_ev)
      deallocate(out_mem)
    end do
    self%pop_meta_written = .true.
    call trace%end()
  END SUBROUTINE write_populations

  !===============================================================================
  !> A stop request: the file <run_dir>/stop_nlte.flag, looked for by the last
  !> nlte_patch of each iteration.  It ends the iterations the way convergence
  !> does, so every depth patch finishes with this iteration's solve and writes
  !> its departures.  The output pass then runs, unless the file's first word is
  !> 'exit'.  The file is removed once read.
  !===============================================================================
  SUBROUTINE check_stop_request (self)
    class(nlte_patch_t), intent(in) :: self
    character(len=256) :: file
    character(len=16) :: word, nit
    logical :: exists
    integer :: lunit, iostat
    !-----------------------------------------------------------------------------
    file = trim(io%run_dir)//'stop_nlte.flag'
    inquire (file=trim(file), exist=exists)
    if (.not. exists) return
    word = ''
    open (newunit=lunit, file=trim(file), form='formatted', status='old', iostat=iostat)
    if (iostat == 0) then
      read (lunit, *, iostat=iostat) word
      close (lunit, status='delete')
    endif
    stop_exit = trim(word) == 'exit'
    stop_requested = .true.
    converged = .true.
    write (nit, '(i0)') self%iter
    if (stop_exit) then
      call console_warn (ch_m3d, 'stop requested after iteration '//trim(nit)// &
        ': departures are written and the output pass is skipped')
    else
      call console_warn (ch_m3d, 'stop requested after iteration '//trim(nit)// &
        ': departures are written, then the output pass runs')
    endif
  END SUBROUTINE check_stop_request

  !====================================================================================
  SUBROUTINE project_populations(self)
    class(nlte_patch_t), intent(inout), target :: self
    class(atom_type), pointer :: atom
    class(composition_type), pointer :: comp
    class(atom_data_t), pointer :: adata
    class(departure_pointer_t), pointer :: dep
    real(8), dimension(:,:,:), allocatable :: norm, totn
    real(8) :: vec(3), tmp
    integer :: i,j,k,l, n(3), imu, nmu, x, y, z
    logical :: iterated, has_departures
    !---------------------------------------------------------------------------------
    nmu = size(rt_out%vec, 2)
    n = self%n
    allocate(norm(n(1),n(2),n(3)))
    allocate(totn(n(1),n(2),n(3)))

    do i=1, size(self%departures, 1)
      do imu=1, nmu
        vec = rt_out%vec(:,imu)
        dep => self%departures(i,1+imu)
        allocate(dep%mem, source=self%departures(i,1)%mem)
        call os%mem_register(storage_size(dep%mem), size(dep%mem), label='departures')
        ! the level energies too: line-list lines of the species find their levels
        ! by energy on every ray-aligned copy (linelist_mod, find_closest)
        if (.not. allocated(dep%ev)) allocate(dep%ev, source=self%departures(i,1)%ev)
        dep%nlevel = self%departures(i,1)%nlevel
        do l=1, size(dep%mem, 4)
          call self%rotate_cube4(vec, dep%mem(:,:,:,l), do_log=.true.)
        end do
      end do
    end do

    do i=1, self%natom
      atom => self%atoms(i)
      adata => atom%data
      iterated = allocated(adata%b)

      do imu=1, nmu
        vec = rt_out%vec(:,imu)
        adata => atom%rotated_data(imu)
        comp => self%m3d_parent%rotated_composition(imu)
        dep => self%departures(atom%idep, 1+imu)
        call atom%init_atom_data(comp, adata, iterate=.false.)

        if (iterated) then
          totn = 0
          call os%allocate(adata%b, shape(adata%nstar), label='b')

          do l=1, atom%nlevel
            do z=1, n(3)
              do y=1, n(2)
                do x=1, n(1)
                  adata%b(x,y,z,l) = dep%mem(x,y,z,l)
                  totn(x,y,z) = totn(x,y,z) + adata%b(x,y,z,l) * adata%nstar(x,y,z,l)
                end do
              end do
            end do
          end do

          do z=1, n(3)
            do y=1, n(2)
              do x=1, n(1)
                norm(x,y,z) = adata%totn(x,y,z) / totn(x,y,z)
              end do
            end do
          end do

          do l=1, atom%nlevel
            do z=1, n(3)
              do y=1, n(2)
                do x=1, n(1)
                  tmp = adata%b(x,y,z,l) * norm(x,y,z)
                  adata%b(x,y,z,l) = merge(tmp, 0.d0, tmp==tmp) ! set to zero if NaN
                end do
              end do
            end do
          end do
        endif
      end do
    end do
  END SUBROUTINE project_populations

  !===============================================================================
  SUBROUTINE save_convergence_report (self)
    class(nlte_patch_t) :: self
    integer :: lunit, i
    character(len=256) :: conv_file
    integer, save :: itimer=0
    logical, save :: first_time = .true.
    !-----------------------------------------------------------------------------
    if (first_time .and. allocated(rank_lvlerr)) then
      first_time = .false.
      call trace%begin('nlte_patch_t%save_convergence_report', itimer=itimer)
      write(conv_file, '(A,A,I3.3,A)') trim(io%run_dir), 'convergence_report_', self%rank,'.txt'
      open(newunit=lunit, file=conv_file, form='formatted')
      write(lunit, '(A6, A12)') 'lvl', 'maxconv'
      do i=1, size(rank_lvlerr)
        write(lunit, '(I6, ES12.4)') i, rank_lvlerr(i)
      end do
      close(lunit)
      call trace%end(itimer)
    endif
  END SUBROUTINE save_convergence_report

  !===============================================================================
  SUBROUTINE print_post_update_info (self)
    class(nlte_patch_t):: self
    class(atom_data_t), pointer :: adata
    character(len=64) :: text
    character(len=1) :: ng_mark
    integer :: i, flux_pos
    real(8) :: global
    logical, save :: first_time = .true.
    !-----------------------------------------------------------------------------
    if (console_on(ch_m3d, 2, all_ranks=.true.)) then
      if (self%iter == 0) then
        write(text,'(A,I2,A,I2,A)') 'LTE Iteration .'
      elseif (converged) then
        write(text,'(A)') 'Converged!'
      else
        write(text,'(A,I2,A,I2,A)') 'Iteration ', self%iter, '/', self%maxiter, '.'
      endif
      call io%header(text)
      if (self%update_skipped) then
        call self%print_to_stdout(text='update skipped')
      elseif (self%m3d_patch_t%check_in_bounds()) then
        call self%print_to_stdout(arr=memcheck(1:4))
      else
        call self%print_to_stdout()
      endif
    endif
    if (associated(self%atoms) .and. console_on(ch_m3d, 0, all_ranks=.true.)) then
      if (.not. console_on(ch_m3d, 2, all_ranks=.true.)) then
        if (first_time) then
          first_time = .false.
          write(console_unit, '(A8, A10, A5, A14)') 'iter', 'global', 'lvl', 'pos'
        endif
        ! an asterisk after the iteration number marks an Ng-accelerated iteration
        ng_mark = ' '
        do i=1, self%natom
          if (self%atoms(i)%ng_applied(self%iter)) ng_mark = '*'
        end do
        write(console_unit, '(I3, A1, "/", I3, F8.2, " %", I5, " [", 3I4, "]")') &
          self%iter, ng_mark, self%maxiter, abs(self%global*100), elvl, epos
      else
        write(console_unit,'(A, ES12.4, A10, ES12.4, A, I4, A, "[", 3I4, "] ", A)') &
          '  local =', self%local,  ', global =', self%global, &
          ', lvl =', elvl, ', pos=', epos, self%atoms(eatom)%id
        flux_pos = global_flux_pos%update(+0)
        if (flux_pos >= 0) then
          write(console_unit,'(A, ES12.4, A, I4, A, F10.2, A)') &
            '  flux global =', global_flux_err%update(0.d0), ', i_nu =', flux_pos, &
            ', lambda =', global_flux_lambda%update(0.d0), ' AA'
        endif
        write(console_unit,'(1X, A80)') &
          '................................................................................'
      endif
    endif
  END SUBROUTINE print_post_update_info

  !===============================================================================
  SUBROUTINE get_osmet_nu (self, sp)
    class(nlte_patch_t), target :: self
    class(atom_type), pointer :: atom
    class(atomic_spectrum_type) :: sp
    integer :: i, j, fs, fe
    !-----------------------------------------------------------------------------
    if (.not. osmet%initialised) return
    do j=1, self%natom
      atom => self%atoms(j)
      do i=1, atom%ncont
        associate(cont => atom%trans%cont(i))
          fs = searchsorted(osmet%lam, cont%lam_min, outbounds_ok=.true.)+1
          fe = searchsorted(osmet%lam, cont%lam_max, outbounds_ok=.true.)
          osmet%active(fs:fe) = .true.
        end associate
      end do
    end do
    call sp%reset(lam=pack(osmet%lam, osmet%active))
  END SUBROUTINE get_osmet_nu

  !===============================================================================
  SUBROUTINE exclude_double_continua (self, atom)
    class(nlte_patch_t), target :: self
    type(atom_type) :: atom
    class(bftrans_type), pointer :: trans
    integer :: i, lower, kc
    logical :: ok
    logical, save :: first_time=.true.
    !-----------------------------------------------------------------------------
    do i=1, size(composition%trans)
      trans => composition%trans(i)
      if (.not. trans%active) cycle
      select case(trans%label)
       case ('H+H', 'H2+', 'H- FF', 'H FF', 'H2-', 'H- BF')
        ! do not exclude
       case default
        ok = .false.
        if (trans%level%g > 0) then
          if (uppercase(trans%ion%id(1:2)) == atom%id) then
            if (trans%molecular) then
              ok = .true.
            else
              do kc=1, atom%ncont
                if (.not. atom%trans%cont(kc)%active) cycle
                lower = atom%trans%cont(kc)%i

                if (atom%ev(lower)==0) then
                  ok = trans%level%xion == atom%ev(lower)
                else
                  ok = abs(trans%level%xion / atom%ev(lower) - 1) < 1d-2
                endif
                if (ok) exit
              end do
            endif
            if (ok) then
              trans%active = .false.
              if (first_time) then
                if (console_on(ch_m3d, 2)) write(console_unit,'(2a)') 'excluding transition ', trim(trans%label)
              endif
              exit
            endif
          endif
        endif
      end select
    end do
    if (first_time) first_time=.false.
  END SUBROUTINE exclude_double_continua

  !===============================================================================
  ! Prevent a model atom's bound-bound transitions from being counted a second time
  ! through the bulk line list.  exclude_from_line_list (the default) drops the atom's
  ! species from the line list wholesale via in_linelist, which both the metal and
  ! the hydrogen kernel honour.  Otherwise only the records that the atom
  ! reproduces are deactivated, so line-list lines the atom does not carry survive; the
  ! match is on rest frequency and gf, taken from the atom's own SoA copy (built on the
  ! output transitions in pre_init_atom) so both sides carry the same conventions --
  ! isotope fraction divided out and fine-structure components split.
  SUBROUTINE exclude_double_lines (self)
    class(nlte_patch_t), target :: self
    class(massive_type), pointer :: species
    class(element_type), pointer :: element
    class(molec_type), pointer :: molec
    integer :: j, ia
    !-----------------------------------------------------------------------------
    if (line_list_precedence) then
      if (console_on(ch_m3d, 1)) write(console_unit,'(a)') &
        'output pass: lines of the atoms'' species come from the line list'
      return
    endif
    do ia=1, self%natom
      associate(atom => self%atoms(ia))
        select type (species => atom%species)
         class is (element_type)
          element => species
          if (exclude_from_line_list) then
            element%ions%in_linelist = .false.
          else
            do j=lbound(element%ions,1), ubound(element%ions,1)
              call deactivate_duplicates(element%ions(j), self%output_atoms(ia)%trans%ll)
            end do
          endif
         class is (molec_type)
          molec => species
          if (exclude_from_line_list) molec%in_linelist = .false.
        end select
        if (exclude_from_line_list) then
          if (console_on(ch_m3d, 1)) write(console_unit,'(3a)') 'excluding all ', trim(atom%id), ' lines from line list'
        endif
      end associate
    end do
  contains
    !---------------------------------------------------------------------------
    ! Deactivate every line of ptc's bulk line list that atom_ll also carries.
    SUBROUTINE deactivate_duplicates (ptc, atom_ll)
      class(particle_type), intent(inout) :: ptc
      type(linelist_type), intent(in) :: atom_ll
      integer :: iso, l, k
      !-------------------------------------------------------------------------
      if (.not. allocated(ptc%ll)) return
      if (atom_ll%nline == 0) return
      do iso=lbound(ptc%ll,1), ubound(ptc%ll,1)
        associate (ll => ptc%ll(iso))
          do l=1, ll%nline
            if (.not. ll%active(l)) cycle
            do k=1, atom_ll%nline
              if (abs(atom_ll%nu0(k) / ll%nu0(l) - 1) < 1e-6 .and. &
                  abs(atom_ll%loggf(k) - ll%loggf(l)) < 1e-4) then
                ll%active(l) = .false.
                if (console_on(ch_m3d, 2)) write(console_unit,'(A, F10.2)') 'excluding linelist line at ', CC / ll%nu0(l) * 1e8
                exit
              endif
            end do
          end do
        end associate
      end do
    END SUBROUTINE deactivate_duplicates
  END SUBROUTINE exclude_double_lines

  !===============================================================================
  SUBROUTINE nlte_patch_nbor_relations(self, another, needed, needs_me, download)
    class(nlte_patch_t) :: self
    class(task_t), pointer :: another         ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed,
    !    does it need me, should I download?
    ! ----------------------------------------------------------------------------
    call trace%begin ('nlte_patch_t%nbor_relations')
    if (another%kind == 'rt_column') then
      ! every column sweeps this slab: an update needs the rates of all of
      ! them, and each of their iterations needs the populations
      needed = .true.; needs_me = .true.; download = .false.
      call trace%end()
      return
    endif
    call self%m3d_patch_t%nbor_relations(another, needed, needs_me, download)

    if (needed .or. needs_me) then
      if (another%id == self%m3d_parent%id) then
        ! the composition patch needs this to be able to
        ! update again before final spectrum is printed
      else
        needed = .false. ; needs_me = .false.
      endif
    endif

    call trace%end()
  END SUBROUTINE nlte_patch_nbor_relations

  !===============================================================================
  FUNCTION check_in_bounds (self) result (cib)
    class(nlte_patch_t) :: self
    logical :: cib
    !-----------------------------------------------------------------------------
    cib = self%m3d_patch_t%check_in_bounds()
    if (associated(self%atoms)) then
      cib = cib .and. lvlcheck > 0 .and. lvlcheck <= self%atoms(atomcheck)%nlevel
    endif
  END FUNCTION check_in_bounds

  !===============================================================================
  SUBROUTINE init_nbors_hash (self, nbors)
    class(nlte_patch_t):: self
    class(nbor_t), pointer:: nbors
    integer, save:: itimer=0
    integer :: key(4), i_nu, i_mu
    !-----------------------------------------------------------------------------
    call self%m3d_patch_t%init_nbors_hash(nbors)
    key(1) = 0
    key(4) = hash_string('rt_column')
    do i_mu=-1, 2
      do i_nu=-n_windows_out, 2*n_windows
        key(2) = i_nu
        key(3) = i_mu
        call self%check_nbor(key, nbors)
      end do
    end do
  END SUBROUTINE init_nbors_hash

  !===============================================================================
  SUBROUTINE save_abundance_to_file (self, atom)
    class(nlte_patch_t), target :: self
    class(atom_type) :: atom
    integer :: lunit, i
    !-----------------------------------------------------------------------------
    if (self%ipos(3)==0) then
      open(newunit=lunit, file=trim(io%run_dir)//'abund.txt', form='FORMATTED', access='append')
      write(lunit, '(A3, F8.4)') atom%id, atom%abnd
      close(lunit)
    endif
  END SUBROUTINE save_abundance_to_file

  !===============================================================================
  SUBROUTINE check_mem_NaN (self)
    class(nlte_patch_t):: self
    integer :: i,j,k,v,t,w
    integer, save :: itimer
    !-----------------------------------------------------------------------------
    if (.not. check_NaN) return
    call trace%begin ('nlte_patch_t%check_mem_NaN', itimer=itimer)
    do w=1, self%natom
      associate(atom=>self%atoms(w))
        if (allocated(atom%data%b)) then
          call check_array(atom%data%b, trim(atom%id)//' departure coefficients', check_zero=.false., abort=.true.)
        endif
      end associate
    end do
    call trace%end(itimer)
  END SUBROUTINE check_mem_NaN

END MODULE
