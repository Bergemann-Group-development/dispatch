MODULE atom_patch_mod
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
                        fmin, cbbesc, cbfesc, cbbhsc, cbfhsc, ng_step, ng_start
    USE trans_mod, only: line_type, cont_type, transition_type
    USE const_mod, only: CC
    USE ion_mod, only: ion_type, particle_type, debye_shield, massive_type
    USE element_mod, only: element_type
    USE molecule_mod, only: molec_type
    USE bbtrans_mod, only: bbtrans_type
    USE composition_mod, only: composition_type
    USE m3d_patch_mod, only: m3d_patch_t, composition, maxiter, key_link_t, &
    global_maxerr, convlim, convmax, atom_update_counter, comoving_frame, global_pos, &
    decouple_continuum, overlapping_bb, overlapping_bf, converged, premature_finish, overrelax, &
    rotate_atmos, long_vec, rank_lvlerr, overlapping_bb_in_output, long_scheme, n_nu, detailed_absmet
    USE bftrans_mod, only: bftrans_type
    USE osmet_mod, only: osmet_type, osmet
    USE nlte_patch_mod, only: nlte_patch_t, n_cntrbf, cntrbf_lines
    USE m3d_debug_mod
    USE shared_mod
    USE, intrinsic :: ieee_arithmetic
    implicit none
  
    !-----------------------------------------------------------------------------
  
    type, public, extends(m3d_patch_t) :: atom_patch_t
      type(atom_type), dimension(:), pointer :: output_atom, bb_atom, bf_atom
      class(atomic_spectrum_type), pointer :: output_sp=>null()
      class(atomic_spectrum_type), pointer :: bf_sp=>null()
      class(atomic_spectrum_type), pointer :: bb_sp=>null()
      real(8) :: maxerr=huge(1.0), prev_global=huge(1.d0), prev_prev_global=huge(1.d0)
      real(8), dimension(:), allocatable :: lvlerr
      real(8) :: local, sor_w
      integer :: natom=0
      logical :: update_skipped=.false., atom_allocated=.false., molec
      logical :: overrelax=.false.
    contains
      procedure :: pre_init => atom_patch_pre_init
      procedure :: execute
      procedure :: finish
      procedure :: project_populations
      procedure :: nbor_relations => atom_patch_nbor_relations
      procedure :: pre_init_atom
      procedure :: pre_init_spectrum
      procedure :: get_osmet_nu
      procedure :: exclude_double_continua
      procedure :: exclude_double_lines
      procedure :: print_post_update_info
      procedure :: check_in_bounds
      procedure :: init_nbors_hash
      procedure :: save_abundance_to_file
      procedure :: save_convergence_report
      !procedure :: load_mem
      procedure :: allocate_mem
      procedure :: deallocate_mem
      procedure :: check_mem_NaN
    end type atom_patch_t

    integer, save :: epos(3), elvl, eatom=1
    logical, save :: exclude_from_line_list=.true.
    logical, save :: save_populations=.true.
    integer, save :: save_rates=0
    character(len=256), save :: precomputed_depart=''
    logical, public, save :: external_sp=.false.
    !-----------------------------------------------------------------------------
  
  contains
    
  !===============================================================================
  SUBROUTINE atom_patch_pre_init (self)
    class(atom_patch_t), target:: self
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%atom_patch_pre_init', itimer=itimer)
    self%kind = 'atom_patch'
    call self%m3d_patch_t%pre_init()
    if (.not. self%atom_allocated) then
      self%maxiter = -1 ! prevents addition to task list
      call trace%end(itimer)
      return
    endif
    call self%set_nv_and_nw(sum(self%atom%nlevel), 2)
    self%static_mem=.true.

    call self%post_pre_init()
    call trace%end(itimer)
  END SUBROUTINE atom_patch_pre_init

  !===============================================================================
  SUBROUTINE pre_init_atom (self)
    class(atom_patch_t) :: self
    type(atom_type), save, dimension(:), allocatable :: atom
    type(atom_trans_t), save, dimension(:), pointer :: dc_trans, out_trans
    class(element_type), pointer :: element
    class(molec_type), pointer :: molecule
    class(ion_type), pointer :: ion
    class(massive_type), pointer :: species
    character(len=256), save, dimension(:), allocatable :: atom_file
    real(8), dimension(:), allocatable :: abundance
    real(8), save :: qnorm=10.0
    integer, save :: itimer=0, natom=0, n_in_lines=0
    integer, save, allocatable :: in_lines(:)
    logical, save :: first_time=.true., exclude_trace_cont=.true., use_atom_abnd=.true.
    logical, save :: relative_abundance = .false.
    logical :: pandas_format
    integer :: i, j, iostat
    namelist /atom_params/ atom_file, qnorm, convlim, exclude_trace_cont, level_partf, &
    use_atom_abnd, logf_cut, red_cut, blue_cut, n_in_lines, abundance, relative_abundance, &
    exclude_from_line_list, precomputed_depart, use_classic_col, kramer_bf, overlapping_bf, &
    overlapping_bb, detailed_absmet, n_cntrbf, convmax, premature_finish, &
    save_populations, nq_fac, nq_min, nq_max, cbbesc, cbfesc, cbbhsc, cbfhsc, &
    disable_bf, disable_bb, overlapping_bb_in_output, save_rates, fmin
    namelist /line_mask/ in_lines, cntrbf_lines
    !-----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%pre_init_atom', itimer=itimer)
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

      natom = count(atom_file /= '')
      if (natom == 0) then
        print *, 'WARNING: no atom specified'
        call trace%end(itimer)
        return
      else
         ! this is the atom to be read
        allocate(atom(natom))

        allocate(out_trans(natom))
        if (decouple_continuum) then
          allocate(self%bb_atom(natom))
          allocate(self%bf_atom(natom))
          allocate(dc_trans(natom))
        endif

        atom(atomcheck)%check = .true.
        call truncate_array(atom_file, natom)
        call truncate_array(abundance, natom)
      endif

      if (n_in_lines+n_cntrbf > 0) then
        allocate(in_lines(n_in_lines), cntrbf_lines(n_cntrbf))
        rewind (io%input)
        read (io%input, line_mask, iostat=iostat)
        if (iostat > 0) then
          rewind (io%input)
          read (io%input, line_mask)
        end if
        write (io%namelist, line_mask)
      endif

      write (*,*) ' - Reading atoms'
      do i=1, natom
        associate(atom_file=>atom_file(i), atom=>atom(i))
          j = rindex(atom_file, '/')+1
          pandas_format = .not. startswith(atom_file(j:len_trim(atom_file)), 'atom.')
          if (pandas_format) then
            level_partf = .false.
            call atom%read_pandas(atom_file, qnorm, in_lines=in_lines)
            call atom%set_colliders(composition)
          else
            call atom%read(atom_file, qnorm, in_lines=in_lines, read_coll=.true.)
          endif
          atom%name = 'input'
          atom%verbose = self%verbose

          call composition%element_by_id(atom%id(1:3), element, verbose=.false.)
          if (associated(element)) then
            species => element
            element%ions%in_opac = .true.
            atom%id = element%id
          else
            call composition%molecule_by_id(atom%id(1:3), molecule, verbose=.false.)
            if (associated(molecule)) then
              species => molecule
              molecule%in_opac = .true.
              atom%id = molecule%id
              do j=1, molecule%nions
                ion => molecule%ions(j)%p
                ion%in_opac = .true.
              end do
            endif
          endif
      
          call atom%prepare_damping(species)

          if (abundance(i) > -99) atom%abnd = abundance(i)
          if (use_atom_abnd .or. abundance(i) > -99) then
            call composition%overwrite_abundance(species, atom%abnd, atom%awgt, relative_abundance)
          else
            atom%abnd = atom%abnd + composition%FeH
          endif
          
          if (exclude_trace_cont) then
            call self%exclude_double_continua(atom)
          endif
          call self%exclude_double_lines(atom, species)

          ! here we create copies of the line and cont arrays,
          ! which will have altered ired and iblue for example
          ! they will be pointed at later
          call out_trans(i)%allocate(bb=atom%trans%line)
          if (decouple_continuum) then
            call dc_trans(i)%allocate(bb=atom%trans%line, bf=atom%trans%cont)
          endif

          if (pandas_format) then
            call os%copy_file(trim(atom_file), trim(io%run_dir)//'atoms/'//trim(atom%id)//'.txt')
          else
            call os%copy_file(trim(atom_file), trim(io%run_dir)//'atoms/'//'atom.'//trim(atom%id))
          endif
        end associate
      end do
      allocate(rank_lvlerr(sum(atom%nlevel, 1)), source=0.d0)
      write (io%namelist, atom_params)
    endif

    if (natom == 0) then
      call trace%end(itimer)
      return
    endif

    ! this is the patch atom, which contains local data
    allocate(self%atom(natom))

    ! this is the output atom, which contains local data
    ! shared with self%atom, but with different transitions
    allocate(self%output_atom(natom))

    ! here we create a copy of atom, but the transitions will
    ! still be linked to the original atom
    self%atom = atom
    self%atom%name = 'main'
    self%n_nu = n_nu

    if (decouple_continuum) then
      allocate(self%bb_atom(natom))
      allocate(self%bf_atom(natom))
    endif

    do i=1, natom
      ! this pointer will be referenced by other atoms with same local data
      allocate(self%atom(i)%data)
      if (rotate_atmos) then
        allocate(self%atom(i)%rotated_data(size(long_vec, 2)))
      endif

      ! here we create a copy of the patch atom, which will share the same
      ! atom%data array, but we link it to the output transitions
      self%output_atom(i) = self%atom(i)
      self%output_atom(i)%trans => out_trans(i)
      self%output_atom(i)%name = 'output'

      if (decouple_continuum) then
        self%bb_atom(i) = self%atom(i)
        self%bf_atom(i) = self%atom(i)
        allocate(self%bb_atom(i)%trans) ! this will detach the original pointer and 
        allocate(self%bf_atom(i)%trans) ! create a new one, to be pointed elsewhere
        self%bb_atom(i)%trans%line => dc_trans(i)%line
        self%bf_atom(i)%trans%cont => dc_trans(i)%cont
        self%bb_atom(i)%name = 'bound-bound'
        self%bf_atom(i)%name = 'bound-free'
      endif
    end do
    self%atom_allocated = .true.
    self%natom = natom
    call trace%end(itimer)
  END SUBROUTINE pre_init_atom

  !===============================================================================
  SUBROUTINE pre_init_spectrum (self, output_sp)
    class(atom_patch_t) :: self
    type(atomic_spectrum_type), target :: output_sp
    type(atomic_spectrum_type), save, target :: sp, bb_sp, bf_sp
    logical, save :: first_time=.true.
    logical :: save_sp, saved_already
    integer :: i
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%pre_init_spectrum', itimer=itimer)
    if (first_time) then
      first_time = .false.
      external_sp = output_sp%nnu > 0 
      save_sp = external_sp .or. (associated(self%atom) .and. size(long_vec,2)>0)
      if (associated(self%atom)) then
        if (decouple_continuum) then
          if (detailed_absmet .and. overlapping_bf) call self%get_osmet_nu(bf_sp)
          call setup_spectrum(bf_sp, 'bf', self%bf_atom, bf=.true., bb=.false., sort=overlapping_bf, set_active=.true.)
          call setup_spectrum(bb_sp, 'bb', self%bb_atom, bf=.false., bb=.true., sort=overlapping_bb, set_active=.true.)
        else
          if (detailed_absmet) call self%get_osmet_nu(sp)
          call setup_spectrum(sp, 'all', self%atom, bb=.true., bf=.true., sort=.true., set_active=.true.)
        endif

        if (.not. external_sp) then
          call setup_spectrum(output_sp, 'out', self%output_atom, bb=.true., bf=.false., sort=overlapping_bb_in_output, set_active=.false.)
          saved_already = .true.
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
    subroutine setup_spectrum(sp, name, atoms, bb, bf, sort, set_active)
      class(atomic_spectrum_type), target :: sp
      character(len=*), intent(in) :: name
      class(atom_type), intent(in), dimension(:) :: atoms
      integer, dimension(:,:), pointer, contiguous :: nac, nal
      logical, intent(in) :: bb, bf, sort, set_active
      integer :: i
      !-----------------------------------------------------------------------------
      sp%name = name
      ! allocate list of active continua and lines
      do i=1, size(atoms)
        call atoms(i)%trans%append_transitions(sp, bf=bf, bb=bb)
      end do
      if (sort) call sp%sort()
      call os%allocate(nac, [sp%nnu, size(atoms)], zero=.true.)
      call os%allocate(nal, [sp%nnu, size(atoms)], zero=.true.)
      call os%allocate(sp%na, sp%nnu, zero=.true.)
      do i=1, size(atoms)
        sp%nac => nac(:,i)
        sp%nal => nal(:,i)
        call atoms(i)%trans%set_ired_iblue(sp, bf=bf, bb=bb, set_active=set_active)
      end do
      sp%na = sum(min(nal+nac, 1), dim=2)
      if (save_sp) then
        call sp%save(name)
      endif
      if (set_active) then
        write(*,'(a,I0,a)') 'sp_from_atom: frequency grid has ', sp%nnu, ' points.'
        write(*,'(a,I0,a)') 'sp_from_atom: ', sum(nac+nal), ' active transition-frequency points'
      endif
      call os%deallocate(nac)
      call os%deallocate(nal)
    end subroutine setup_spectrum
  END SUBROUTINE pre_init_spectrum

  !===============================================================================
  SUBROUTINE execute (self)
    class(atom_patch_t) :: self
    class(atom_type), pointer :: atom
    class(atom_data_t), pointer :: adata
    class(composition_type), pointer :: comp
    class(line_type), pointer :: line
    class(cont_type), pointer :: cont
    class(element_type), pointer :: element
    class(molec_type), pointer :: molec
    real(4), dimension(:,:,:), allocatable :: cij, cji
    integer :: i, j, k, l, p(3), m(3), ipos
    integer, save:: count=0
    real(8) :: global
    logical :: do_ng, locally_converged, found
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin (trim(self%kind)//'%update', itimer=itimer)
    do_ng = ng_step > 0 .and. self%iter >= ng_start .and. mod(self%iter-ng_start, ng_step) == 0
    comp => self%m3d_parent%composition
    self%maxerr = 0
    locally_converged = premature_finish
    m = self%n

    if (save_rates>0) then
      allocate(cij(m(1),m(2),m(3)), source=0.0)
      allocate(cji(m(1),m(2),m(3)), source=0.0)
    endif

    do k=1, self%natom
      atom => self%atom(k)
      adata => atom%data

      if (self%iter == 0) then
        call atom%init_atom_data(comp, adata, self%maxiter>0, found)
        if (.not. found) call self%save_abundance_to_file(atom)
        !if (precomputed_depart /= '') then
        !  call self%load_mem(precomputed_depart, self%kind)
        !endif
      elseif (.not. adata%converged) then
        call atom%normalize_profiles(self%iter)
        call atom%update_populations(comp, self%local, self%pcheck, cij, cji)

        if (save_rates==k) then
          call self%save3(real(adata%gamma_ij(:,:,:,linecheck)), label='rij')
          call self%save3(real(adata%gamma_ji(:,:,:,linecheck)), label='rji')
          call self%save3(cij, label='cij')
          call self%save3(cji, label='cji')
        endif
        
        call os%deallocate(adata%gamma_ij, label='atom%gamma_ij')
        call os%deallocate(adata%gamma_ji, label='atom%gamma_ji')

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
          call atom%ng_acceleration()
        endif
      endif
      if (.not. (maxiter==0 .or. long_scheme=='none')) then
        call atom%compute_cont_factors(atom%data)
      endif
      locally_converged = locally_converged .and. adata%converged
    end do

    if (locally_converged) then
      call self%finish()
    endif

    !$omp critical (maxerr_reset)
    self%prev_prev_global = self%prev_global
    self%prev_global = global_maxerr%update(0.d0)
    ipos = global_pos%update(+0)
    if (ipos == self%ipos(3) .or. abs(self%maxerr) > abs(self%prev_global)) then
      if (ipos /= self%ipos(3)) call global_pos%reset(self%ipos(3))
      call global_maxerr%reset(self%maxerr)
    endif
    !$omp end critical (maxerr_reset)
    
    if (self%iter > 0) then
      count = atom_update_counter%update(+1)
    endif
    if (count >= self%task_list%dims(3)) then
      global = global_maxerr%update(0.d0)
      if ((abs(global) > 0) .and. abs(global) < convlim .and. abs(self%maxerr) < convlim) then
        !write(stdout, *) self%ipos(3)+1, 'Converged!', global
        if (self%iter < self%maxiter) converged = .true.
        call self%finish()
      endif
    endif

    if (self%print_update_info) then
      p = self%pcheck
      l = self%linecheck
      adata => self%atom(atomcheck)%data
      memcheck(1) = adata%n(p(1),p(2),p(3),lvlcheck)
      if (allocated(adata%gamma_ij)) then
        if (self%linecheck<=atom%nrad) then
          memcheck(2) = adata%gamma_ij(p(1),p(2),p(3),l)
          memcheck(3) = adata%gamma_ji(p(1),p(2),p(3),l)
          !memcheck(2) = memcheck(3) / memcheck(2)
          !memcheck(3) = adata%n(p(1),p(2),p(3),atom%line(l)%i) / adata%n(p(1),p(2),p(3),atom%line(l)%j)
          if (l<=atom%nline) then
            memcheck(4) = adata%wphi(p(1),p(2),p(3),l)
          endif
        endif
      else
        memcheck(2) = 10.0**(atom%abnd-12.0)
        memcheck(3) = comp%toth(p(1),p(2),p(3))
      endif
    endif

    if (self%is_finished() .or. self%iter==self%maxiter) then
      do i=1, self%natom
        call os%deallocate(self%atom(i)%data%wphi, label='wphi')
      end do
    endif
    call trace%end(itimer)
  END SUBROUTINE execute

  !===============================================================================
  SUBROUTINE finish (self)
    class(atom_patch_t):: self
    class(atom_type), pointer :: atom
    class(atom_data_t), pointer :: adata
    class(composition_type), pointer :: comp
    real(4), dimension(:,:,:,:,:,:), allocatable :: out_mem
    integer :: i, j, l
    integer, save :: count=0
    !-----------------------------------------------------------------------------
    if (self%is_finished()) return
    l = 1
    do i=1, self%natom
      atom => self%atom(i)
      adata => atom%data
      call self%project_populations(atom, out_mem)
      if (save_populations) then
        call self%save6(out_mem, label='pop_'//trim(atom%id))
      endif
      deallocate(out_mem)
      do j=1, atom%nlevel
        rank_lvlerr(l) = merge(rank_lvlerr(l), adata%lvlerr(j), abs(rank_lvlerr(l)) > abs(adata%lvlerr(j)))
        l = l +1
      end do
    end do
    call self%m3d_patch_t%finish()
    count = count + 1
    if (count == self%task_list%dims(3)) then
      call self%save_convergence_report()
    endif
  END SUBROUTINE finish

  !====================================================================================
  SUBROUTINE project_populations(self, atom, out_mem)
    class(atom_patch_t), intent(inout) :: self
    class(atom_type), intent(in) :: atom
    class(composition_type), pointer :: comp
    class(atom_data_t), pointer :: adata
    real(4), intent(out), dimension(:,:,:,:,:,:), allocatable :: out_mem
    real(8), dimension(:,:,:), allocatable :: norm, totn
    real(8) :: vec(3), tmp
    integer :: i,j,k,l, n(3), imu, nmu
    logical :: iterated
    !---------------------------------------------------------------------------------
    adata => atom%data
    n = shape(adata%totn)
    nmu = merge(size(atom%rotated_data), 0, associated(atom%rotated_data))
    allocate(out_mem(n(1), n(2), n(3), atom%nlevel, 2, 1+nmu))

    if (nmu > 0) then
      allocate(norm, totn, mold=adata%totn)
    endif

    do l=1, atom%nlevel
    do k=1, n(3)
    do j=1, n(2)
    do i=1, n(1)
      tmp = adata%nstar(i,j,k,l)
      tmp = merge(ieee_value(tmp, ieee_quiet_nan), tmp, tmp==0)
      out_mem(i,j,k,l,1,1) = real(adata%nstar(i,j,k,l))
      out_mem(i,j,k,l,2,1) = real(adata%n(i,j,k,l) / tmp)
    end do
    end do
    end do
    end do

    do imu=1, nmu
      vec = long_vec(:,imu)
      adata => atom%rotated_data(imu)
      comp => self%m3d_parent%rotated_composition(imu)
      iterated = .not. associated(adata%n, adata%nstar)
      call atom%init_atom_data(comp, adata, iterated)
      if (iterated) then
        totn = 0
        do l=1, atom%nlevel
          do k=1, n(3)
          do j=1, n(2)
          do i=1, n(1)
            tmp = adata%nstar(i,j,k,l)
            tmp = merge(ieee_value(tmp, ieee_quiet_nan), tmp, tmp==0)
            adata%n(i,j,k,l) = adata%n(i,j,k,l) / tmp
          end do
          end do
          end do
          call self%rotate_cube8(vec, adata%n(:,:,:,l), do_log=.true.)
          do k=1, n(3)
          do j=1, n(2)
          do i=1, n(1)
            adata%n(i,j,k,l) = adata%n(i,j,k,l) * adata%nstar(i,j,k,l)
            totn(i,j,k) = totn(i,j,k) + adata%n(i,j,k,l)
          end do
          end do
          end do
        end do
        do k=1, n(3)
        do j=1, n(2)
        do i=1, n(1)
          norm(i,j,k) = 1 / totn(i,j,k)
        end do
        end do
        end do
        do l=1, atom%nlevel
          do k=1, n(3)
          do j=1, n(2)
          do i=1, n(1)
            tmp = adata%n(i,j,k,l) * adata%totn(i,j,k) * norm(i,j,k)
            adata%n(i,j,k,l) = merge(tmp, 0.d0, tmp==tmp)
            tmp = adata%nstar(i,j,k,l)
            tmp = merge(ieee_value(tmp, ieee_quiet_nan), tmp, tmp==0)
            out_mem(i,j,k,l,1,1+imu) = real(adata%nstar(i,j,k,l))
            out_mem(i,j,k,l,2,1+imu) = real(adata%n(i,j,k,l) / tmp)
          end do
          end do
          end do
        end do
        call atom%compute_cont_factors(adata)
      endif
    end do
  END SUBROUTINE project_populations

  !===============================================================================
  SUBROUTINE save_convergence_report (self)
    class(atom_patch_t) :: self
    integer :: lunit, i
    character(len=256) :: conv_file
    integer, save :: itimer=0
    logical, save :: first_time = .true.
    !-----------------------------------------------------------------------------
    if (first_time) then
      first_time = .false.
      call trace%begin('atom_patch_t%save_convergence_report', itimer=itimer)
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
    class(atom_patch_t):: self
    class(atom_data_t), pointer :: adata
    character(len=64) :: text
    integer :: i, ipos
    real(8) :: global
    !-----------------------------------------------------------------------------
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
    global = global_maxerr%update(0.d0)
    ipos = global_pos%update(+0)
    write(stdout,'(A, ES12.4, A10, ES12.4, A, I4, A, "[", 3I4, "] ", A)') &
      '  local =', self%local,  ', global =', global, &
      ', lvl =',elvl, ', pos=', epos, self%atom(eatom)%id
    if (self%overrelax) then
      write(stdout,'(A, 3(ES12.4))') '  omega =', self%sor_w, self%prev_prev_global, self%prev_global
    endif

    if (self%verbose > 1) then
      write(stdout,'(1X, A80)') &
      '................................................................................'
    endif
  END SUBROUTINE print_post_update_info

  !===============================================================================
  SUBROUTINE get_osmet_nu (self, sp)
    class(atom_patch_t), target :: self
    class(atom_type), pointer :: atom
    class(atomic_spectrum_type) :: sp
    integer :: i, j, fs, fe
    !-----------------------------------------------------------------------------
    if (.not. osmet%initialised) return
    do j=1, self%natom
      atom => self%atom(j)
      do i=1, atom%ncont
        associate(cont => atom%trans%cont(i))
          fs = searchsorted(osmet%nu, cont%nu_min, outbounds_ok=.true.)+1
          fe = searchsorted(osmet%nu, cont%nu_max, outbounds_ok=.true.)
          osmet%active(fs:fe) = .true.
        end associate
      end do
    end do
    sp%nnu = count(osmet%active)
    allocate(sp%nu(sp%nnu))
    sp%nu = pack(osmet%nu, osmet%active)
  END SUBROUTINE get_osmet_nu

  !===============================================================================
  SUBROUTINE exclude_double_continua (self, atom)
    class(atom_patch_t), target :: self
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
      case ('(H+H)', 'H2+', 'H- FF', 'H FF', 'H2-', 'H- BF')
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
                print *, 'excluding transition ', trans%label
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
  SUBROUTINE exclude_double_lines (self, atom, species)
    class(atom_patch_t), target :: self
    type(atom_type) :: atom
    class(massive_type), pointer :: species
    class(element_type), pointer :: element
    class(molec_type), pointer :: molec
    class(ion_type), pointer :: ion
    class(bbtrans_type), pointer :: ll_line, line
    integer :: i, j, l
    logical, save :: first_time=.true.
    !-----------------------------------------------------------------------------
    select type (species)
    class is (element_type)
      element => species
      if (all(element%ions%nline == 0)) return
      if (exclude_from_line_list) then
        element%ions%in_linelist = .false.
      else
        do j=lbound(element%ions,1), ubound(element%ions,1) - 1
          ion => element%ions(j)
          if (ion%nline == 0) cycle
  
          do l=1, ion%nline
            ll_line => ion%line(l)
            
            do i=1, atom%nline
              line => atom%trans%line(i)
              if (abs(line%nu0 / ll_line%nu0 -1) < 1e-6 &
                  .and. abs(line%Aji / ll_line%Aji - 1) < 1e-4) then
                call ll_line%deactivate()
                print '(A, F10.2)', 'excluding linelist line at ', CC / ll_line%nu0 *1e8
              endif
            end do
          end do
        end do
      endif
    class is (molec_type)
      molec => species
      if (molec%nline == 0) return
      if (exclude_from_line_list) molec%in_linelist = .false.
    end select
    if (exclude_from_line_list) then
        print *, 'excluding all ', trim(atom%id), ' lines from line list'
    endif
  END SUBROUTINE exclude_double_lines
  
  !===============================================================================
  SUBROUTINE atom_patch_nbor_relations(self, another, needed, needs_me, download)
    class(atom_patch_t) :: self
    class(task_t), pointer :: another         ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                             !    does it need me, should I download?
    ! ----------------------------------------------------------------------------
    call trace%begin ('atom_patch_t%nbor_relations')
    call self%m3d_patch_t%nbor_relations(another, needed, needs_me, download)
  
    if (needed .or. needs_me) then
      select type (another)
      class is (nlte_patch_t)
        ! is needed and needs_me
      class default
        if (another%id == self%m3d_parent%id) then
          ! the composition patch needs this to be able to
          ! update again before final spectrum is printed
        else
          needed = .false. ; needs_me = .false.
        endif
      end select
    endif
  
    call trace%end()
  END SUBROUTINE atom_patch_nbor_relations

  !===============================================================================
  FUNCTION check_in_bounds (self) result (cib)
    class(atom_patch_t) :: self
    logical :: cib
    !-----------------------------------------------------------------------------
    cib = self%m3d_patch_t%check_in_bounds()
    cib = cib .and. lvlcheck > 0 .and. lvlcheck <= self%atom(atomcheck)%nlevel
  END FUNCTION check_in_bounds

  !===============================================================================
  SUBROUTINE init_nbors_hash (self, nbors)
    class(atom_patch_t):: self
    class(nbor_t), pointer:: nbors
    integer, save:: itimer=0
    integer :: key(4), i_nu, i_mu
    !-----------------------------------------------------------------------------
    call self%m3d_patch_t%init_nbors_hash(nbors)
    key = self%key
    key(4) = hash_string('nlte_patch')
    do i_mu=-1, 2
      do i_nu=-self%n_nu, 2*self%n_nu
        key(2) = i_nu
        key(3) = i_mu
        call self%check_nbor(key, nbors)
      end do
    end do
  END SUBROUTINE init_nbors_hash

  !===============================================================================
  SUBROUTINE save_abundance_to_file (self, atom)
    class(atom_patch_t), target :: self
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
  SUBROUTINE allocate_mem(self, m_in)
    class(atom_patch_t) :: self
    integer, intent(in), optional :: m_in(6)
    integer :: i, m(4), nl
    ! ----------------------------------------------------------------------------
    ! This patch has no memory data to be allocated. Everything is in atom%data.
    return
  END SUBROUTINE allocate_mem

  !===============================================================================
  SUBROUTINE deallocate_mem(self)
    class(atom_patch_t) :: self
    integer :: i
    ! ----------------------------------------------------------------------------
    return
  END SUBROUTINE deallocate_mem

  !===============================================================================
  !SUBROUTINE load_mem (self, folder, label)
  !  class(atom_patch_t), intent(inout) :: self
  !  character(len=256), intent(in) :: folder
  !  character(len=*), intent(in) :: label
  !  character(len=256) :: filename, sfolder, metafile, file
  !  real(4), dimension(:,:,:,:,:), allocatable :: mem
  !  integer :: lunit, file_size, mem_size, iostat, prec, oldm(5), newm(5), iv, iw, k, i
  !  logical, save :: first_time=.true.
  !  ! ----------------------------------------------------------------------------
  !  write(sfolder, '(A, A)') trim(folder), '/save/'
  !  write(filename, '(A, A, 2(I3.3), A)') trim(label), '_', self%key(1:2), '.bin'
  !  write(metafile, '(A, A)') trim(label), '_meta.txt'
  !  call os%check_file(trim(sfolder)//filename)
  !  call os%check_file(trim(sfolder)//metafile)
  !  open(newunit=lunit, file=trim(sfolder)//metafile, form='formatted', iostat=iostat)
  !  read(lunit,*) ! skip header
  !  do while (iostat==0)
  !    read(lunit, *, iostat=iostat) file, prec, oldm
  !    if (file == filename) exit
  !  end do
  !  close(lunit)
  !  if (iostat /= 0) then
  !    print *, 'file not found in metafile! Stopping...'
  !    stop
  !  endif
  !  newm = [self%n, self%atom%nlevel, 2]
  !  open(newunit=lunit, file=trim(sfolder)//filename, access='stream')
  !  if (all(oldm == newm)) then
  !    allocate(mem(oldm(1), oldm(2), oldm(3), oldm(4), oldm(5)))
  !    read(lunit) mem
  !    do i=1, self%atom%nlevel
  !      self%atom%data%n(:,:,:,i) = dble(mem(:,:,:,i,2)) * self%atom%data%nstar(:,:,:,i)
  !      where(self%atom%data%n(:,:,:,i) /= self%atom%data%n(:,:,:,i)) self%atom%data%n(:,:,:,i) = 0
  !    end do
  !  !elseif (all(oldm(3:5) == newm(3:5))) then
  !  !  if (first_time) then
  !  !    first_time = .false.
  !  !    write(stdout, '(A, 2(I6,"x",I4))') 'Interpolating patch memory when loading', oldm(1:2), newm(1:2)
  !  !  endif
  !  !  call os%allocate(mem2, oldm)
  !  !  read(lunit) mem2
  !  !
  !  !  !mem2(:,:,:,:,2) = mem2(:,:,:,:,2) / mem2(:,:,:,:,1)
  !  !
  !  !  do iw=1, newm(5)
  !  !  do iv=1, newm(4)
  !  !  do  k=1, newm(3)
  !  !    call interp2d_periodic(mem2(:,:,k,iv,iw), mem(:,:,k,iv,iw))
  !  !  end do
  !  !  end do
  !  !  end do
  !  !
  !  !  mem(:,:,:,:,2) = mem(:,:,:,:,2) * mem(:,:,:,:,1)
  !  !
  !  !  call os%deallocate(mem2)
  !  else
  !    write(stdout, *) 'Cannot load patch save ', trim(sfolder)//trim(filename), &
  !    'Memory sizes do not match!'
  !    print *, oldm
  !    print *, newm
  !    stop
  !  endif
  !  close(lunit)
  !END SUBROUTINE load_mem

  !===============================================================================
  SUBROUTINE check_mem_NaN (self)
    class(atom_patch_t):: self
    integer :: i,j,k,v,t,w
    integer, save :: itimer
    !-----------------------------------------------------------------------------
    call trace%begin ('atom_patch_t%check_mem_NaN', itimer=itimer)
    do w=1, self%natom
      associate(n=>self%atom(w)%data%n)
        if (associated(self%atom(w)%data%n)) then
          do v=lbound(n, 4), ubound(n, 4)
          do k=lbound(n, 3), ubound(n, 3)
          do j=lbound(n, 2), ubound(n, 2)
          do i=lbound(n, 1), ubound(n, 1)
            if (n(i,j,k,v) /= n(i,j,k,v)) then
              write(stdout, '(A,A,3(I4),I5,I3,I3)') self%m3d_patch_info(), '  NaN', i,j,k,v
              stop
            elseif (n(i,j,k,v) > inf4) then
              write(stdout, '(A,A,3(I4),I5,I3,I3)') self%m3d_patch_info(), '  Infinity', i,j,k,v
              stop
            endif
          end do
          end do
          end do
          end do
        endif
      end associate
    end do
    call trace%end(itimer)
  END SUBROUTINE check_mem_NaN

  END MODULE
