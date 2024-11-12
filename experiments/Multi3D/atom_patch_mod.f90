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
                        red_cut, blue_cut, logf_cut, atom_data_t, m1d_legacy_mode, &
                        nq_fac, nq_min, nq_max, disable_bf, disable_bb, level_partf
    USE coll_mod, only: cbbesc, cbfesc, cbbhsc, cbfhsc
    USE trans_mod, only: line_type, cont_type, transition_type
    USE const_mod, only: CC
    USE ion_mod, only: ion_type, particle_type, debye_shield
    USE element_mod, only: element_type
    USE molecule_mod, only: molec_type
    USE bbtrans_mod, only: bbtrans_type
    USE composition_mod, only: composition_type
    USE m3d_patch_mod, only: m3d_patch_t, composition, maxiter, check_NaN, cswitch, key_link_t, &
    global_maxerr, convlim, convmax, atom_update_counter, comoving_frame, global_pos, ng_step, ng_start, &
    decouple_continuum, overlapping_bb, overlapping_bf, line_blanketing, converged, premature_finish, &
    rotate_atmos, long_vec, rank_lvlerr, overlapping_bb_in_output, long_scheme, n_nu
    USE bftrans_mod, only: bftrans_type
    USE hprofl_mod, only: hydrogen_BPO
    USE osmet_mod, only:osmet_type
    USE nlte_patch_mod, only:nlte_patch_t, n_cntrbf, cntrbf_lines
    USE m3d_debug_mod
    USE shared_mod
    USE, intrinsic :: ieee_arithmetic
    implicit none
  
    !-----------------------------------------------------------------------------
  
    type, public, extends(m3d_patch_t) :: atom_patch_t
      type(atom_type) :: output_atom, bb_atom, bf_atom
      class(atomic_spectrum_type), pointer :: output_sp=>null()
      class(atomic_spectrum_type), pointer :: bf_sp=>null()
      class(atomic_spectrum_type), pointer :: bb_sp=>null()
      real(8) :: maxerr=huge(1.0)
      real(8), dimension(:), allocatable :: lvlerr
      logical :: update_skipped=.false., atom_allocated=.false., molec
    contains
      procedure :: pre_init => atom_patch_pre_init
      procedure :: execute
      procedure :: finish
      procedure :: nbor_relations => atom_patch_nbor_relations
      procedure :: pre_init_atom
      procedure :: pre_init_spectrum
      procedure :: init_atom_data
      procedure :: update_populations
      procedure :: get_osmet_nu
      procedure :: exclude_double_continua
      procedure :: exclude_double_lines
      procedure :: print_post_update_info
      procedure :: check_in_bounds
      procedure :: init_nbors_hash
      procedure :: save_abundance_to_file
      procedure :: save_convergence_report
      procedure :: load_mem
      procedure :: allocate_mem
      procedure :: deallocate_mem
      procedure :: check_mem_NaN
      procedure :: ng_acceleration
    end type atom_patch_t

    integer, save :: epos(3), elvl
    logical, save :: exclude_from_line_list=.true.
    logical, save :: detailed_absmet=.false.
    logical, save :: pandas_format=.false.
    logical, save :: save_populations=.true.
    logical, save :: save_rates=.false.
    character(len=256), save :: precomputed_depart=''
    character(len=256), public :: atom_file=''
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
    call self%set_nv_and_nw(self%atom%nlevel, 2)
    self%static_mem=.true.

    call self%post_pre_init()
    call trace%end(itimer)
  END SUBROUTINE atom_patch_pre_init

  !===============================================================================
  SUBROUTINE pre_init_atom (self)
    class(atom_patch_t) :: self
    type(atom_type), save, target :: atom, output_atom
    type(line_type), save, dimension(:), pointer, contiguous :: dc_line, out_line
    type(cont_type), save, dimension(:), pointer, contiguous :: dc_cont, out_cont
    real(8), save :: qnorm=10.0, abundance=-99.9
    integer, save :: itimer=0, n_in_lines=0
    integer, save, allocatable :: in_lines(:)
    logical, save :: first_time=.true., exclude_trace_cont=.true., use_atom_abnd=.true.
    integer :: i, iostat
    namelist /atom_params/ atom_file, qnorm, convlim, exclude_trace_cont, level_partf, &
    use_atom_abnd, logf_cut, red_cut, blue_cut, n_in_lines, abundance, &
    exclude_from_line_list, precomputed_depart, use_classic_col, kramer_bf, overlapping_bf, &
    overlapping_bb, hydrogen_BPO, detailed_absmet, n_cntrbf, convmax, premature_finish, &
    pandas_format, save_populations, nq_fac, nq_min, nq_max, cbbesc, cbfesc, cbbhsc, cbfhsc, &
    disable_bf, disable_bb, overlapping_bb_in_output, save_rates
    namelist /line_mask/ in_lines, cntrbf_lines
    !-----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%pre_init_atom', itimer=itimer)
    if (first_time) then
      first_time = .false.

      rewind (io%input)
      read (io%input, atom_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, atom_params)
      end if
      call io%namelist_warning ('atom_params', iostat)

      if (atom_file == '') then
        print *, 'WARNING: no atom specified'
        call trace%end(itimer)
        return
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

      write (*,*) ' - Reading atom'
      if (pandas_format) then
        call atom%read_pandas(atom_file, qnorm, in_lines=in_lines)
      else
        call atom%read(atom_file, qnorm, in_lines=in_lines, read_coll=.true.)
      endif
      atom%name = 'input'
      atom%verbose = self%verbose
      allocate(rank_lvlerr(atom%nlevel), source=0.d0)

      call os%copy_file(trim(atom_file), trim(io%run_dir)//'atom.txt')
      if (abundance > -99) atom%abnd = abundance

      if (decouple_continuum) then
        allocate(dc_line(atom%nline))
        allocate(dc_cont(atom%ncont))
        dc_line = atom%line
        dc_cont = atom%cont
      endif

      allocate(out_line(atom%nline))
      allocate(out_cont(atom%ncont))
      out_line = atom%line
      out_cont = atom%cont

      write (io%namelist, atom_params)
    endif

    if (atom_file == '') then
      call trace%end(itimer)
      return
    endif

    allocate(self%atom)
    self%atom = atom
    self%atom%name = 'main'
    self%n_nu = n_nu

    ! this pointer will be referenced by other atoms with same local data
    allocate(self%atom%data)
    if (rotate_atmos) then
      allocate(self%atom%rotated_data(size(long_vec, 2)))
    endif

    if (use_atom_abnd .or. abundance > -99) then
      call composition%overwrite_abundance(atom%id, atom%abnd, atom%awgt)
    else
      self%atom%abnd = self%atom%abnd + composition%FeH
    endif
    
    if (exclude_trace_cont) then
      call self%exclude_double_continua()
    endif
    call self%exclude_double_lines()

    if (decouple_continuum) then
      self%bb_atom = self%atom
      self%bf_atom = self%atom
      call self%bb_atom%link_transitions(null(), dc_line)
      call self%bf_atom%link_transitions(dc_cont, null())
      self%bb_atom%name = 'bound-bound'
      self%bf_atom%name = 'bound-free'
    endif

    self%output_atom = self%atom
    call self%output_atom%link_transitions(out_cont, out_line)
    self%output_atom%name = 'output'

    self%atom_allocated = .true.
    call trace%end(itimer)
  END SUBROUTINE pre_init_atom

  !===============================================================================
  SUBROUTINE pre_init_spectrum (self, output_sp)
    class(atom_patch_t) :: self
    type(atomic_spectrum_type), target :: output_sp
    type(atomic_spectrum_type), save, target :: sp, bb_sp, bf_sp
    logical, save :: first_time=.true.
    logical :: save_sp
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%pre_init_spectrum', itimer=itimer)
    if (first_time) then
      first_time = .false.
      external_sp = output_sp%nnu > 0 
      save_sp = external_sp .or. (associated(self%atom) .and. size(long_vec,2)>0)
      if (associated(self%atom)) then
        if (decouple_continuum) then
          bf_sp%name = 'bound-free'
          if (detailed_absmet .and. overlapping_bf) then
             call self%get_osmet_nu(bf_sp, composition%osmet)
          endif
          call self%bf_atom%create_spectrum(bf_sp, bf=.true., sort=overlapping_bf)
          call self%bf_atom%set_ired_iblue(bf_sp, bf=.true., set_active=.true.)
          if (save_sp) call bf_sp%save('bf')

          bb_sp%name = 'bound-bound'
          call self%bb_atom%create_spectrum(bb_sp, bb=.true., sort=overlapping_bb)
          call self%bb_atom%set_ired_iblue(bb_sp, bb=.true., set_active=.true.)
          if (save_sp) call bb_sp%save('bb')

        else
          sp%name = 'all'
          if (detailed_absmet) call self%get_osmet_nu(sp, composition%osmet)
          call self%atom%create_spectrum(sp, bb=.true., bf=.true.)
          call self%atom%set_ired_iblue(sp, bb=.true., bf=.true., set_active=.true.)
          if (save_sp) call sp%save('all')
        endif

        if (external_sp) then
          call self%output_atom%set_ired_iblue(output_sp, bb=.true., bf=.true.) 
        else
          output_sp%name = 'output based on atom'
          call self%output_atom%create_spectrum(output_sp, bb=.true., sort=overlapping_bb_in_output)
          call self%output_atom%set_ired_iblue(output_sp, bb=.true., bf=.true.) 
        endif
      elseif (external_sp) then
        call self%output_atom%set_ired_iblue(output_sp, bb=.true., bf=.true.) 
      !else
      !  stop 'Neither atom nor spectral range specified!'
      endif
      if (save_sp) then
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
  END SUBROUTINE pre_init_spectrum

  !===============================================================================
  SUBROUTINE init_atom_data (self, comp, adata)
    class(atom_patch_t), target :: self
    class(atom_type), pointer :: atom
    class(composition_type), target :: comp
    class(atom_data_t), target :: adata
    class(element_type), pointer :: element
    class(molec_type), pointer :: molec
    real(8), dimension(:,:,:), pointer, contiguous :: dxi
    real(4) :: nmax, nlim
    integer, save:: itimer=0
    integer :: i, m(4)
    !-----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%init_atom_data', itimer=itimer)
    atom => self%atom
    dxi => comp%dxi

    m = [shape(comp%temp), atom%nlevel]
    if (.not. associated(adata%nstar)) call os%allocate(adata%nstar, m, label='nstar')
    if (.not. associated(adata%totn))  call os%allocate(adata%totn, m(1:3), label='totn')
    if (.not. associated(adata%n))     call os%allocate(adata%n, m, label='n')

    call comp%get_totn(atom%id, atom%abnd, adata%totn, element, molec)
    if (.not. (associated(element) .or. associated(molec))) then
      call self%save_abundance_to_file()
    endif

    if (associated(molec)) atom%tmax_update = comp%tmolim
    if (associated(element)) call atom%check_ionisation_thresholds(element)
    call atom%compute_lte_populations(adata%nstar, adata%totn, real(comp%temp), real(comp%saha), real(dxi), element, molec)
    call trace%end(itimer)
  END SUBROUTINE init_atom_data

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
    real(8), dimension(:,:,:,:), pointer, contiguous :: gamma_ij, gamma_ji, nstar
    real(4), dimension(:,:,:,:), pointer, contiguous :: wphi
    integer :: i, j, l, p(3), ii, jj, m(4)
    integer, save:: count=0
    real(8) :: global, gij, vec(3)
    logical :: do_ng
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin (trim(self%kind)//'%update', itimer=itimer)
    atom => self%atom
    adata => atom%data

    if (self%iter == 0) then
      self%composition => self%m3d_parent%composition
      call self%init_atom_data(self%composition, adata)
      if (precomputed_depart /= '') then
        call self%load_mem(precomputed_depart, self%kind)
      elseif (maxiter > 0) then
        do i=1, atom%nlevel
          adata%n(:,:,:,i) = adata%nstar(:,:,:,i) 
        end do
      else
        call os%deallocate(adata%n, 'n')
        adata%n => adata%nstar
      endif
    else
      wphi => adata%wphi
      nstar => adata%nstar
      gamma_ij => adata%gamma_ij
      gamma_ji => adata%gamma_ji

      if (self%iter == 1) then
        !self%dtime = self%dt_fixed

        ! Normalize wphi
        do i=1, atom%nline
          if (.not. atom%line(i)%active) cycle
          wphi(:,:,:,i) = 1 / wphi(:,:,:,i)
        end do
      endif

      if (adata%converged) then
        self%update_skipped = .true.
      else
        do i=1, atom%nline
          line => atom%line(i)
          if (.not. line%active) cycle
          gij = atom%g(line%i) / atom%g(line%j)
          gamma_ij(:,:,:,i) = gamma_ij(:,:,:,i) * wphi(:,:,:,i)
          gamma_ji(:,:,:,i) = gamma_ji(:,:,:,i) * wphi(:,:,:,i)
          gamma_ji(:,:,:,i) = gamma_ji(:,:,:,i) * gij
        end do
        do i=i, atom%nrad
          cont => atom%cont(i-atom%nline)
          if (.not. cont%active) cycle
          gamma_ji(:,:,:,i) = gamma_ji(:,:,:,i) * nstar(:,:,:,cont%i) / nstar(:,:,:,cont%j)
        end do

        if (ng_step > 0) then
          do i=1, atom%nlevel
            adata%oldn(:,:,:,i,1) = adata%oldn(:,:,:,i,2)
            adata%oldn(:,:,:,i,2) = adata%oldn(:,:,:,i,3)
            adata%oldn(:,:,:,i,3) = adata%n(:,:,:,i)
          end do
        endif

        call self%update_populations(gamma_ij, gamma_ji)
        adata%converged = self%maxerr < convmax
        self%update_skipped = .false.
        do_ng = ng_step > 0 .and. mod(self%iter-ng_start, ng_step) == 0 .and. self%iter >= ng_start
        if (adata%converged .and. premature_finish) then
          call self%finish()
        elseif (do_ng) then
          call self%ng_acceleration()
        endif
      endif
    endif

    if (.not. self%update_skipped) then
      if (.not. (maxiter==0 .and. long_scheme=='none')) call atom%compute_cont_factors(atom%data)
    endif
    
    if (self%iter > 0) then
      count = atom_update_counter%update(+1)
    endif
    if (count >= self%task_list%dims(3)) then
      global = global_maxerr%update(0.d0)
      if ((global > 0) .and. global < convlim .and. self%maxerr < convlim) then
        !write(stdout, *) self%ipos(3)+1, 'Converged!', global
        if (self%iter < self%maxiter) converged = .true.
        call self%finish()
      endif
    endif

    if (self%print_update_info) then
      p = self%pcheck
      l = self%linecheck
      memcheck(1) = adata%n(p(1),p(2),p(3),self%lvlcheck)
      if (associated(adata%gamma_ij)) then
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
        memcheck(3) = self%composition%toth(p(1),p(2),p(3))
      endif
    endif

    if (associated(adata%gamma_ij)) then
      call os%deallocate(adata%gamma_ij, label='atom%gamma_ij')
      call os%deallocate(adata%gamma_ji, label='atom%gamma_ji')
    endif

    if (self%is_finished() .or. self%iter==self%maxiter) then
      if (associated(adata%wphi)) then
        call os%deallocate(adata%wphi, 'wphi')
      endif
    endif
    call trace%end(itimer)
  END SUBROUTINE execute

  !===============================================================================
  SUBROUTINE ng_acceleration (self)
    class(atom_patch_t):: self
    integer :: i, x, y, z, k, epos(4), nvec
    real(8), dimension(:), allocatable :: n1, n2, n3, n4, n5, d0, d1, d2, e
    real(8) :: d01, d02, d11, d12, d21, d22, a, b, den, emax, enew
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin ('atom_patch_t%ng_acceleration', itimer=itimer)
    nvec = self%atom%nlevel * product(self%n(1:2))
    allocate(e(nvec))
    allocate(n1, n2, n3, n4, n5, d0, d1, d2, mold=e)
    emax = 1e-3
    do z=1, self%n(3)
      i = 0
      do y=1, self%n(2)
      do x=1, self%n(1)
      do k=1, self%atom%nlevel
        i = i+1
        n4(i) = self%atom%data%n(x,y,z,k)
        n1(i) = self%atom%data%oldn(x,y,z,k,1) / n4(i)
        n2(i) = self%atom%data%oldn(x,y,z,k,2) / n4(i)
        n3(i) = self%atom%data%oldn(x,y,z,k,3) / n4(i)
      end do
      end do
      end do

      d0 = 1. - n3
      d1 = n3 - n2
      d2 = n2 - n1

      d1 = d0 - d1
      d2 = d0 - d2

      d01 = dot_product(d0, d1)
      d02 = dot_product(d0, d2)
      d11 = dot_product(d1, d1)
      d22 = dot_product(d2, d2)
      d12 = dot_product(d1, d2)
      d21 = d12

      den = d11*d22 - d12*d21

      if (abs(den) > 0) then
        a = (d01*d22 - d02*d21) / den
        b = (d02*d11 - d01*d21) / den

        e = 1.d0 - a - b + a*n3 + b*n2
      else
        e = 1.d0
      endif

      e = max(e, 0.001d0)

      i = 0
      do y=1, self%n(2)
      do x=1, self%n(1)
      do k=1,self%atom%nlevel
        i = i+1
        self%atom%data%n(x,y,z,k) = e(i)*n4(i)
        enew = e(i) - 1.d0

        if (abs(enew) > emax) then
          emax = enew
          epos = [k, x, y, z]
        endif
      end do
      end do
      end do
    end do
    print '(A, ES14.7, 4(I6))', 'NG acceleration ', emax, epos
    call trace%end(itimer)
  END SUBROUTINE ng_acceleration

  !===============================================================================
  SUBROUTINE finish (self)
    class(atom_patch_t):: self
    class(atom_type), pointer :: atom
    class(atom_data_t), pointer :: adata
    class(composition_type), pointer :: comp
    real(4), dimension(:,:,:,:,:,:), allocatable :: out_mem
    real(8), dimension(:,:,:), allocatable :: norm, totn
    real(8) :: vec(3), tmp
    integer :: i, j, k, l, imu, n(3)
    integer, save :: count=0
    !-----------------------------------------------------------------------------
    atom => self%atom
    n = self%n

    if (rotate_atmos) then
      allocate(norm, totn, mold=atom%data%totn)
      allocate(out_mem(n(1), n(2), n(3), atom%nlevel, 2, 1+self%n_out))
    else
      allocate(out_mem(n(1), n(2), n(3), atom%nlevel, 2, 1))
    endif

    adata => atom%data
    do l=1, atom%nlevel
    do k=1, self%n(3)
    do j=1, self%n(2)
    do i=1, self%n(1)
      tmp = adata%nstar(i,j,k,l)
      tmp = merge(ieee_value(tmp, ieee_quiet_nan), tmp, tmp==0)
      out_mem(i,j,k,l,1,1) = real(adata%nstar(i,j,k,l))
      out_mem(i,j,k,l,2,1) = real(adata%n(i,j,k,l) / tmp)
    end do
    end do
    end do
    end do

    if (rotate_atmos) then
      do imu=1, self%n_out
        vec = long_vec(:,imu)
        adata => atom%rotated_data(imu)
        comp => self%m3d_parent%rotated_composition(imu)
        call self%init_atom_data(comp, adata)
        if (associated(self%atom%data%n, self%atom%data%nstar)) then
          call os%deallocate(adata%n, 'n')
          adata%n => adata%nstar
        else
          totn = 0
          do l=1, atom%nlevel
            do k=1, self%n(3)
            do j=1, self%n(2)
            do i=1, self%n(1)
              tmp = atom%data%nstar(i,j,k,l)
              tmp = merge(ieee_value(tmp, ieee_quiet_nan), tmp, tmp==0)
              adata%n(i,j,k,l) = atom%data%n(i,j,k,l) / tmp
            end do
            end do
            end do
            call self%rotate_cube8(vec, adata%n(:,:,:,l), do_log=.true.)
            do k=1, self%n(3)
            do j=1, self%n(2)
            do i=1, self%n(1)
              adata%n(i,j,k,l) = adata%n(i,j,k,l) * adata%nstar(i,j,k,l)
              totn(i,j,k) = totn(i,j,k) + adata%n(i,j,k,l)
            end do
            end do
            end do
          end do
          do k=1, self%n(3)
          do j=1, self%n(2)
          do i=1, self%n(1)
            norm(i,j,k) = 1 / totn(i,j,k)
          end do
          end do
          end do
          do l=1, atom%nlevel
            do k=1, self%n(3)
            do j=1, self%n(2)
            do i=1, self%n(1)
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
        endif
        if (.not. (maxiter==0 .and. long_scheme=='none')) call atom%compute_cont_factors(adata)
      end do
    endif

    if (save_populations) then
      call self%save6(out_mem, label='atom_patch')
    endif
    call self%m3d_patch_t%finish()

    do i=1, atom%nlevel
      rank_lvlerr(i) = max(rank_lvlerr(i), self%lvlerr(i))
    end do
    count = count + 1
    if (count == self%task_list%dims(3)) then
      call self%save_convergence_report()
    endif
  END SUBROUTINE finish

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
    integer :: i, l, p(3), ipos
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
    write(stdout,'(A, ES12.4, A10, ES12.4, A, I4, A, "[", 3I4, "]")') &
      '  local =', self%maxerr,  ', global =', global, &
      ', lvl =',elvl, ', pos=', epos
    write(stdout,'(1X, A80)') &
    '................................................................................'
  END SUBROUTINE print_post_update_info

  !===============================================================================
  SUBROUTINE update_populations (self, gamma_ij, gamma_ji)
    class(atom_patch_t), target :: self
    class(atom_type), pointer :: atom
    class(transition_type), pointer :: trans
    integer, dimension(3) :: ep
    real(8), dimension(:), pointer :: Gij, Gji
    real(8), dimension(:), pointer :: n, nstar
    real(8), dimension(:,:,:,:), target :: gamma_ij, gamma_ji
    real(4), dimension(:,:,:), allocatable :: cij, cji
    integer :: x, y, z, i, j, ii, ipos, il, nx, ny, nz, it, jt
    real(8) :: temp, ne, global, totn, nh1
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%update_populations', itimer=itimer)
    atom => self%atom
    self%lvlerr = 0.0d0
    self%maxerr = 0.0d0
    nx = self%n(1)
    ny = self%n(2)
    nz = self%n(3)

    if (save_rates) then
      allocate(cij(nx,ny,nz), source=0.0)
      allocate(cji(nx,ny,nz), source=0.0)
      if (linecheck <= atom%nline) then
        trans => atom%line(linecheck)
      else
        trans => atom%cont(linecheck-atom%nline)
      endif
      i = trans%i
      j = trans%j
    endif

    do z=1, nz
      do y=1, ny
        do x=1, nx
          ne    = self%composition%ne(x,y,z)
          temp  = self%composition%temp(x,y,z)
          nh1   = self%composition%elements(1)%ions(1)%levels(1)%n(x,y,z)
          totn  =  atom%data%totn  (x,y,z)
          n     => atom%data%n     (x,y,z,:)
          nstar => atom%data%nstar (x,y,z,:)
          Gij   => gamma_ij   (x,y,z,:)
          Gji   => gamma_ji   (x,y,z,:)

          if (totn > 0) then
            call atom%solve_stat_equilibrium(temp, ne, totn, n, nstar, nh1, Gij, Gji, cswitch, check_NaN, z)

            ! keep track of maximum population change
            do ii=1,atom%nlevel
              if (abs(atom%data%e(ii)) > self%lvlerr(ii)) then
                  self%lvlerr(ii) = abs(atom%data%e(ii))
                  if (self%lvlerr(ii) > self%maxerr) then
                    self%maxerr = self%lvlerr(ii)
                    ep = [x, y, z + self%zs]
                    il = ii
                  endif
              end if
            end do

            if (save_rates) then
              cij(x,y,z) = atom%data%c(j, i)
              cji(x,y,z) = atom%data%c(i, j)
            endif
          endif
        end do
      end do
    end do

    if (save_rates) then
      call self%save3(real(gamma_ij(:,:,:,linecheck)), label='rij')
      call self%save3(real(gamma_ji(:,:,:,linecheck)), label='rji')
      call self%save3(cij, label='cij')
      call self%save3(cji, label='cji')
    endif

    if (self%maxerr /= self%maxerr) then
      call shared%abort('NaN in convergence error')
    endif

    !$omp critical (maxerr_reset)
    global = global_maxerr%update(0.d0)
    ipos = global_pos%update(+0)
    if (ipos == self%ipos(3) .or. self%maxerr > global) then
      if (ipos /= self%ipos(3)) call global_pos%reset(self%ipos(3))
      call global_maxerr%reset(self%maxerr)
      epos = ep
      elvl = il
    endif
    !$omp end critical (maxerr_reset)
    call trace%end(itimer)
  END SUBROUTINE update_populations

  !===============================================================================
  SUBROUTINE get_osmet_nu (self, sp, osmet)
    class(atom_patch_t), target :: self
    class(atomic_spectrum_type) :: sp
    class(osmet_type) :: osmet
    integer :: i, fs, fe
    !-----------------------------------------------------------------------------
    if (.not. osmet%initialised) return
    do i=1, self%atom%ncont
      associate(cont => self%atom%cont(i))
        fs = searchsorted(osmet%nu, cont%nu_min, outbounds_ok=.true.)+1
        fe = searchsorted(osmet%nu, cont%nu_max, outbounds_ok=.true.)
        osmet%active(fs:fe) = .true.
      end associate
    end do
    sp%nnu = count(osmet%active)
    allocate(sp%nu(sp%nnu))
    sp%nu = pack(osmet%nu, osmet%active)
  END SUBROUTINE get_osmet_nu

  !===============================================================================
  SUBROUTINE exclude_double_continua (self)
    class(atom_patch_t), target :: self
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
          if (uppercase(trans%ion%id(1:2)) == self%atom%id) then
            if (trans%molecular) then
              ok = .true.
            else
              do kc=1, self%atom%ncont
                if (.not. self%atom%cont(kc)%active) cycle
                lower = self%atom%cont(kc)%i

                if (self%atom%ev(lower)==0) then
                  ok = trans%level%xion == self%atom%ev(lower)
                else
                  ok = abs(trans%level%xion / self%atom%ev(lower) - 1) < 1d-2
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
  SUBROUTINE exclude_double_lines (self)
    class(atom_patch_t), target :: self
    class(element_type), pointer :: element
    class(molec_type), pointer :: molec
    class(ion_type), pointer :: ion
    class(bbtrans_type), pointer :: ll_line, line
    integer :: i, j, l
    logical, save :: first_time=.true.
    !-----------------------------------------------------------------------------
    call composition%element_by_id(trim(self%atom%id), element, verbose=.false.)
    if (.not. associated(element)) then
      call composition%molecule_by_id(trim(self%atom%id), molec, verbose=.false.)
      if (.not. associated(molec)) then
        return
      else
        if (molec%nline == 0) return
      endif
    else
      if (all(element%ions%nline == 0)) return
    endif

    if (exclude_from_line_list) then
      if (associated(element)) then
        element%ions%in_linelist = .false.
      else
        molec%in_linelist = .false.
      endif
      if (first_time) then
        first_time = .false.
        print *, 'excluding all ', trim(self%atom%id), ' lines from line list'
      endif
    else
      do j=lbound(element%ions,1), ubound(element%ions,1) - 1
        ion => element%ions(j)
        if (ion%nline == 0) cycle

        do l=1, ion%nline
          ll_line => ion%line(l)
          
          do i=1, self%atom%nline
            line => self%atom%line(i)
            if (abs(line%nu0 / ll_line%nu0 -1) < 1e-6 &
                .and. abs(line%Aji / ll_line%Aji - 1) < 1e-4) then
              call ll_line%deactivate()
              print '(A, F10.2)', 'excluding linelist line at ', CC / ll_line%nu0 *1e8
            endif
          end do
        end do
      end do
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
    cib = cib .and. self%lvlcheck > 0 .and. self%lvlcheck <= self%atom%nlevel
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
  SUBROUTINE save_abundance_to_file (self)
    class(atom_patch_t), target :: self
    integer :: lunit, i
    logical :: first_time=.true.
    !-----------------------------------------------------------------------------
    if (first_time) then
      first_time = .false.
      open(newunit=lunit, file=trim(io%run_dir)//'abund.txt', form='FORMATTED', access='append')
      write(lunit, '(A3, F8.4)') self%atom%id, self%atom%abnd
      close(lunit)
    endif
  END SUBROUTINE save_abundance_to_file

  !===============================================================================
  SUBROUTINE allocate_mem(self, m_in)
    class(atom_patch_t) :: self
    integer, intent(in), optional :: m_in(6)
    integer :: m(4)
    ! ----------------------------------------------------------------------------
    associate(atom=>self%atom)
      allocate(atom%data%C(atom%nlevel, atom%nlevel))
      allocate(atom%data%W(atom%nlevel, atom%nlevel))
      allocate(atom%data%e(atom%nlevel))
      allocate(self%lvlerr(atom%nlevel), source=0.d0)
      call os%allocate(atom%data%nstar, [self%n, atom%nlevel], label='nstar')
      call os%allocate(atom%data%n,     [self%n, atom%nlevel], label='n')
      call os%allocate(atom%data%totn, self%n, label='totn')
      if (self%maxiter > 0) then
        call os%allocate(atom%data%wphi,  [self%n, atom%nline],  label='wphi', zero=.true.)
        if (ng_step > 0) then
          call os%allocate(atom%data%oldn,  [self%n, atom%nlevel, 3],  label='oldn')
        endif
      endif
    end associate
  END SUBROUTINE allocate_mem

  !===============================================================================
  SUBROUTINE deallocate_mem(self)
    class(atom_patch_t) :: self
    ! ----------------------------------------------------------------------------
    associate(atom=>self%atom)
      deallocate(atom%data%C)
      deallocate(atom%data%W)
      deallocate(atom%data%e)
      deallocate(atom%data%oldn)
      deallocate(self%lvlerr)
      call os%deallocate(self%atom%data%nstar, label='nstar')
      call os%deallocate(self%atom%data%n,     label='n')
      call os%deallocate(atom%data%totn,       label='totn')
      if (self%maxiter > 0) then
        call os%deallocate(atom%data%wphi, label='wphi')
        if (ng_step > 0) then
          call os%deallocate(atom%data%oldn, label='oldn')
        endif
      endif
    end associate
  END SUBROUTINE deallocate_mem

  !===============================================================================
  SUBROUTINE load_mem (self, folder, label)
    class(atom_patch_t), intent(inout) :: self
    character(len=256), intent(in) :: folder
    character(len=*), intent(in) :: label
    character(len=256) :: filename, sfolder, metafile, file
    real(4), dimension(:,:,:,:,:), pointer, contiguous :: mem
    integer :: lunit, file_size, mem_size, iostat, prec, oldm(5), newm(5), iv, iw, k, i
    logical, save :: first_time=.true.
    ! ----------------------------------------------------------------------------
    write(sfolder, '(A, A)') trim(folder), '/save/'
    write(filename, '(A, A, 2(I3.3), A)') trim(label), '_', self%key(1:2), '.bin'
    write(metafile, '(A, A)') trim(label), '_meta.txt'
    call os%check_file(trim(sfolder)//filename)
    call os%check_file(trim(sfolder)//metafile)
    open(newunit=lunit, file=trim(sfolder)//metafile, form='formatted', iostat=iostat)
    read(lunit,*) ! skip header
    do while (iostat==0)
      read(lunit, *, iostat=iostat) file, prec, oldm
      if (file == filename) exit
    end do
    close(lunit)
    if (iostat /= 0) then
      print *, 'file not found in metafile! Stopping...'
      stop
    endif
    newm = [self%n, self%atom%nlevel, 2]
    open(newunit=lunit, file=trim(sfolder)//filename, access='stream')
    if (all(oldm == newm)) then
      call os%allocate(mem, newm)
      read(lunit) mem
      do i=1, self%atom%nlevel
        self%atom%data%n(:,:,:,i) = dble(mem(:,:,:,i,2)) * self%atom%data%nstar(:,:,:,i)
        where(self%atom%data%n(:,:,:,i) /= self%atom%data%n(:,:,:,i)) self%atom%data%n(:,:,:,i) = 0
      end do
      call os%deallocate(mem)
    !elseif (all(oldm(3:5) == newm(3:5))) then
    !  if (first_time) then
    !    first_time = .false.
    !    write(stdout, '(A, 2(I6,"x",I4))') 'Interpolating patch memory when loading', oldm(1:2), newm(1:2)
    !  endif
    !  call os%allocate(mem2, oldm)
    !  read(lunit) mem2
    !
    !  !mem2(:,:,:,:,2) = mem2(:,:,:,:,2) / mem2(:,:,:,:,1)
    !
    !  do iw=1, newm(5)
    !  do iv=1, newm(4)
    !  do  k=1, newm(3)
    !    call interp2d_periodic(mem2(:,:,k,iv,iw), mem(:,:,k,iv,iw))
    !  end do
    !  end do
    !  end do
    !
    !  mem(:,:,:,:,2) = mem(:,:,:,:,2) * mem(:,:,:,:,1)
    !
    !  call os%deallocate(mem2)
    else
      write(stdout, *) 'Cannot load patch save ', trim(sfolder)//trim(filename), &
      'Memory sizes do not match!'
      print *, oldm
      print *, newm
      stop
    endif
    close(lunit)
  END SUBROUTINE load_mem

  !===============================================================================
  SUBROUTINE check_mem_NaN (self)
    class(atom_patch_t):: self
    integer :: i,j,k,v,t,w
    integer, save :: itimer
    !-----------------------------------------------------------------------------
    call trace%begin ('atom_patch_t%check_mem_NaN', itimer=itimer)
    associate(n=>self%atom%data%n)
      if (associated(self%atom%data%n)) then
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
    call trace%end(itimer)
  END SUBROUTINE check_mem_NaN

  END MODULE
