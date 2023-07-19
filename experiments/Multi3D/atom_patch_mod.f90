MODULE atom_patch_mod
    USE os_mod
    USE list_mod
    USE halo_mod
    USE gpatch_mod
    USE link_mod
    USE task_mod
    USE omp_timer_mod
    !USE index_mod
    !USE mesh_mod
    USE io_mod
    !USE aux_mod
    USE timer_mod
    USE scaling_mod
    USE trace_mod
    USE mpi_mod
    !USE omp_mod
    !USE omp_timer_mod
    USE math_mod
    USE bits_mod
    !USE units_mod
    !USE scratch_mod
    USE handyfuncs_mod
    USE atom_cont_patch_mod, only: atom_cont_patch_t
    USE angle_rad_patch_mod, only: rad_patch_t
    USE atom_rate_patch_mod, only: atom_rate_t
    USE atom_mod, only: atom_type, atomic_spectrum_type, use_classic_col, kramer_bf, debye_shield, &
                        red_cut, blue_cut, logf_cut, local_data_t, m1d_legacy_mode, nq_fac
    USE trans_mod, only: line_type, cont_type
    USE mconst, only: CC
    USE ion_mod, only: ion_type, particle_type
    USE element_mod, only: element_type
    USE molecule_mod, only: molec_type
    USE bbtrans_mod, only: bbtrans_type
    USE composition_mod, only: composition_type
    USE m3d_patch_mod, only: m3d_patch_t, maxiter, check_NaN, cswitch, key_link_t, &
    global_maxerr, convlim, atom_update_counter, comoving_frame, global_pos, &
    decouple_continuum, overlapping_bb, overlapping_bf, line_blanketing, converged
    USE bftrans_mod, only: bftrans_type
    USE hprofl_mod, only: hydrogen_BPO
    USE osmet_mod, only:osmet_type
    USE atom_line_patch_mod, only:atom_line_t, external_sp
    USE net_opac_patch_mod, only:net_opac_t
    USE m3d_debug_mod
    implicit none
  
    !-----------------------------------------------------------------------------
  
    type, public, extends(m3d_patch_t) :: atom_patch_t
      type(atom_type) :: output_atom, bb_atom, bf_atom
      class(atomic_spectrum_type), pointer :: output_sp, bf_sp, bb_sp
      real(8) :: maxerr=huge(1.0)
      real(8), dimension(:), allocatable :: lvlerr
      logical :: update_skipped=.false., atom_allocated=.false., molec
    contains
      procedure :: pre_init => atom_patch_pre_init
      procedure :: execute
      procedure :: nbor_relations => atom_patch_nbor_relations
      procedure :: create_sub_patches
      procedure :: pre_init_atom
      procedure :: pre_init_spectrum
      procedure :: init_atom_data
      procedure :: update_populations
      procedure :: check_skip
      procedure :: get_osmet_nu
      procedure :: exclude_double_continua
      procedure :: exclude_double_lines
      procedure :: print_post_update_info
      procedure :: check_in_bounds
      procedure :: save_abundance_to_file
      procedure :: load_mem
    end type atom_patch_t

    integer, save :: epos(3), elvl
    logical, save :: level_partf=.true.
    logical, save :: exclude_from_line_list=.true.
    logical, save :: detailed_absmet=.false.
    character(len=128), save :: precomputed_depart=''
    character(len=128), public :: atom_file=''

    !-----------------------------------------------------------------------------
  
  contains
    
  !===============================================================================
  SUBROUTINE atom_patch_pre_init (self)
    class(atom_patch_t), target:: self
    integer, save :: itimer=0
    integer :: n_nu, n_mu
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

    n_nu = self%n_nu
    n_mu = self%n_mu

    if (self%maxiter > 0) then
      if (decouple_continuum) then
        call self%create_sub_patches([1, n_nu],   [0, n_mu], self%bb_atom, self%bb_sp)
        call self%create_sub_patches([1, n_nu]+n_nu, [0, 0], self%bf_atom, self%bf_sp)
      else
        call self%create_sub_patches([1, n_nu], [0, n_mu], self%atom, self%sp)
      endif
    endif

    call self%create_sub_patches([-n_nu, -1], [-self%n_out, 0], self%output_atom, self%output_sp)

    call self%post_pre_init()
    self%atom%fcheck = self%fcheck
    self%atom%pcheck = self%pcheck
    self%atom%linecheck = self%linecheck
    call trace%end(itimer)
  END SUBROUTINE atom_patch_pre_init

  !===============================================================================
  SUBROUTINE create_sub_patches (self, nu, mu, atom, sp)
    class(atom_patch_t) :: self
    integer, intent(in) :: nu(2), mu(2)
    class(atom_type), intent(in), target :: atom
    class(atomic_spectrum_type), intent(in), pointer :: sp
    class(atom_rate_t), pointer :: rate
    class(atom_cont_patch_t), pointer :: cont
    class(atom_line_t), pointer :: line
    integer :: n, m
    logical :: cf, dc, new_cont, new_line, new_rate
    !-----------------------------------------------------------------------------
    if (.not. associated(sp)) return
    if (sp%nnu == 0) return
    do n=nu(1), nu(2)
        cf = comoving_frame .and. n>0
        dc = decouple_continuum
      do m=mu(1), mu(2)
        new_cont = m==0 .and. (merge(n>self%n_nu, .true., dc) .or. n<0) .and. associated(atom%cont)
        new_line = (m<0 .or. merge(m==0, m/=0, cf)) .and. associated(atom%line)
        new_rate = merge(m==0, m>0, cf) .or. (dc .and. m==0 .and. n>self%n_nu)

        if (new_cont) then
          allocate(cont)
          call cont%pre_init_and_link(self, i_mu=m, i_nu=n, sp=sp, atom=atom)
        endif
        if (new_line) then
          allocate(line)
          call line%pre_init_and_link(self, i_mu=m, i_nu=n, sp=sp, atom=atom)
        endif
        if (new_rate) then
          allocate(rate)
          call rate%pre_init_and_link(self, i_mu=m, i_nu=n, sp=sp, atom=atom)
        endif
      end do
    end do
  END SUBROUTINE create_sub_patches

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
    integer :: iostat
    namelist /atom_params/ atom_file, qnorm, convlim, exclude_trace_cont, level_partf, &
    use_atom_abnd, logf_cut, red_cut, blue_cut, n_in_lines, abundance, debye_shield, &
    exclude_from_line_list, precomputed_depart, use_classic_col, kramer_bf, overlapping_bf, &
    overlapping_bb, nq_fac, hydrogen_BPO, detailed_absmet
    namelist /line_mask/ in_lines
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
      write (io%namelist, atom_params)

      if (atom_file == '') then
        print *, 'WARNING: no atom specified'
        call trace%end(itimer)
        return
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

      write (*,*) ' - Reading atom'
      call atom%read(atom_file, qnorm, in_lines=in_lines, read_coll=.true.)
      atom%name = 'input'

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
    endif

    if (atom_file == '') then
      call trace%end(itimer)
      return
    endif

    allocate(self%atom)
    self%atom = atom
    self%atom%name = 'main'

    ! this pointer will be referenced by other atoms with same local data
    allocate(self%atom%data)

    if (use_atom_abnd .or. abundance > -99) then
      call self%composition%overwrite_abundance(atom%id, atom%abnd, atom%awgt)
    else
      self%atom%abnd = self%atom%abnd + self%composition%FeH
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
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%pre_init_spectrum', itimer=itimer)
    if (first_time) then
      first_time = .false.
      external_sp = output_sp%nnu > 0 
      if (associated(self%atom)) then
        if (decouple_continuum) then
          bf_sp%name = 'bound-free'
          call self%bf_atom%create_spectrum(bf_sp, bf=.true., sort=overlapping_bf)
          call self%bf_atom%set_ired_iblue(bf_sp, bf=.true., set_active=.true.)
          call bf_sp%save('bf')

          bb_sp%name = 'bound-bound'
          call self%bb_atom%create_spectrum(bb_sp, bb=.true., sort=overlapping_bb)
          call self%bb_atom%set_ired_iblue(bb_sp, bb=.true., set_active=.true.)
          call bb_sp%save('bb')

        else
          sp%name = 'all'
          if (detailed_absmet) call self%get_osmet_nu(sp, self%composition%osmet)
          call self%atom%create_spectrum(sp, bb=.true., bf=.true.)
          call self%atom%set_ired_iblue(sp, bb=.true., bf=.true., set_active=.true.)
          call sp%save('all')
        endif

        if (external_sp) then
          call self%output_atom%set_ired_iblue(output_sp, bb=.true., bf=.true.) 
        else
          output_sp%name = 'output based on atom'
          call self%output_atom%create_spectrum(output_sp, bb=.true.)
          call self%output_atom%set_ired_iblue(output_sp, bb=.true., bf=.true.) 
        endif
      endif
      call output_sp%save('out')
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
  SUBROUTINE init_atom_data (self)
    class(atom_patch_t), target :: self
    class(atom_type), pointer :: atom
    class(element_type), pointer :: element
    class(molec_type), pointer :: molec
    real(4), dimension(:,:,:), pointer, contiguous :: dxi, totn, df, saha
    real(4), dimension(:,:,:), pointer :: temp, ne
    real(4), dimension(:,:,:,:), pointer, contiguous :: n, nstar
    real(4), dimension(:,:,:,:), pointer :: nh
    integer, save:: itimer=0
    integer :: i
    !-----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%init_atom_data', itimer=itimer)
    atom => self%atom

    allocate(atom%C(atom%nlevel, atom%nlevel))
    allocate(atom%W(atom%nlevel, atom%nlevel))
    allocate(atom%e(atom%nlevel))
    allocate(atom%oldn(atom%nlevel))
    allocate(self%lvlerr(atom%nlevel))

    ne   => self%m3d_parent%mem(:,:,:,self%i_ne,            self%it,1)
    temp => self%m3d_parent%mem(:,:,:,self%i_temp,          self%it,1)
    nh   => self%m3d_parent%mem(:,:,:,self%i_nh1:self%i_nh6,self%it,1)

    atom%data%nstar => self%mem(:,:,:,:,self%it,1)
    atom%data%n     => self%mem(:,:,:,:,self%it,2)
    n     => atom%data%n    
    nstar => atom%data%nstar

    call os%allocate(atom%data%totn, self%n, label='totn')
    totn => atom%data%totn

    element => null()
    molec   => null()

    call self%composition%get_totn(atom%id, atom%abnd, totn, element, molec)
    if (.not. (associated(element) .or. associated(molec))) then
      call self%save_abundance_to_file()
    endif

    if (associated(molec)) atom%tmax_update = self%composition%tmolim

    if (level_partf) nullify(element)

    call os%allocate(saha, self%n)
    call os%allocate(dxi, self%n)
    saha = real(self%composition%saha)
    if (debye_shield) then
      dxi = real(self%composition%dxi)
    else
      dxi = 0
    endif

    if (precomputed_depart /= '') then
      call self%load_mem(precomputed_depart, self%out_mem, self%kind)
      do i=1, atom%nlevel
        where (totn > 0)
          n(:,:,:,i) = n(:,:,:,i) / nstar(:,:,:,i)
        end where
      end do
      call atom%compute_lte_populations(nstar, totn, temp, saha, dxi, element, molec)
      do i=1, atom%nlevel
        n(:,:,:,i) = n(:,:,:,i) * nstar(:,:,:,i) 
      end do
    else
      call atom%compute_lte_populations(nstar, totn, temp, saha, dxi, element, molec)
      do i=1, atom%nlevel
        n(:,:,:,i) = nstar(:,:,:,i) 
      end do
    endif

    call atom%compute_cont_factors(temp, dxi, saha)
    call os%deallocate(dxi)

    call trace%end(itimer)
  END SUBROUTINE init_atom_data

  !===============================================================================
  SUBROUTINE execute (self)
    class(atom_patch_t) :: self
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: nbtask
    class(line_type), pointer :: line
    real(8), dimension(:,:,:,:), pointer, contiguous :: gamma_ij, gamma_ji
    real(4), dimension(:,:,:,:), pointer, contiguous :: wphi
    integer :: i
    integer, save:: count=0
    real(8) :: global
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin (trim(self%kind)//'%update', itimer=itimer)
    if (self%iter == 0) then
      call self%init_atom_data()
    else
      wphi => self%atom%data%wphi
      gamma_ij => self%atom%data%gamma_ij
      gamma_ji => self%atom%data%gamma_ji

      if (self%iter == 1) then
        self%dtime = self%dt_fixed

        ! Normalize wphi
        do i=1, self%atom%nline
          wphi(:,:,:,i) = 1 / wphi(:,:,:,i)
        end do
      endif

      !if (.not. self%atom%stiff) then
        do i=1, self%atom%nline
          gamma_ij(:,:,:,i) = gamma_ij(:,:,:,i) * wphi(:,:,:,i)
          gamma_ji(:,:,:,i) = gamma_ji(:,:,:,i) * wphi(:,:,:,i)
        end do

        call self%update_populations(gamma_ij, gamma_ji)
        self%update_skipped = .false.
      !else
      !  self%update_skipped = .true.
      !endif
    endif
    
    if (self%iter > 0) then
      count = atom_update_counter%update(+1)
    endif
    if (count >= self%task_list%dims(3)) then
      global = global_maxerr%update(0.d0)
      if ((global > 0) .and. global < convlim .and. self%maxerr < convlim) then
        write(stdout, *) self%ipos(3)+1, 'Converged!', count
        if (self%iter < self%maxiter) converged = .true.
        call self%finish()
      endif
    endif

    if (associated(self%atom%data%gamma_ij)) then
      call os%deallocate(self%atom%data%gamma_ij, label='atom%gamma_ij')
      call os%deallocate(self%atom%data%gamma_ji, label='atom%gamma_ji')
    endif
    call trace%end(itimer)
  END SUBROUTINE execute

  !===============================================================================
  SUBROUTINE print_post_update_info (self)
    class(atom_patch_t):: self
    class(local_data_t), pointer :: adata
    character(len=64) :: text
    integer :: i, l, p(3), ipos
    real(8) :: global
    !-----------------------------------------------------------------------------
    if (self%iter == 0) then
      write(text,'(A,I2,A,I2,A)') 'LTE Iteration .'
    else
      write(text,'(A,I2,A,I2,A)') 'Iteration ', self%iter, '/', self%maxiter, '.'
    endif
    call io%header(text)
    if (self%update_skipped) then
      call self%print_to_stdout(text='update skipped')
    elseif (self%m3d_patch_t%check_in_bounds()) then
      i = self%lvlcheck
      l = self%linecheck
      p = self%pcheck
      adata => self%atom%data
      !arr(1) = self%lvlerr(i)
      memcheck(1) = adata%n(p(1),p(2),p(3),i)
      memcheck(2) = adata%totn(p(1),p(2),p(3))
      if (associated(adata%gamma_ij) .and. l<=self%atom%nrad) then
        memcheck(3) = adata%gamma_ij(p(1),p(2),p(3),l)
        memcheck(4) = adata%gamma_ji(p(1),p(2),p(3),l)
        if (l<=self%atom%nline) then
          memcheck(5) = adata%wphi(p(1),p(2),p(3),l)
        endif
      else
        memcheck(3) = 10.0**(self%atom%abnd-12.0)
        memcheck(4) = self%composition%toth(p(1),p(2),p(3))
      endif
      call self%print_to_stdout(arr=memcheck(1:5))
    else
      call self%print_to_stdout()
    endif
    global = global_maxerr%update(0.d0)
    ipos = global_pos%update(+0)
    write(stdout,'(A19, ES19.7, A10, ES15.7, A9, I2, 2X, 4I4)') 'local err =', self%maxerr, &
      ', global =', global, ', ipos = ', ipos+1, elvl, epos
    write(stdout,'(1X, A80)') &
    '................................................................................'
  END SUBROUTINE print_post_update_info

  !===============================================================================
  SUBROUTINE check_skip (self)
    class(atom_patch_t) :: self
    class(nbor_t), pointer :: nbor
    class(task_t), pointer :: nbtask
    integer :: iter, atom_iter
    real(8) :: time
    ! ----------------------------------------------------------------------------
    time=huge(1.d0)
    nbor => self%link%nbor
    do while (associated(nbor))
      nbtask => nbor%task
      select type(nbtask)
      class is (atom_rate_t)
        if (.not. nbtask%out_ray) then
          time = min(time, nbtask%rad%time)
        endif
      end select
      nbor => nbor%next
    end do
    iter = max(0, floor(time - self%time))
    atom_iter = max(0, min(self%maxiter - self%iter - 1, iter))
    call self%skip_iterations(atom_iter)
    nbor => self%link%nbor
    do while (associated(nbor))
      if (nbor%needs_me) then
        nbtask => nbor%task
        select type(nbtask)
        class is (m3d_patch_t)
          if (.not. nbtask%out_ray .and. .not. nbtask%is_finished()) then
            call nbtask%skip_iterations(iter)
          endif
        end select
      endif
      nbor => nbor%next
    end do
  END SUBROUTINE check_skip

  !===============================================================================
  SUBROUTINE update_populations (self, gamma_ij, gamma_ji)
    class(atom_patch_t), target :: self
    class(atom_type), pointer :: atom
    integer, dimension(3) :: ep
    real(4), pointer :: temp, ne, totn
    real(4), dimension(:), pointer :: nh
    real(8), dimension(:), pointer :: Gij, Gji
    real(4), dimension(:), pointer :: n, nstar
    real(8), dimension(:,:,:,:), target :: gamma_ij, gamma_ji
    integer :: i, j, k, ii, ipos, il
    integer, save :: itimer=0
    real(8) :: global
    ! ----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%update_populations', itimer=itimer)
    atom => self%atom
    self%lvlerr = 0.0d0
    self%maxerr = 0.0d0

    do k=1, self%n(3)
      do j=1, self%n(2)
        do i=1, self%n(1)
          ne    => self%m3d_parent%mem(i,j,k,self%i_ne,            self%it,1)
          temp  => self%m3d_parent%mem(i,j,k,self%i_temp,          self%it,1)
          nh    => self%m3d_parent%mem(i,j,k,self%i_nh1:self%i_nh6,self%it,1)
          totn  => atom%data%totn  (i,j,k)
          n     => atom%data%n     (i,j,k,:)
          nstar => atom%data%nstar (i,j,k,:)
          Gij   => gamma_ij   (i,j,k,:)
          Gji   => gamma_ji   (i,j,k,:)

          if (totn > 0) then
            call atom%solve_stat_equilibrium(temp, ne, totn, n, nstar, nh, Gij, Gji, cswitch, check_NaN)

            ! keep track of maximum population change
            do ii=1,atom%nlevel
              if (abs(atom%e(ii)) > self%lvlerr(ii)) then
                  self%lvlerr(ii) = abs(atom%e(ii))
                  if (self%lvlerr(ii) > self%maxerr) then
                    self%maxerr = self%lvlerr(ii)
                    ep = [i, j, k + self%zs]
                    il = ii
                  endif
              end if
            end do
          endif
        end do
      end do
    end do

    global = global_maxerr%update(0.d0)
    ipos = global_pos%update(+0)
    if (ipos == self%ipos(3) .or. self%maxerr > global) then
      if (ipos /= self%ipos(3)) call global_pos%reset(self%ipos(3))
      call global_maxerr%reset(self%maxerr)
      epos = ep
      elvl = il
    endif
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
    logical :: same
    logical, save :: first_time=.true.
    !-----------------------------------------------------------------------------
    do i=1, size(self%composition%trans)
      trans => self%composition%trans(i)
      if (.not. trans%active) cycle
      select case(trans%label)
      case ('(H+H)', 'H2+', 'H- FF', 'H FF', 'H2-', 'H- BF')
        ! do not exclude
      case default
        if (trans%level%g > 0) then
          if (uppercase(trans%ion%id(1:2)) == self%atom%id) then
            do kc=1, self%atom%ncont
              if (.not. self%atom%cont(kc)%active) cycle
              lower = self%atom%cont(kc)%i

              if (self%atom%ev(lower)==0) then
                same = trans%level%xion == self%atom%ev(lower)
              else
                same = abs(trans%level%xion / self%atom%ev(lower) - 1) < 1d-2
              endif

              if (same) then
                trans%active = .false.
                if (first_time) then
                  print *, 'excluding transition ', trans%label
                endif
                exit
              endif
            end do
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
    call self%composition%element_by_id(trim(self%atom%id), element)
    if (.not. associated(element)) then
      call self%composition%molecule_by_id(trim(self%atom%id), molec)
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
      class is (atom_line_t)
        needed = .false.
        if (associated(another%atom)) then
          if (another%atom%id /= self%atom%id) then
            needs_me = .false.
          endif
        else
          needs_me = .false.
        endif
      class is (atom_cont_patch_t)
        needed = .false.
        if (another%atom%id /= self%atom%id) then
          needs_me = .false.
        endif
      class is (atom_rate_t)
        if (another%atom%id /= self%atom%id) then
          needed = .false.; needs_me = .false.
        endif
      class is (net_opac_t)
        needed = .false.
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
    cib = .false.
    
    if (self%m3d_patch_t%check_in_bounds()) then
      cib = self%lvlcheck > 0 .and. self%lvlcheck <= self%atom%nlevel
    elseif (self%verbose > 2) then
      cib = .true.
    endif
  END FUNCTION check_in_bounds

  !===============================================================================
  SUBROUTINE save_abundance_to_file (self)
    class(atom_patch_t), target :: self
    integer :: lunit, i
    !-----------------------------------------------------------------------------
    if (self%ipos(3)==0) then
      open(newunit=lunit, file=trim(io%run_dir)//'abund.txt', form='FORMATTED', access='append')
      write(lunit, '(A3, F8.4)') self%atom%id, self%atom%abnd
      close(lunit)
    endif
  END SUBROUTINE save_abundance_to_file

  !===============================================================================
  SUBROUTINE load_mem (self, folder, mem, label)
    class(atom_patch_t), intent(inout) :: self
    character(len=128), intent(in) :: folder
    character(len=*), intent(in) :: label
    real(4), intent(inout), dimension(:,:,:,:,:) :: mem
    real(4), dimension(:,:,:,:,:), pointer, contiguous :: mem2
    character(len=128) :: filename, sfolder, metafile, file
    integer :: lunit, file_size, mem_size, iostat, prec, oldm(5), newm(5), iv, iw, k
    logical, save :: first_time=.true.
    ! ----------------------------------------------------------------------------
    call self%check_access(self)
    write(sfolder, '(A, A)') trim(folder), 'save/'
    write(filename, '(A, A, 3(I3.3), A)') trim(label), '_', self%key(1:3), '.bin'
    write(metafile, '(A, A, 3(I3.3), A)') trim(label), '_meta.txt'
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
    newm = shape(mem)
    open(newunit=lunit, file=trim(sfolder)//filename, access='stream')
    if (all(oldm == newm)) then
      read(lunit) mem
    elseif (all(oldm(3:5) == newm(3:5))) then
      if (first_time) then
        first_time = .false.
        write(stdout, '(A, 2(I6,"x",I4))') 'Interpolating patch memory when loading', oldm(1:2), newm(1:2)
      endif
      call os%allocate(mem2, oldm)
      read(lunit) mem2

      mem2(:,:,:,:,2) = mem2(:,:,:,:,2) / mem2(:,:,:,:,1)

      do iw=1, newm(5)
      do iv=1, newm(4)
      do  k=1, newm(3)
        call interp2d_periodic(mem2(:,:,k,iv,iw), mem(:,:,k,iv,iw))
      end do
      end do
      end do

      mem(:,:,:,:,2) = mem(:,:,:,:,2) * mem(:,:,:,:,1)

      call os%deallocate(mem2)
    else
      write(stdout, *) 'Cannot load patch save ', trim(sfolder)//trim(filename), &
      'Memory sizes do not match!'
      print *, oldm
      print *, newm
      stop
    endif
    close(lunit)
  END SUBROUTINE load_mem

  END MODULE
