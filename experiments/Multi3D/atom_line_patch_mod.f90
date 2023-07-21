MODULE atom_line_patch_mod
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
    !USE math_mod
    USE bits_mod
    !USE units_mod
    !USE scratch_mod
    USE mconst
    USE handyfuncs_mod
    USE spectrum_mod, only: set_weights
    USE ion_mod, only: ion_type
    USE atom_mod, only: atom_type
    USE trans_mod, only: line_type
    USE m3d_halo_mod, only:lock_mem, unlock_mem
    USE m3d_patch_mod, only: m3d_patch_t, maxiter, comoving_frame
    USE ray_patch_mod, only: ray_patch_t
    USE op_math_mod, only: searchsorted
    USE hprofl_mod, only: hydrogen_profile, prepare_hydrogen_damping
    USE m3d_debug_mod
    implicit none
  
    !-----------------------------------------------------------------------------
  
    type, public, extends(ray_patch_t) :: atom_line_t
      type(ray_patch_t), pointer :: net_opac
      type(line_type), dimension(:), contiguous, pointer :: line
      integer, dimension(:), contiguous, pointer :: idx_map
      integer :: nline, ired, iblue
    contains
      procedure :: pre_init => atom_line_pre_init
      procedure :: init => atom_line_init
      procedure :: execute
      procedure :: allocate_mem
      procedure :: deallocate_mem
      procedure :: nbor_relations
      procedure :: create_line_indices
      procedure :: create_line_frequencies
      procedure :: interpolate_lines
      procedure :: print_post_update_info
      procedure :: check_in_bounds
      procedure :: integrate_wphi
      procedure :: init_nbors_hash
      procedure :: add_to_net_opac
    end type atom_line_t

    logical, public, save :: external_sp=.false.

    !-----------------------------------------------------------------------------
  
  contains

  !===============================================================================
  SUBROUTINE atom_line_pre_init (self)
    class(atom_line_t), target :: self
    class(line_type), pointer :: line
    integer :: i
    integer, save :: itimer=0
    !---------------------------------------------------------------------------
    call trace%begin('atom_line_t%pre_init', itimer=itimer)
    self%kind = 'atom_line'
    call self%ray_patch_t%pre_init()
    self%ng = 0
    call self%post_pre_init()
    call trace%end(itimer)
  END SUBROUTINE atom_line_pre_init

  !===============================================================================
  SUBROUTINE atom_line_init (self)
    class(atom_line_t) :: self
    class(line_type), pointer :: line
    integer :: i, fs, fe, nv
    !---------------------------------------------------------------------------
    if (comoving_frame .and. .not. self%out_ray) then
      call self%frequency_padding()
    endif
    call self%create_line_indices()
    call self%create_line_frequencies()

    nv = sum(self%line%nq)
    call self%set_nv_and_nw(nv, 1)
    call self%m3d_patch_t%init()

    fs = 1
    do i=1, self%nline
      line => self%line(i)
      fe = fs + line%nq - 1
      line%phi => self%mem(:,:,:,fs:fe,1,1)
      fs = fe + 1
    end do
  END SUBROUTINE atom_line_init

  !===============================================================================
  SUBROUTINE create_line_indices (self)
    class(atom_line_t) :: self
    class(line_type), pointer :: line
    integer :: i, nline, fs, fe, idx
    logical :: skip
    !---------------------------------------------------------------------------
    self%ired  = self%fe
    self%iblue = self%fs

    call os%allocate(self%idx_map, self%atom%nline, zero=.true.)
    nline = 0
    do i=1, self%atom%nline
      line => self%atom%line(i)
      if (.not. line%active) cycle
      call line%check_range(self%fs, self%fe, fs, fe, skip)
      if (skip) cycle
      nline = nline + 1 
      self%idx_map(i) = nline
      self%ired  = min(self%ired,  fs)
      self%iblue = max(self%iblue, fe)
    end do

    self%nline = nline
    allocate(self%line(self%nline))
    call os%mem_register(storage_size(self%line), self%nline, label='atom%line')

    do i=1, self%atom%nline
      idx = self%idx_map(i)
      if (idx==0) cycle
      self%line(idx) = self%atom%line(i)
    end do
  END SUBROUTINE create_line_indices

  !===============================================================================
  SUBROUTINE create_line_frequencies (self)
    class(atom_line_t) :: self
    class(line_type), pointer :: line
    integer :: i, f, fs, fe
    real(8) :: nus, nue, pre_pad_nus, pre_pad_nue
    real(8), dimension(:), allocatable :: nu, wq
    !---------------------------------------------------------------------------
    if (self%sp%sorted) then
      nus = self%sp%nu(self%fs)
      nue = self%sp%nu(self%fe)
      pre_pad_nus = self%sp%nu(self%pre_pad_fs)
      pre_pad_nue = self%sp%nu(self%pre_pad_fe)
    endif

    do i=1, self%nline
      line => self%line(i)
      if (external_sp .and. self%out_ray) then
        line%nq = line%nnu
      else
        if (.not. self%sp%sorted) then
          nus = merge(self%sp%nu(self%fs), line%nu_input(1), self%fs>line%ired)
          nue = merge(self%sp%nu(self%fe), line%nu_input(line%nq), self%fe<line%iblue)
          pre_pad_nus = nus
          pre_pad_nue = nue
        endif
        fs = max(searchsorted(line%nu_input, nus, outbounds_ok=.true.), 1)
        fe = min(searchsorted(line%nu_input, nue, outbounds_ok=.true.)+1, size(line%nu_input))

        line%nu => line%nu_input(fs:fe)
        line%wq => line%wq(fs:fe)
        line%nnu = fe-fs+1

        call os%allocate(line%q , line%nnu, 'line%q')
        line%q = line%nu / line%nu0 - 1

        ! this is needed for computing wphi
        line%fs = max(searchsorted(line%nu, pre_pad_nus, outbounds_ok=.true.), 1)
        line%fe = min(searchsorted(line%nu, pre_pad_nue, outbounds_ok=.true.), line%nnu)
      endif
    end do
  END SUBROUTINE create_line_frequencies

  !===============================================================================
  SUBROUTINE execute (self)
    class(atom_line_t) :: self
    class(line_type), pointer :: line
    class(ion_type), pointer :: He1
    real(4), dimension(:,:,:), pointer, contiguous :: tothi, ne, temp, adamp
    real(4), dimension(:,:,:), pointer, contiguous :: vproj
    real(4), dimension(:,:,:), pointer, contiguous :: inv_df, ni, nj, vbar2, stark_pe
    real(4), dimension(:,:,:,:), pointer, contiguous :: phi
    real(4), dimension(:,:,:,:), pointer :: nh, chi, eta
    integer :: i, f, p(3), fs, fe, iw, it
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin (trim(self%kind)//'%update', itimer=itimer)
    p = self%pcheck
    call self%check_access(self%m3d_parent%m3d_parent)
    ne   => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_ne,            self%it,1)
    temp => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_temp,          self%it,1)
    nh   => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_nh1:self%i_nh6,self%it,1)
    call self%composition%ion_by_id('He I', He1)

    call os%allocate(tothi,       self%n)
    call os%allocate(adamp,       self%n)
    call os%allocate(vbar2, self%n)
    call os%allocate(stark_pe,    self%n)
    call os%allocate(vproj,       self%n)
    call os%allocate(inv_df,      self%n)
    
    call self%comp_vproj(vproj)

    tothi = sum(nh(:,:,:,1:5), 4) ! total neutral hydrogen density
    vbar2 = self%atom%compute_vbar2(temp)
    inv_df = self%atom%inverse_doppler_factor(temp)
    stark_pe = ne * temp**one_sixth

    do i=1, self%nline
      line => self%line(i)
      select case (line%profile_type)
      case ('voigt')
        adamp = 0
        call line%damping(adamp, ne, tothi, vbar2, stark_pe)
        if (self%atom%id == 'H') then
          call line%dampingH(adamp, temp, ne, nh(:,:,:,6))
        endif
        call line%profile(line%q, vproj, inv_df, adamp, line%phi)

      case ('hlinop')
        call hydrogen_profile(line, line%q, temp, ne, nh(:,:,:,1), vproj, inv_df, real(He1%n), line%phi)
      end select

      if (self%print_update_info) then
        if (i==self%linecheck) then
          memcheck(2) = adamp(p(1),p(2),p(3))
          memcheck(3) = vproj(p(1),p(2),p(3))
        endif
      endif

      if (.not. self%out_ray .and. self%iter==0) then
        call os%allocate(line%wphi, self%n, zero=.true., label=trim(self%kind)//'%wphi')
        do f=line%fs, line%fe
          line%wphi = line%wphi + line%phi(:,:,:,f) * line%wq(f)
        end do
        if (self%print_update_info) then
          if (i==self%linecheck) then
            memcheck(4) = line%wphi(p(1),p(2),p(3))
          endif
        endif
      endif

      do f=1,line%nnu
        line%phi(:,:,:,f) = line%phi(:,:,:,f) * (line%q(f) + 1) * over_4pi
      end do
    enddo

    call os%deallocate(tothi)
    call os%deallocate(adamp)
    call os%deallocate(vbar2)
    call os%deallocate(stark_pe)
    call os%deallocate(inv_df)
    call os%deallocate(vproj)

    !call self%add_to_net_opac()

    if (.not. self%out_ray .and. self%iter==0) then
      call self%integrate_wphi()
    endif
    call trace%end(itimer)
  END SUBROUTINE execute

  !===============================================================================
  SUBROUTINE add_to_net_opac (self)
    class(atom_line_t) :: self
    real(4), dimension(:,:,:,:), pointer :: nh, chi, eta
    integer :: i, f, p(3), fs, fe, iw, it
    logical :: compressed
    !-----------------------------------------------------------------------------
    call lock_mem(self%net_opac, 'add_to_net_opac')
    call self%net_opac%check_mem(locked=.true.)
    do iw=1, self%net_opac%nw
      do it=1, self%net_opac%nt
        do f=1, self%net_opac%nv
          self%net_opac%mem(:,:,:,f,it,iw) = 0
        end do
      end do
    end do
    self%net_opac%nullified = .true.

    it = self%net_opac%it
    chi => self%net_opac%authority_mem(:,:,:,:,it,1)
    eta => self%net_opac%authority_mem(:,:,:,:,it,2)
    call self%interpolate_lines(chi, eta, LTE=.false.)
    if (self%out_ray) then
      chi => self%net_opac%authority_mem(:,:,:,:,it,3)
      eta => self%net_opac%authority_mem(:,:,:,:,it,4)
      call self%interpolate_lines(chi, eta, LTE=.true.)
    endif
    ! do not unlock the task here. Net_opac is locked until its update
  END SUBROUTINE add_to_net_opac

  !===============================================================================
  SUBROUTINE integrate_wphi (self)
    class(atom_line_t) :: self
    class(line_type), pointer :: line
    class(atom_type), pointer :: atom
    real(4), dimension(:,:,:), pointer, contiguous :: wphi
    integer :: i, k
    !-----------------------------------------------------------------------------
    call lock_mem(self%m3d_parent, 'integrate_wphi')
    atom => self%atom
    if (.not. associated(atom%data%wphi)) then
      call os%allocate(atom%data%wphi, [self%n, atom%nline], zero=.true., label=trim(self%m3d_parent%kind)//'%wphi')
    endif
    do i=1, self%nline
      k = self%line(i)%ntrans
      wphi => atom%data%wphi(:,:,:,k)
      wphi = wphi + self%wmu * self%line(i)%wphi
      call os%deallocate(self%line(i)%wphi, label=trim(self%kind)//'%wphi')
    end do
    call unlock_mem(self%m3d_parent, 'integrate_wphi')
  END SUBROUTINE integrate_wphi

  !===============================================================================
  SUBROUTINE interpolate_lines (self, chi, eta, LTE)
    class(atom_line_t) :: self
    class(line_type), pointer :: line
    real(4), dimension(:,:,:,:), intent(inout) :: chi, eta
    logical, intent(in) :: LTE
    real(8), dimension(:), contiguous, pointer :: hnu_4pi, nu, q
    real(4), dimension(:,:,:), pointer, contiguous :: ni, nj, chi_fac, eta_fac
    real(4), dimension(:,:,:,:), contiguous, pointer :: scratch_phi, phi
    integer :: i, fcs, fce, fs, fe, f, newq, ff, i1, i2, m(4), nnu, p(3), x,y,z
    logical :: skip
    integer, save :: itimer=0
    !---------------------------------------------------------------------------
    call trace%begin('atom_line_t%interpolate_lines', itimer=itimer)
    fcs = max(self%ired,  self%fs)
    fce = min(self%iblue, self%fe)

    nu => self%sp%nu(fcs:fce)

    m = shape(chi)
    m(4) = size(nu)
    call os%allocate(scratch_phi, m, label='scratch_phi')
    call os%allocate(chi_fac, self%n)
    call os%allocate(eta_fac, self%n)

    do i=1, self%nline
      line => self%line(i)
      call line%check_range(fcs,fce,fs,fe,skip)
      !if (skip) print *, 'line not active. Something is wrong!'

      nu => self%sp%nu(fs:fe)
      nnu = size(nu)
      phi => scratch_phi(:,:,:,1:nnu)
  
      call line%interpolate_phi(nu, phi)

      if (LTE) then
        ni => self%atom%data%nstar(:,:,:,line%i)
        nj => self%atom%data%nstar(:,:,:,line%j)
      else
        ni => self%atom%data%n(:,:,:,line%i)
        nj => self%atom%data%n(:,:,:,line%j)
      endif

      eta_fac = line%Aji * nj
      chi_fac = line%Bij * ni - line%Bji * nj
  
      do f=1,nnu
        ff = fs - self%fs + f
        do z=1,m(3)
        do y=1,m(2)
        do x=1,m(1)
          eta(x,y,z,ff) = eta(x,y,z,ff) + phi(x,y,z,f) * eta_fac(x,y,z)
          chi(x,y,z,ff) = chi(x,y,z,ff) + phi(x,y,z,f) * chi_fac(x,y,z)
        end do
        end do
        end do
      end do
    end do

    call os%deallocate(scratch_phi, label='scratch_phi')
    call os%deallocate(chi_fac)
    call os%deallocate(eta_fac)
    call trace%end(itimer)
  END SUBROUTINE interpolate_lines

  !===============================================================================
  SUBROUTINE nbor_relations(self, another, needed, needs_me, download)
    class(atom_line_t) :: self
    class(task_t), pointer :: another        ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                            ! does it need me, should I download?
    class(task_t), pointer :: other
    logical :: default
    ! ----------------------------------------------------------------------------
    call trace%begin ('atom_line_t%nbor_relations')
    call self%ray_patch_t%nbor_relations(another, needed, needs_me, download)

    if (needed .or. needs_me) then
      select type(another)
      class is (atom_line_t)
        ! do nothing
      class is (m3d_patch_t)
        if (another%id == self%m3d_parent%id) then
          needs_me = .false.
        elseif (another%kind == 'net_opac') then
          select type (another)
          class is (ray_patch_t)
            self%net_opac => another
          end select
          needed = .false.
        elseif(another%kind == 'atom_rate') then
          needed = .false.
        else
          needed = .false.; needs_me = .false.
        endif
      class default
        needed = .false.; needs_me = .false.
      end select
    endif
    call trace%end()
  END SUBROUTINE nbor_relations

  !===============================================================================
  SUBROUTINE allocate_mem(self, allocator, m_in, label)
    class(atom_line_t) :: self
    class(m3d_patch_t), intent(in), target :: allocator
    class(line_type), pointer :: line
    integer, intent(in), optional :: m_in(6)
    character(len=*), optional :: label
    integer :: i, fs, fe
    ! ----------------------------------------------------------------------------
    call self%m3d_patch_t%allocate_mem(allocator, m_in)
    fs = 1
    do i=1, self%nline
      line => self%line(i)
      fe = fs + line%nq - 1
      line%phi => self%mem(:,:,:,fs:fe,1,1)
      fs = fe + 1
    end do
  END SUBROUTINE allocate_mem

  !===============================================================================
  SUBROUTINE deallocate_mem(self, deallocator, print)
    class(atom_line_t) :: self
    class(m3d_patch_t), intent(in), target :: deallocator
    logical, optional :: print
    ! ----------------------------------------------------------------------------
    call self%m3d_patch_t%deallocate_mem(deallocator, print)
    if (self%is_finished() .and. self%is_clear(bits%busy)) then
      if (associated(self%line)) deallocate(self%line)
      if (associated(self%idx_map)) call os%deallocate(self%idx_map)
    endif    
  END SUBROUTINE deallocate_mem

  !===============================================================================
  SUBROUTINE print_post_update_info (self)
    class(atom_line_t) :: self
    real(8), dimension(:), contiguous, pointer :: nu
    real(4), dimension(:,:,:,:), contiguous, pointer :: phi
    integer :: p(3), i
    !-----------------------------------------------------------------------------
    i = self%linecheck
    p = self%pcheck
    nu => self%sp%nu(self%fcheck:self%fcheck)
    call os%allocate(phi, [self%n, 1])
    call self%line(i)%interpolate_phi(nu, phi)
    memcheck(1) = phi(p(1),p(2),p(3),1)
    if (self%out_ray) then
      call self%print_to_stdout(arr=memcheck(1:3))
    else
      call self%print_to_stdout(arr=memcheck(1:4))
    endif
    call os%deallocate(phi)
  END SUBROUTINE print_post_update_info
  
  !===============================================================================
  FUNCTION check_in_bounds (self) result (cib)
    class(atom_line_t) :: self
    integer :: i, fs, fe
    logical :: skip, cib
    !-----------------------------------------------------------------------------+
    cib = .false.
    if (.not. associated(self%idx_map)) return
    if (self%m3d_patch_t%check_in_bounds()) then
      if (self%linecheck <= self%atom%nline) then
        if (self%fcheck >= self%fs .and. self%fcheck <= self%fe) then
          i = self%idx_map(self%linecheck)
          if (i > 0) then
            call self%line(i)%check_range(self%fcheck,self%fcheck,fs,fe,skip)
            if (.not. skip) then
              cib=.true.
            endif
          endif
          self%linecheck = i
        endif
      endif
    endif
  END FUNCTION check_in_bounds

  !===============================================================================
  SUBROUTINE init_nbors_hash (self, nbors)
    class(atom_line_t):: self
    class(nbor_t), pointer:: nbors
    !.............................................................................
    class(task_t), pointer:: task
    class(nbor_t), pointer:: nbor
    integer, save:: itimer=0
    integer :: i, key(4)
    !-----------------------------------------------------------------------------
    call trace%begin('atom_line_t%init_nbors_hash', itimer=itimer)
    self%interest%l = self%authority%l - self%ds
    self%interest%u = self%authority%u + self%ds
    call self%m3d_patch_t%init_nbors_hash(nbors)
    call self%ray_patch_t%init_nbors_hash(nbors)

    key(1:3) = self%key(1:3)
    key(4) = hash_string('net_opac')
    call self%check_nbor(key, nbors)

    key(4) = hash_string('atom_rate')
    call self%check_nbor(key, nbors)
    call trace%end (itimer)
  END SUBROUTINE init_nbors_hash

  END MODULE
