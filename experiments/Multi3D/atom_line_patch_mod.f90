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
    USE atom_mod, only: atom_type
    USE trans_mod, only: line_type
    USE m3d_patch_mod, only: m3d_patch_t, maxiter, hash_string, comoving_frame
    USE ray_patch_mod, only: ray_patch_t
    USE op_math_mod, only: get_interval
    implicit none
  
    !-----------------------------------------------------------------------------
  
    type, public, extends(ray_patch_t) :: atom_line_t
      class(atom_type), pointer :: atom
      type(line_type), dimension(:), contiguous, pointer :: line
      real(4), dimension(:,:,:), contiguous, pointer :: doppler_factor
      integer, dimension(:), contiguous, pointer :: idx_map
      integer :: nline, ired, iblue
    contains
      procedure :: pre_init => atom_line_pre_init
      procedure :: init => atom_line_init
      procedure :: update => atom_line_update
      procedure :: deallocate_mem => atom_line_deallocate_mem
      procedure :: nbor_relations => atom_line_nbor_relations
      procedure :: create_line_indices
      procedure :: create_line_frequencies
      procedure :: interpolate_lines
      procedure :: print_post_update_info
      procedure :: check_in_bounds
      procedure :: integrate_wphi
      procedure :: init_nbors_hash
    end type atom_line_t

    logical, public, save :: second_sp=.false.

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

    nv = 0
    do i=1, self%nline
      line => self%line(i)
      nv = nv + line%nq
    end do

    call self%set_nv_and_nw(nv, 1)
    call self%m3d_patch_t%init()

    fs = 1
    do i=1, self%nline
      line => self%line(i)
      fe = fs + line%nq - 1
      line%phi => self%mem(:,:,:,fs:fe,1,1)
      fs = fe + 1
    end do
    self%print_update_info = self%check_in_bounds()
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
    real(8), allocatable :: nu(:)
    !---------------------------------------------------------------------------
    nus = self%sp%nu(self%fs)
    nue = self%sp%nu(self%fe)

    pre_pad_nus = self%sp%nu(self%pre_pad_fs)
    pre_pad_nue = self%sp%nu(self%pre_pad_fe)

    do i=1, self%nline
      line => self%line(i)
      if (second_sp .and. self%out_ray) then
        line%nq = line%nnu
        !call line%freql(self%atom%qnorm, symmetric=.true.)
      else
        if (allocated(line%nu)) deallocate(line%nu)
        if (allocated(line%q))  deallocate(line%q)
        line%nnu = line%nq ! reset number of frequencies to input
        call line%freql(self%atom%qnorm, symmetric=.true.)

        call get_interval(line%nu, nus, nue, fs, fe, limit=.true.)

        ! this is needed for computing wphi
        call get_interval(line%nu, pre_pad_nus, pre_pad_nue, line%fs, line%fe, limit=.true.)

        line%fs = line%fs - fs + 1
        line%fe = line%fe - fs + 1
        
        line%nnu = fe-fs+1
        allocate(nu(line%nnu))
        nu = line%nu(fs:fe)
        deallocate(line%nu)
        allocate(line%nu(line%nnu))
        line%nu = nu
        deallocate(nu)

        line%nq = line%nnu
        allocate(line%q(line%nq))
        line%q = line%nu/line%nu0 - 1
      endif
    end do
  END SUBROUTINE create_line_frequencies

  !===============================================================================
  SUBROUTINE atom_line_update (self)
    class(atom_line_t) :: self
    class(line_type), pointer :: line
    real(4), dimension(:,:,:), pointer, contiguous :: tothi, ne, temp, adamp, doppler_factor
    real(4), dimension(:,:,:), pointer, contiguous :: inv_df, ni, nj, unsold_temp, stark_pe, vproj
    real(4), dimension(:,:,:,:), pointer, contiguous :: phi
    real(4), dimension(:,:,:,:), pointer :: nh
    real(8), dimension(:), allocatable :: wnu
    integer :: i, f, p(3)
    integer, save:: itimer=0
    !-----------------------------------------------------------------------------
    if (self%is_finished()) return
    call trace%begin('atom_line_t%atom_line_update', itimer=itimer)

    call self%check_access(self%m3d_parent%m3d_parent)
    ne   => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_ne,            self%it,1)
    temp => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_temp,          self%it,1)
    nh   => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_nh1:self%i_nh6,self%it,1)
    doppler_factor => self%atom%doppler_factor

    call os%allocate(tothi,       self%n)
    call os%allocate(adamp,       self%n)
    call os%allocate(unsold_temp, self%n)
    call os%allocate(stark_pe,    self%n)
    call os%allocate(vproj,       self%n)
    call os%allocate(inv_df,      self%n)
    
    call self%comp_vproj(vproj)
    inv_df = 1 / doppler_factor

    tothi = sum(nh(:,:,:,1:5), 4) ! total neutral hydrogen density
    unsold_temp =  8.0*BK / PI * (1.0/MH + 1.0/self%atom%awgt) * temp
    stark_pe = ne * bk * temp**one_sixth

    do i=1, self%nline
      line => self%line(i)
      adamp = 0
      call line%damping(adamp, temp, ne, tothi, unsold_temp, stark_pe)
      if (self%atom%id == 'H') then
        call line%dampingH(adamp, temp, ne, nh(:,:,:,6))
      endif
      
      call line%profile(line%q, vproj, inv_df, adamp, line%phi)

      if (.not. self%out_ray) then
        call os%allocate(line%wphi, self%n, zero=.true., label=trim(self%kind)//'%wphi')
        call self%sp%set_weights(wnu, line%q(line%fs:line%fe), line%fs, line%fe)
        do f=line%fs, line%fe
          line%wphi = line%wphi + line%phi(:,:,:,f) * wnu(f)
        end do
        deallocate(wnu)
      endif

      do f=1,line%nnu
        line%phi(:,:,:,f) = line%phi(:,:,:,f) * over_4pi * (line%q(f) + 1)
      end do
    enddo

    call os%deallocate(tothi)
    call os%deallocate(adamp)
    call os%deallocate(unsold_temp)
    call os%deallocate(stark_pe)
    call os%deallocate(vproj)
    call os%deallocate(inv_df)

    if (.not. self%out_ray) then
      call self%integrate_wphi()
    endif
    call self%post_update()
    call trace%end(itimer)
  END SUBROUTINE atom_line_update

  !===============================================================================
  SUBROUTINE integrate_wphi (self)
    class(atom_line_t) :: self
    class(line_type), pointer :: line
    integer :: i, k
    !-----------------------------------------------------------------------------
    call self%m3d_parent%mem_lock(1)%set('integrate_wphi')
    do i=1, self%nline
      k = self%line(i)%ntrans
      line => self%atom%line(k)
      if (.not. associated(line%wphi)) then
        call os%allocate(line%wphi, self%n, zero=.true., label=trim(self%m3d_parent%kind)//'%wphi')
      endif
      line%wphi = line%wphi + self%wmu * self%line(i)%wphi
      call os%deallocate(self%line(i)%wphi, label=trim(self%kind)//'%wphi')
    end do
    call self%m3d_parent%mem_lock(1)%unset('integrate_wphi')
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
    !---------------------------------------------------------------------------
    fcs = max(self%ired,  self%fs)
    fce = min(self%iblue, self%fe)

    nu => self%sp%nu(fcs:fce)

    m = shape(chi)
    m(4) = size(nu)
    call os%allocate(scratch_phi, m, label='scratch_phi')
    call os%allocate(chi_fac, self%n)
    call os%allocate(eta_fac, self%n)

    !if (.not. associated(self%mem)) then
    !  print *, 'mem not allocated. Something is wrong!'
    !endif

    do i=1, self%nline
      line => self%line(i)
      call line%check_range(fcs,fce,fs,fe,skip)
      !if (skip) print *, 'line not active. Something is wrong!'

      nu => self%sp%nu(fs:fe)
      nnu = size(nu)
      phi => scratch_phi(:,:,:,1:nnu)
  
      call line%interpolate_phi(nu, phi)

      !if (self%fcheck >= fs .and. self%fcheck <= fe) then
      !  f = self%fcheck-fs+1
      !  ff = self%fcheck-self%fs+1
      !  p = self%pcheck
      !  print '(8(ES15.7))', chi(p(1),p(2),p(3),ff) + phi(p(1),p(2),p(3),f) * line%chi_fac(p(1),p(2),p(3)), &
      !    phi(p(1),p(2),p(3),f) * line%chi_fac(p(1),p(2),p(3)), phi(p(1),p(2),p(3),f), &
      !    self%atom%n(p(1),p(2),p(3),line%i), self%atom%n(p(1),p(2),p(3),line%j)
      !endif

      if (LTE) then
        ni => self%atom%nstar(:,:,:,line%i)
        nj => self%atom%nstar(:,:,:,line%j)
      else
        ni => self%atom%n(:,:,:,line%i)
        nj => self%atom%n(:,:,:,line%j)
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
  END SUBROUTINE interpolate_lines

  !===============================================================================
  SUBROUTINE atom_line_nbor_relations(self, another, needed, needs_me, download)
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
      class is (m3d_patch_t)
        if (another%id == self%m3d_parent%id) then
          needs_me = .false.
        elseif (another%kind == 'net_opac') then
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
  END SUBROUTINE atom_line_nbor_relations

  !===============================================================================
  SUBROUTINE atom_line_deallocate_mem(self, deallocator, print)
    class(atom_line_t) :: self
    class(m3d_patch_t), intent(in), target :: deallocator
    logical, optional :: print
    ! ----------------------------------------------------------------------------
    if (self%is_finished() .and. self%is_clear(bits%busy)) then
      call self%m3d_patch_t%deallocate_mem(deallocator, print)
      if (associated(self%line)) deallocate(self%line)
      if (associated(self%idx_map)) call os%deallocate(self%idx_map)
    endif    
  END SUBROUTINE atom_line_deallocate_mem

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
    call self%print_to_stdout(arr=phi(p(1),p(2),p(3),:))
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
    call self%m3d_patch_t%init_nbors_hash(nbors)

    key(1:3) = self%key(1:3)
    key(4) = hash_string('net_opac')
    call self%check_nbor(key, nbors)

    key(4) = hash_string('atom_rate')
    call self%check_nbor(key, nbors)
    call trace%end (itimer)
  END SUBROUTINE init_nbors_hash

  END MODULE
