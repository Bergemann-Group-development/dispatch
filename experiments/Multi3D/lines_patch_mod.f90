module lines_patch_mod
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
    USE atom_var, only: atom_type
    USE trans_mod, only: line_type
    USE m3d_patch_mod, only: m3d_patch_t, maxiter
    USE ray_patch_mod, only: ray_patch_t
    implicit none
  
    !-----------------------------------------------------------------------------
  
    type, public, extends(ray_patch_t) :: lines_patch_t
      class(atom_type), pointer :: atom
      type(line_type), dimension(:), contiguous, pointer :: line
      real(4), dimension(:,:,:), contiguous, pointer :: doppler_factor
      integer, dimension(:), contiguous, pointer :: idx_map
      integer :: nline, ired, iblue
    contains
      procedure :: pre_init => lines_patch_pre_init
      procedure :: init => lines_patch_init
      procedure :: update => lines_patch_update
      !procedure :: post_update => linelist_post_update
      !procedure :: finish => linelist_finish
      procedure :: nbor_relations => lines_patch_nbor_relations
      procedure :: create_line_indices
      procedure :: create_line_frequencies
      procedure :: interpolate_lines
    end type lines_patch_t

    !-----------------------------------------------------------------------------
  
  contains

  !===============================================================================
  SUBROUTINE lines_patch_pre_init (self)
    class(lines_patch_t), target :: self
    class(line_type), pointer :: line
    integer :: i, nv
    !---------------------------------------------------------------------------
    call self%ray_patch_t%pre_init()
    call self%append_to_task_list()
    call self%create_line_indices()
    call self%create_line_frequencies()

    nv = 0
    do i=1, self%nline
      line => self%line(i)
      nv = nv + line%nq
    end do

    call self%set_nv_and_nw(nv, 1)
    self%ng = 0
  END SUBROUTINE lines_patch_pre_init

  !===============================================================================
  SUBROUTINE lines_patch_init (self)
    class(lines_patch_t) :: self
    class(line_type), pointer :: line
    integer :: i, fs, fe
    !---------------------------------------------------------------------------
    call self%ray_patch_t%init()

    fs = 1
    do i=1, self%nline
      line => self%line(i)
      fe = fs + line%nq - 1
      line%phi => self%mem(:,:,:,fs:fe,1,1)
      fs = fe + 1
    end do
  END SUBROUTINE lines_patch_init

  !===============================================================================
  SUBROUTINE create_line_indices (self)
    class(lines_patch_t) :: self
    class(line_type), pointer :: line
    integer :: i, nline, fs, fe, idx
    logical :: skip
    !---------------------------------------------------------------------------
    self%ired  = self%fe
    self%iblue = self%fs

    call os%allocate(self%idx_map, [self%atom%nline], zero=.true.)
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
    do i=1, self%atom%nline
      idx = self%idx_map(i)
      if (idx==0) cycle

      self%line(idx) = self%atom%line(i)
    end do
  END SUBROUTINE create_line_indices

  !===============================================================================
  SUBROUTINE create_line_frequencies (self)
    class(lines_patch_t) :: self
    class(line_type), pointer :: line
    integer :: i, f, fs, fe
    real(8) :: nus, nue
    real(8), allocatable :: nu(:)
    !---------------------------------------------------------------------------
    nus = self%sp%nu(self%fs)
    nue = self%sp%nu(self%fe)

    do i=1, self%nline
      line => self%line(i)
      if (allocated(line%nu)) deallocate(line%nu)
      if (allocated(line%q))  deallocate(line%q)
      line%nnu = line%nq ! reset number of frequencies to input
      allocate(line%nu(line%nnu))
      call line%freql(self%atom%qnorm, symmetric=.true.)

      fs = 1
      fe = line%nq
      do f=1, line%nq
        if (line%nu(f) < nus) then
          fs = f
        elseif (line%nu(f) > nue) then
          fe = f
          exit
        endif
      end do
      
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
    end do
  END SUBROUTINE create_line_frequencies

  !===============================================================================
  SUBROUTINE lines_patch_update (self)
    use spectrum_mod, only: set_weights
    class(lines_patch_t) :: self
    class(line_type), pointer :: line
    real(4), dimension(:,:,:), pointer, contiguous :: tothi, ne, temp, adamp, doppler_factor, vel_dop, ni, nj
    real(4), dimension(:,:,:,:), pointer, contiguous :: phi
    real(4), dimension(:,:,:,:), pointer :: nh
    real(8), dimension(:), allocatable :: wnu
    integer :: i, f
    integer, save:: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('lines_patch_t%lines_patch_update', itimer=itimer)
    call os%allocate(tothi, self%n)
    call os%allocate(adamp, self%n)

    ne   => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_ne,            self%it,1)
    temp => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_temp,          self%it,1)
    nh   => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_nh1:self%i_nh6,self%it,1)
    doppler_factor => self%atom%doppler_factor
    tothi = sum(nh(:,:,:,1:5), 4) ! total neutral hydrogen density

    if (self%iter == 0) then
      call self%comp_vproj(self%vproj)
      self%vproj = (1 + self%vproj) / doppler_factor
    endif

    do i=1, self%nline
      line => self%line(i)
      adamp = line%damping(temp, ne, nh(:,:,:,6), tothi, doppler_factor, self%atom%awgt)
      !print *, i, 176, line%estimate_qmax(temp(1,1,176), ne(1,1,176), nh(1,1,176,6), tothi(1,1,176), self%atom%awgt)
      !print *, temp(1,1,176), ne(1,1,176), nh(1,1,176,6), tothi(1,1,176)
      line%adamp => adamp
      call line%profile(line%nu, self%vproj, doppler_factor, line%phi)

      if (self%iter == 0) then
        call os%allocate(line%wphi, self%n, zero=.true.)
        call set_weights(wnu, line%nu, 1, line%nnu)

        do f=1,line%nnu
          line%wphi = line%wphi + line%phi(:,:,:,f) * wnu(f)
          line%phi(:,:,:,f) = line%phi(:,:,:,f) * h_4pi * line%nu(f)
        end do
        deallocate(wnu)
        call os%allocate(line%chi_fac, self%n)
        call os%allocate(line%eta_fac, self%n)
      endif

      ni => self%atom%n(:,:,:,line%i)
      nj => self%atom%n(:,:,:,line%j)

      line%eta_fac = line%Aji * nj
      line%chi_fac = line%Bij * ni - line%Bji * nj
    enddo

    call os%deallocate(tothi)
    call os%deallocate(adamp)
    call self%post_update()
    call trace%end(itimer)
  END SUBROUTINE lines_patch_update

  !===============================================================================
  SUBROUTINE interpolate_lines (self, chi, eta)
    class(lines_patch_t) :: self
    class(line_type), pointer :: line
    real(4), dimension(:,:,:,:), intent(inout) :: chi, eta
    real(8), dimension(:), contiguous, pointer :: hnu_4pi, nu, q
    real(4), dimension(:,:,:,:), contiguous, pointer :: scratch_phi, phi
    integer :: i, fcs, fce, fs, fe, f, newq, ff, i1, i2, m(4), nnu, p(3)
    logical :: skip
    !---------------------------------------------------------------------------
    fcs = max(self%ired,  self%fs)
    fce = min(self%iblue, self%fe)

    nu => self%sp%nu(fcs:fce)

    m = shape(chi)
    m(4) = size(nu)
    call os%allocate(scratch_phi, m)

    do i=1, self%nline
      line => self%line(i)

      call line%check_range(fcs,fce,fs,fe,skip)
      if (skip) print *, 'line not active. Something is wrong!'

      nu => self%sp%nu(fs:fe)
      nnu = size(nu)
      phi => scratch_phi(:,:,:,1:nnu)
  
      call line%interpolate_phi(nu, phi)

      !if (self%fcheck >= fs .and. self%fcheck <= fe) then
      !  f = self%fcheck-fs+1
      !  ff = self%fcheck-self%fs+1
      !  p = self%pcheck
      !  print '(8(ES15.7))', chi(p(1),p(2),p(3),ff) + phi(p(1),p(2),p(3),f) * scr2(p(1),p(2),p(3)), &
      !    phi(p(1),p(2),p(3),f) * scr2(p(1),p(2),p(3)), phi(p(1),p(2),p(3),f), ni(p(1),p(2),p(3)), &
      !    nj(p(1),p(2),p(3)), nu(f)
      !endif
  
      do f=1,nnu
        ff = fs - self%fs + f
        eta(:,:,:,ff) = eta(:,:,:,ff) + phi(:,:,:,f) * line%eta_fac
        chi(:,:,:,ff) = chi(:,:,:,ff) + phi(:,:,:,f) * line%chi_fac
      end do
    end do

    call os%deallocate(scratch_phi)
  END SUBROUTINE interpolate_lines

  !===============================================================================
  SUBROUTINE lines_patch_nbor_relations(self, another, needed, needs_me, download)
    class(lines_patch_t) :: self
    class(task_t), pointer :: another        ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                            ! does it need me, should I download?
    class(task_t), pointer :: other
    logical :: default
    ! ----------------------------------------------------------------------------
    call trace%begin ('lines_patch_t%nbor_relations')
    call self%ray_patch_t%nbor_relations(another, needed, needs_me, download)

    if (needed .or. needs_me) then
      select type(another)
      class is (m3d_patch_t)
        if (another%id == self%m3d_parent%id) then
          needs_me = .false.
        elseif (another%kind == 'net_opac') then
          needed = .false.
        elseif(another%kind == 'rate_patch') then
          needed = .false.
        else
          needed = .false.; needs_me = .false.
        endif
      class default
        needed = .false.; needs_me = .false.
      end select
    endif
    call trace%end()
  END SUBROUTINE lines_patch_nbor_relations

  end module
