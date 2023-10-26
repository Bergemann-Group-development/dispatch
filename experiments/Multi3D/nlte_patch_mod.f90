MODULE nlte_patch_mod
  !USE rt_integral_mod
  !USE radau_mod
  !USE boundaries_mod
  !USE eos_mod
  !USE list_mod
  USE os_mod
  USE halo_mod
  USE gpatch_mod
  USE link_mod
  USE task_mod
  !USE index_mod
  !USE mesh_mod
  USE io_mod
  !USE aux_mod
  USE timer_mod
  USE scaling_mod
  USE trace_mod
  !USE mpi_mod
  !USE omp_mod
  !USE omp_timer_mod
  !USE math_mod
  USE bits_mod
  USE os_mod
  !USE units_mod
  !USE scratch_mod
  USE patch_mod, only: patch_t
  USE mmath, only: cc_kernel
  USE op_math_mod, only: pchip_mono
  USE handyfuncs_mod
  USE atom_mod, only:atom_type
  USE spectrum_mod, only:set_weights
  USE ray_patch_mod, only:ray_patch_t
  USE m3d_halo_mod, only: lock_mem, unlock_mem
  USE m3d_patch_mod, only: all_vec=>vec, all_wmu=>wmu, use_vel, m3d_patch_t, &
  absmet_in_output, pre_update, long_wmu, long_vec, ray_offset, pivot_z, dlam
  USE ray_mod, only: one_dimensional, ray_type
  USE trans_mod, only:line_type, cont_type
  USE bbtrans_mod, only:voigt_function2
  USE m3d_debug_mod
  implicit none
  private

  !-----------------------------------------------------------------------------
  type, public, extends(ray_patch_t) :: nlte_patch_t
    real(4), dimension(:,:,:,:), allocatable :: Jnu, Psi, PsiJ
    real(4), dimension(:), allocatable :: twohnu3_c2, h4pi_nu
    real(4), dimension(:,:,:,:,:,:), pointer :: mem_remap
    real(8), dimension(:), pointer :: nu
    integer :: mus, mue, n_set, nf, n_mu
  contains
    procedure :: pre_init => nlte_patch_pre_init
    procedure :: execute
    procedure :: lte_background_opacities
    procedure :: add_bound_free_opacities
    procedure :: bound_bound_damping
    procedure :: add_bound_bound_opacities
    procedure :: bound_bound_rates
    procedure :: bound_free_rates
    procedure :: back_to_work
    procedure :: eta_to_Snu
    procedure :: set_upstream_boundary
    procedure :: integrate_intensity
    procedure :: integrate_rates
    procedure :: work_to_mem
    procedure :: save_intensity
    procedure :: rotate_cube
    procedure :: short_char
    procedure :: long_char
    procedure :: get_vproj
    procedure :: print_post_update_info
    procedure :: nbor_relations => nlte_patch_nbor_relations
    procedure :: init_nbors_hash
    procedure :: allocate_mem
  end type nlte_patch_t

  !-----------------------------------------------------------------------------
  integer, save:: verbose=0                           ! verbosity
  integer :: f, x, y, z, nx, ny, nz, fs, fe, nf
  real(4), dimension(:,:,:,:,:), allocatable :: work, back
  !$omp threadprivate(work, back, x, y, z, nx, ny, nz, f, fs, fe, nf)
    
contains

!===============================================================================
SUBROUTINE nlte_patch_pre_init (self)
  class(nlte_patch_t), target:: self
  integer, save :: itimer=0
  ! ----------------------------------------------------------------------------
  call trace%begin('nlte_patch_t%_pre_init', itimer=itimer)
  self%kind = 'nlte_patch'
  call self%m3d_patch_t%pre_init()

  self%dir = merge(1, -1, abs(self%i_mu)==1)
  self%n_mu = merge(size(long_wmu), size(all_wmu)/2, self%out_ray)
  self%mus = merge(1, self%n_mu+1, self%dir>0)
  self%mue = merge(self%n_mu, size(all_wmu), self%dir>0)

  call self%set_freq_indices()
  self%nf = self%fe-self%fs+1
  self%n_set = merge(3, 1, self%out_ray)

  self%ng(3) = 1
  call self%post_pre_init()

  call trace%end(itimer)
END SUBROUTINE nlte_patch_pre_init

!===============================================================================
SUBROUTINE execute (self)
  class(nlte_patch_t) :: self
  real(4), dimension(:,:,:,:), allocatable :: gamma_ij, gamma_ji, adamp
  real(4), dimension(:,:,:), pointer, contiguous :: temp
  real(4), dimension(:,:,:), allocatable :: inv_df, vproj
  real(4) :: wmu
  real(8) :: vec(3)
  logical :: get_rates, use_linelist, use_atom, bound_bound, bound_free
  integer :: i_mu, i, zs, ze, p(3), ng, i_set, l
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('nlte_patch%execute', itimer=itimer)
  ng = self%ng(3)
  fs = self%fs
  fe = self%fe
  nx = self%patch%n(1)
  ny = self%patch%n(2)
  nz = self%patch%n(3)
  zs = 1 - ng
  ze = nz + ng
  p = self%pcheck

  get_rates = .not. self%out_ray
  use_linelist = self%out_ray .and. self%composition%ll_initialised

  use_atom = associated(self%atom)
  if (use_atom) then
    bound_bound = associated(self%atom%line)
    bound_free = associated(self%atom%cont)
  else
    bound_bound = .false.
    bound_free = .false.
  endif

  self%mem_remap(1:3,fs:fe,1:nx,1:ny,self%mus:self%mue,1:self%n_set) => self%mem
  self%nu(fs:fe) => self%sp%nu(fs:fe)
  temp => self%m3d_parent%mem(:,:,:,self%i_temp,self%it,1)

  call allocate_arrays()

  self%twohnu3_c2 = twoh_c2 * self%nu**3
  self%h4pi_nu = 4*PI / (hh*self%nu)

  inv_df = self%atom%inverse_doppler_factor(temp)

  call self%lte_background_opacities(temp, include_absmet=get_rates.or.absmet_in_output)

  if (use_atom) then
    if (bound_free) call self%add_bound_free_opacities()
    if (bound_bound) call self%bound_bound_damping(adamp, temp, inv_df)
  endif

  do i_set=1, self%n_set
    ! Starting loop over angles
    do i_mu=self%mus, self%mue
      if (get_rates) then
        wmu = all_wmu(i_mu)
        vec = all_vec(:,i_mu)
      else
        wmu = long_wmu(i_mu)
        vec = long_vec(:,i_mu)
      endif
      wmu = wmu * over_4pi
      call self%get_vproj(vproj, vec)
      call self%back_to_work()

      if (i_set < 3) then ! i_set==3 is the continuum only
        if (bound_bound) call self%add_bound_bound_opacities(wmu, inv_df, vproj, adamp, LTE=i_set==2)
        if (use_linelist) then
          call self%composition%add_linelist_opacities(self%nu, work(1,:,:,:,1:nz), work(2,:,:,:,1:nz), back(3,:,:,:,:), vproj, dlam)
        endif
      endif

      call self%eta_to_Snu()
      call self%set_upstream_boundary(i_mu, i_set, ng)

      if (one_dimensional .or. self%out_ray) then
        if (.not. one_dimensional) call self%rotate_cube(vec)
        call self%long_char(vec, get_rates)
      else
        call self%short_char(vec, get_rates)
      endif

      if (self%print_update_info .and. i_mu==self%mus) then
        memcheck(1:3) = work(:, self%fcheck, p(1), p(2), p(3))
        memcheck(2:3) = memcheck(2:3) * self%twohnu3_c2(self%fcheck)
      endif

      call self%work_to_mem(i_mu, i_set)
      
      if (get_rates) then
        if (bound_free) call self%integrate_intensity(wmu)
        if (bound_bound) call self%bound_bound_rates(wmu, inv_df, vproj, adamp, gamma_ij, gamma_ji)
      endif
    end do
  end do

  if (get_rates) then
    if (bound_free) call self%bound_free_rates(gamma_ij, gamma_ji)
    if (use_atom) call self%integrate_rates(gamma_ij, gamma_ji)
  elseif (self%patch%boundaries%is_set(bits%zu)) then
    call self%save_intensity()
  endif
  call deallocate_arrays
  call trace%end(itimer)

contains

  !-----------------------------------------------------------------------------
  subroutine allocate_arrays()
    allocate(back(3,fs:fe,nx,ny,nz), source=0.)
    allocate(work(3,fs:fe,nx,ny,zs:ze), source=0.)
    allocate(self%twohnu3_c2(fs:fe))
    allocate(self%h4pi_nu(fs:fe))
    allocate(inv_df, vproj, mold=temp)
    call mem_register(back,      label='back')
    call mem_register(work,      label='work')
    if (bound_bound) then
      allocate(adamp(nx,ny,nz,self%atom%nline), source=0.)
      call mem_register(adamp, label='adamp')
    endif
    if (get_rates) then
      allocate(self%Psi(fs:fe,nx,ny,nz), source=0.)
      allocate(gamma_ij(nx, ny, nz, self%atom%nrad), source=0.)
      allocate(gamma_ji(nx, ny, nz, self%atom%nrad), source=0.)
      call mem_register(self%Psi, label='Psi')
      call mem_register(gamma_ij, label='gamma_ij')
      call mem_register(gamma_ji, label='gamma_ji')
      if (bound_free) then
        allocate(self%Jnu(fs:fe,nx,ny,nz), source=0.)
        allocate(self%PsiJ(fs:fe,nx,ny,nz), source=0.)
        call mem_register(self%Jnu,  label='Jnu')
        call mem_register(self%PsiJ, label='PsiJ')
      endif
    endif  
  end subroutine allocate_arrays

  !-----------------------------------------------------------------------------
  subroutine deallocate_arrays()
    call mem_unregister(work     , label='work')
    call mem_unregister(back     , label='back')
    deallocate(work, back)
    deallocate(self%twohnu3_c2, self%h4pi_nu)
    if (bound_bound) then
      call mem_unregister(adamp, label='adamp')
      deallocate(adamp)
    endif
    if (get_rates) then
      call mem_unregister(self%Psi, label='Psi')
      call mem_unregister(gamma_ij, label='gamma_ij')
      call mem_unregister(gamma_ji, label='gamma_ji')
      deallocate(gamma_ij, gamma_ji, self%Psi)
      if (bound_free) then
        call mem_unregister(self%Jnu , label='Jnu')
        call mem_unregister(self%PsiJ, label='PsiJ')
        deallocate(self%Jnu, self%PsiJ)
      endif
    endif
  end subroutine deallocate_arrays
END SUBROUTINE execute

!===============================================================================
SUBROUTINE lte_background_opacities (self, temp, include_absmet)
  class(nlte_patch_t) :: self
  real(4), contiguous, dimension(:,:,:) :: temp
  real(4) :: h_kT, hv_kT, bnu
  logical :: include_absmet
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('nlte_patch%lte_background_opacities', itimer=itimer)
  do z=1, nz
  do y=1, ny
  do x=1, nx
    ! computing boltzmann factor
    h_kT = h_k / temp(x,y,z)
    do f=self%fs, self%fe
      hv_kT = self%nu(f) * h_kT
      back(3,f,x,y,z) = exp(-hv_kT)
    end do
  end do
  end do
  end do

  ! getting background LTE absorption
  call self%composition%absorption(self%nu, back(1,:,:,:,:), blz=back(3,:,:,:,:), &
  include_absmet=include_absmet, sorted=self%sp%sorted)

  do z=1, nz
  do y=1, ny
  do x=1, nx
    do f=self%fs, self%fe
      ! The source function is the planck function in LTE, therefore eta=chi*Bnu
      associate(blz=>back(3,f,x,y,z), chi=>back(1,f,x,y,z), eta=>back(2,f,x,y,z))
        if (blz>1e-10) then
          bnu = blz / (1.0 - blz)
        else
          bnu = blz
        endif
        eta = chi * bnu
      end associate
    end do
  end do
  end do
  end do
  call trace%end(itimer)
END SUBROUTINE lte_background_opacities

!===============================================================================
SUBROUTINE add_bound_free_opacities (self)
  class(nlte_patch_t) :: self
  class(cont_type), pointer :: cont
  real(4), contiguous, dimension(:,:,:,:), pointer :: n
  real(4) :: ni, nj, tmp, alpha
  logical :: skip
  integer :: i, fsmin, femax
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('nlte_patch%add_bound_free_opacities', itimer=itimer)
  ! Summing contributions from all bound-free transitions
  n => self%atom%data%n
  fsmin = huge(1); femax = 0
  do i=1, self%atom%ncont
    cont => self%atom%cont(i)
    if (.not. cont%active) cycle
    call cont%check_range(self%fs, self%fe, fs, fe, skip)
    if (skip) cycle
    fsmin = min(fsmin, fs)
    femax = max(femax, fe)

    do z=1, nz
    do y=1, ny
    do x=1, nx
      ni = n(x,y,z,cont%i)
      nj = n(x,y,z,cont%j)
      tmp = self%atom%data%ni_over_nj_lte(x,y,z,i) * nj
      do f=fs, fe
        alpha = cont%alpha(f)
        work(1,f,x,y,z) = work(1,f,x,y,z) + alpha * ni   ! chi
        work(2,f,x,y,z) = work(2,f,x,y,z) + alpha * tmp  ! eta
      end do
    end do
    end do
    end do
  end do

  ! Converting the bf-contributions to chi and eta
  do z=1, nz
  do y=1, ny
  do x=1, nx
    do f=fsmin, femax
      associate(chi=>work(1,f,x,y,z), eta=>work(2,f,x,y,z), blz=>back(3,f,x,y,z))
        eta = eta * blz
        chi = chi - eta
      end associate
    end do

    do f=self%fs, self%fe
      ! adding background (LTE) and bound-free to give continuum opacities
      back(1,f,x,y,z) = back(1,f,x,y,z) + work(1,f,x,y,z)
      back(2,f,x,y,z) = back(2,f,x,y,z) + work(2,f,x,y,z)
    end do
  end do
  end do
  end do
  call trace%end(itimer)
END SUBROUTINE add_bound_free_opacities


!===============================================================================
SUBROUTINE bound_bound_damping (self, adamp, temp, inv_df)
  class(nlte_patch_t) :: self
  class(line_type), pointer :: line
  real(4), contiguous, dimension(:,:,:,:) :: adamp
  real(4), contiguous, dimension(:,:,:) :: temp, inv_df
  real(4), dimension(:,:,:), allocatable :: tothi, vbar2
  real(4), dimension(:,:,:), allocatable :: stark_pe, unsold, vproj
  real(4), dimension(:,:,:), pointer, contiguous :: ne
  real(4), dimension(:,:,:,:), pointer :: nh
  integer :: i
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin (trim(self%kind)//'%bound_bound_damping', itimer=itimer)
  allocate(tothi, vbar2, stark_pe, unsold, vproj, mold=temp)
  ne => self%m3d_parent%mem(:,:,:,self%i_ne,            self%it,1)
  nh => self%m3d_parent%mem(:,:,:,self%i_nh1:self%i_nh6,self%it,1)
  tothi = sum(nh(:,:,:,1:5), 4) ! total neutral hydrogen density
  vbar2 = self%atom%compute_vbar2(temp)
  stark_pe = ne * temp**one_sixth
  unsold = vbar2**0.3 * tothi

  ! computing damping coefficients for all bb-transitions
  do i=1, self%atom%nline
    line => self%atom%line(i)
    if (.not. line%active) cycle
    call line%damping(adamp(:,:,:,i), ne, tothi, vbar2, stark_pe, unsold)
    adamp(:,:,:,i) = adamp(:,:,:,i) * inv_df * over_4pi / line%nu0
  end do
  call trace%end(itimer)
END SUBROUTINE bound_bound_damping

!===============================================================================
SUBROUTINE add_bound_bound_opacities (self, wmu, inv_df, vproj, adamp, LTE)
  class(nlte_patch_t) :: self
  class(line_type), pointer :: line
  real(4), contiguous, dimension(:,:,:) :: inv_df, vproj
  real(4), contiguous, dimension(:,:,:,:) :: adamp
  real(4), contiguous, dimension(:,:,:,:), pointer :: n
  real(4) :: wmu, idf, vpr, a, idf_4pi, ni, nj, Bji_nj, Bij_ni, q, phi, widf, chi_fac
  integer :: i
  logical :: skip, LTE
  integer, save :: itimer=0, itimer1=0
  !-----------------------------------------------------------------------------
  call trace%begin (trim(self%kind)//'%add_bound_bound_opacities', itimer=itimer)
  ! Copying continuum opacities to angle dependent opacities
  if (LTE) then
    n => self%atom%data%nstar
  else
    n => self%atom%data%n
  endif

  do i=1, self%atom%nline
    line => self%atom%line(i)
    if (.not. line%active) cycle
    call line%check_range(self%fs, self%fe, fs, fe, skip)
    if (skip) cycle

    call trace%begin (trim(self%kind)//'%wphi', itimer=itimer1)
    
    ! Saving integrated line profiles wphi for profile normalisation
    if (self%iter == 0 .and. .not. self%out_ray) then
      call lock_mem(self%m3d_parent, 'wphi')
      do z=1, nz
      do y=1, ny
      do x=1, nx
        idf = inv_df(x,y,z)
        vpr = vproj(x,y,z)
        a = adamp(x,y,z,i)
        widf = wmu * idf
  
        do f=fs, fe
          associate(wphi=>self%atom%data%wphi(x,y,z,i))
            q = line%q(f)
            phi = (q - vpr) * idf
            phi = voigt_function2(a, phi) * widf
            wphi = wphi + phi * line%wq(f)
          end associate
        end do 
      end do
      end do
      end do
      call unlock_mem(self%m3d_parent, 'wphi')
    endif

    call trace%end(itimer1)

    ! Adding bound-bound opacities
    do z=1, nz
    do y=1, ny
    do x=1, nx
      ni = n(x,y,z,line%i)
      nj = n(x,y,z,line%j)
      idf = inv_df(x,y,z)
      vpr = vproj(x,y,z)
      a = adamp(x,y,z,i)
      idf_4pi = idf * over_4pi

      Bji_nj = line%Bji * nj
      Bij_ni = line%Bij * ni
      chi_fac = Bij_ni - Bji_nj

      do f=fs, fe
        associate(chi=>work(1,f,x,y,z), eta=>work(2,f,x,y,z))
          q = line%q(f)
          phi = (q - vpr) * idf
          phi = voigt_function2(a, phi)
          phi = phi * idf_4pi * (q + 1)
          chi = chi + phi * chi_fac
          eta = eta + phi * Bji_nj
        end associate
      end do 
    end do
    end do
    end do
  end do

  call trace%end(itimer)
END SUBROUTINE add_bound_bound_opacities

!===============================================================================
SUBROUTINE bound_bound_rates (self, wmu, inv_df, vproj, adamp, gamma_ij, gamma_ji)
  class(nlte_patch_t) :: self
  class(line_type), pointer :: line
  real(4), contiguous, dimension(:,:,:) :: inv_df, vproj
  real(4), contiguous, dimension(:,:,:,:) :: adamp, gamma_ij, gamma_ji
  real(4), contiguous, dimension(:,:,:,:), pointer :: n
  real(4), dimension(:), allocatable :: w_nu
  real(4) :: wmu, idf, vpr, a, idf_4pi, ni, nj, chi_fac, Bji_nj, q, phi
  real(4) :: Bij_ni, g_ij, g_ji, vvv, ufac
  integer :: i, zup, l
  logical :: skip
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin (trim(self%kind)//'%bound_bound_rates', itimer=itimer)
  allocate(w_nu(self%fs:self%fe))
  ! Computing bound-bound radiative rates
  n => self%atom%data%n
  do i=1, self%atom%nline
    line => self%atom%line(i)
    if (.not. line%active) cycle
    call line%check_range(self%fs, self%fe, fs, fe, skip)
    if (skip) cycle

    w_nu(fs:fe) = wmu * self%h4pi_nu(fs:fe) * line%wnu(fs:fe)

    do z=1, nz
      zup = z - self%dir
    do y=1,ny
    do x=1,nx
      ni = n(x,y,z,line%i)
      nj = n(x,y,z,line%j)
      idf = inv_df(x,y,z)
      vpr = vproj(x,y,z)
      a = adamp(x,y,z,i)
      idf_4pi = idf * over_4pi

      Bji_nj = line%Bji * nj
      Bij_ni = line%Bij * ni

      g_ij = 0
      g_ji = 0
      do f=fs, fe
        associate(Inu=>work(3,f,x,y,z), Psi=>self%Psi(f,x,y,z))
          q = line%q(f)
          phi = (q - vpr) * idf
          phi = voigt_function2(a, phi)
          phi = phi * idf_4pi * (q + 1)
          vvv = phi * w_nu(f)
          phi = phi * Psi
          g_ji = g_ji + vvv*(1. + Inu - phi*Bij_ni)
          g_ij = g_ij + vvv*(Inu - phi*Bji_nj)
        end associate
      end do 
      gamma_ji(x,y,z,i) = gamma_ji(x,y,z,i) + g_ji * line%Aji
      gamma_ij(x,y,z,i) = gamma_ij(x,y,z,i) + g_ij * line%Aij
    end do
    end do
    end do
  end do
  call trace%end(itimer)
END SUBROUTINE bound_bound_rates

!===============================================================================
SUBROUTINE bound_free_rates (self, gamma_ij, gamma_ji)
  class(nlte_patch_t) :: self
  class(cont_type), pointer :: cont
  real(4), contiguous, dimension(:,:,:,:) :: gamma_ij, gamma_ji
  real(4), contiguous, dimension(:,:,:,:), pointer :: n
  real(4), dimension(:), allocatable :: alpha_2hnu3_c2, w_alpha, w_nu
  real(4) :: ni, nj, tmp, eta_fac, q, phi
  real(4) :: g_ij, g_ji, vij, vji, ufac
  integer :: i, j
  logical :: skip
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin (trim(self%kind)//'%bound_free_rates', itimer=itimer)
  allocate(alpha_2hnu3_c2(self%fs:self%fe))
  allocate(w_alpha(self%fs:self%fe))
  allocate(w_nu(self%fs:self%fe))
  ! Computing bound-free radiative rates based on half integrated mean intensity
  n => self%atom%data%n
  do j=1, self%atom%ncont
    i = j + self%atom%nline
    cont => self%atom%cont(j)
    if (.not. cont%active) cycle
    call cont%check_range(self%fs,self%fe,fs,fe,skip)
    if (skip) cycle

    w_nu(fs:fe) = 0.5 * self%h4pi_nu(fs:fe) * cont%wnu(fs:fe)
    alpha_2hnu3_c2(fs:fe) = self%twohnu3_c2(fs:fe) * cont%alpha(fs:fe)
    w_alpha(fs:fe) = alpha_2hnu3_c2(fs:fe) * w_nu(fs:fe)
    
    do z=1,nz
    do y=1,ny
    do x=1,nx
      tmp = self%atom%data%ni_over_nj_lte(x,y,z,j)
      ni = n(x,y,z,cont%i)
      nj = n(x,y,z,cont%j) * tmp
      g_ij = 0
      g_ji = 0
      do f=fs, fe
        associate(Jnu=>self%Jnu(f,x,y,z), blz=>back(3,f,x,y,z), Psi=>self%PsiJ(f,x,y,z))
          vij = w_alpha(f)
          vji = vij * blz
          phi = Psi * cont%alpha(f)
          g_ji = g_ji + vji*(1. + Jnu - phi*ni)
          g_ij = g_ij + vij*(Jnu - phi*blz*nj)
        end associate
      end do
      gamma_ji(x,y,z,i) = gamma_ji(x,y,z,i) + g_ji * tmp
      gamma_ij(x,y,z,i) = gamma_ij(x,y,z,i) + g_ij
    end do
    end do
    end do
  end do
  call trace%end(itimer)
END SUBROUTINE bound_free_rates

!===============================================================================
SUBROUTINE set_upstream_boundary (self, i_mu, i_set, ng)
  class(nlte_patch_t) :: self
  integer :: i_mu, i_set, ng, i
  real(4) :: xxold(3), xxnew(ng)
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('nlte_patch%set_upstream_boundary', itimer=itimer)
  z = merge(lbound(work,5), ubound(work,5), self%dir > 0)
  nx = self%patch%n(1)
  ny = self%patch%n(2)
  nz = self%patch%n(3)

  if (self%boundaries%is_set(bits%zu) .and. self%dir < 0) then
    work(:,:,:,:,z) = 0.
  elseif (self%boundaries%is_set(bits%zl) .and. self%dir > 0) then
    xxold = [(i*self%ds(3), i=1, 3)]
    xxnew = [(i*self%ds(3), i=1-ng, 0)]
    do y=1,  ny
    do x=1,  nx
    do f=self%fs, self%fe
      do i=1, 2
        call pchip_mono(xxold, work(i,f,x,y,z+1:z+3), &
                        xxnew, work(i,f,x,y,z:z), sorted=.true.)
        work(i,f,x,y,z) = max(work(i,f,x,y,z), 0.)
      end do
      work(3,f,x,y,z) = work(2,f,x,y,z)
    end do
    end do
    end do
  else
    do f=self%fs, self%fe
      work(1,f,:,:,z) = self%mem_remap(1,f,:,:,i_mu,i_set)
      work(2,f,:,:,z) = self%mem_remap(2,f,:,:,i_mu,i_set)
      work(3,f,:,:,z) = self%mem_remap(3,f,:,:,i_mu,i_set)
    end do
  endif
  call trace%end(itimer)
END SUBROUTINE set_upstream_boundary

!===============================================================================
SUBROUTINE integrate_intensity (self, wmu)
  class(nlte_patch_t) :: self
  real(4) :: wmu, wmu2
  integer :: x, y, z, f
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('nlte_patch%integrate_intensity', itimer=itimer)
  wmu2 = 2*wmu
  do z=1, self%patch%n(3)
  do y=1, self%patch%n(2)
  do x=1, self%patch%n(1)
  do f=self%fs, self%fe
    self%Jnu(f,x,y,z) = self%Jnu(f,x,y,z) + wmu2 * work(3,f,x,y,z)
  end do
  end do
  end do
  end do
  do z=1, self%patch%n(3)
  do y=1, self%patch%n(2)
  do x=1, self%patch%n(1)
  do f=self%fs, self%fe
    self%PsiJ(f,x,y,z) = self%PsiJ(f,x,y,z) + wmu2 * self%Psi(f,x,y,z)
  end do
  end do
  end do
  end do
  call trace%end(itimer)
END SUBROUTINE integrate_intensity

!===============================================================================
SUBROUTINE integrate_rates (self, gamma_ij, gamma_ji)
  class(nlte_patch_t) :: self
  class(atom_type), pointer :: atom
  integer :: m(4), i, k
  real(4), dimension(:,:,:,:), contiguous :: gamma_ij, gamma_ji
  !-----------------------------------------------------------------------------
  call lock_mem(self%m3d_parent, 'integrate_rates')
  atom => self%atom
  if (.not. associated(atom%data%gamma_ij)) then
    m = [self%patch%n, atom%nrad]
    call os%allocate(atom%data%gamma_ij, m, zero=.true., label='atom%gamma_ij')
    call os%allocate(atom%data%gamma_ji, m, zero=.true., label='atom%gamma_ji')
  endif
  do i=1, size(gamma_ij, 4)
    atom%data%gamma_ij(:,:,:,i) = atom%data%gamma_ij(:,:,:,i) + gamma_ij(:,:,:,i)
    atom%data%gamma_ji(:,:,:,i) = atom%data%gamma_ji(:,:,:,i) + gamma_ji(:,:,:,i)
  end do
  call unlock_mem(self%m3d_parent, 'integrate_rates')
END SUBROUTINE integrate_rates

!===============================================================================
SUBROUTINE work_to_mem (self, i_mu, i_set)
  class(nlte_patch_t) :: self
  integer :: f, i_mu, i_set, z
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('nlte_patch%work_to_mem', itimer=itimer)
  z = merge(self%patch%n(3), 1, self%dir > 0)
  do f=self%fs, self%fe
    self%mem_remap(1,f,:,:,i_mu,i_set) = work(1,f,:,:,z)
    self%mem_remap(2,f,:,:,i_mu,i_set) = work(2,f,:,:,z)
    self%mem_remap(3,f,:,:,i_mu,i_set) = work(3,f,:,:,z)
  end do
  call trace%end(itimer)
END SUBROUTINE work_to_mem

!===============================================================================
SUBROUTINE back_to_work (self)
  class(nlte_patch_t) :: self
  !-----------------------------------------------------------------------------
  ! copy back opacities into angle dependent opacities
  do z=1, nz
  do y=1, ny
  do x=1, nx
  do f=self%fs, self%fe
    associate(chi=>work(1,f,x,y,z), eta=>work(2,f,x,y,z))
      chi = back(1,f,x,y,z)
      eta = back(2,f,x,y,z)
    end associate
  end do
  end do
  end do
  end do
END SUBROUTINE back_to_work

!===============================================================================
SUBROUTINE eta_to_Snu (self)
  class(nlte_patch_t) :: self
  !-----------------------------------------------------------------------------
  ! Converting eta into source function Snu=eta/chi
  do z=1, nz
  do y=1, ny
  do x=1, nx
    do f=self%fs, self%fe
      work(2,f,x,y,z) = work(2,f,x,y,z) / work(1,f,x,y,z)
    end do
  end do
  end do
  end do
END SUBROUTINE eta_to_Snu

!===============================================================================
SUBROUTINE rotate_cube (self, vec)
  class(nlte_patch_t) :: self
  real(4), dimension(2) :: delta, deltak
  real(8) :: vec(3)
  real(4) :: bottom, bot2pivot, height, f1, f2, a, b, c, d
  integer :: up1, up2, x1, x2, y1, y2, i
  integer, save :: itimer=0
  ! ----------------------------------------------------------------------------
  call trace%begin('net_opac_t%rotate_cube', itimer=itimer)
  delta = vec(:2) / (vec(3) * self%ds(:2)) ! shift in units of grid cells

  ! bot2pivot is the distance between the bottom of the patch and the pivot plane
  bottom = self%position(3) - 0.5*self%size(3)
  bot2pivot = pivot_z * self%ds(3) - bottom

  do z=1, nz
    height = (z-1)*self%ds(3) ! height of plane as measured from bottom
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

    do i=1, 2
      do y=1, ny
        y1 = modulo(y+up2-1, ny)+1
        y2 = modulo(y+up2,   ny)+1
    
        do x=1, nx
          x1 = modulo(x+up1-1, nx)+1
          x2 = modulo(x+up1,   nx)+1

          associate(plane=>work(i,:,:,:,z), shifted=>work(3,:,:,:,z))
            shifted(:,x,y) = a*plane(:,x1,y1)
            shifted(:,x,y) = b*plane(:,x2,y1) + shifted(:,x,y)
            shifted(:,x,y) = c*plane(:,x1,y2) + shifted(:,x,y)
            shifted(:,x,y) = d*plane(:,x2,y2) + shifted(:,x,y)
          end associate
        end do
      end do
      work(i,:,:,:,z) = work(3,:,:,:,z)
    end do
  end do
  call trace%end(itimer)
END SUBROUTINE rotate_cube

!===============================================================================
SUBROUTINE short_char (self, vec, get_rates)
  class(nlte_patch_t) :: self
  type(ray_type) :: ray
  real(8) :: vec(3)
  logical :: get_rates
  real(4) :: Inuup, xxup, ssup, xxmin, ssmin, Inumin, dsup05, xx, ss, Inu, u0, Psi
  real(8) :: dt, ex ! important
  real(4), dimension(:,:), allocatable :: kernel, work_up, work_min
  integer :: dx, dy, zs, ze, xup, yup, zup, x1, y1, ii, jj
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('nlte_patch%short_char', itimer=itimer)
  call ray%init(vec)
  call ray%set_stencil(self%ds * scaling%l)

  allocate(kernel, source=real(ray%upkern))
  allocate(work_up(3,self%fs:self%fe))
  allocate(work_min(3,self%fs:self%fe))

  dsup05 = 0.5*real(ray%dsup)

  dx = size(kernel, 1)
  dy = size(kernel, 2)

  zs = merge(1, nz, self%dir > 0)
  ze = merge(nz, 1, self%dir > 0)

  do z=zs, ze, self%dir
    zup = z - self%dir
    
    do y=1, ny
      y1 = y + ray%upo2 - 1
  
      do x=1, nx
        x1 = x + ray%upo1 - 1
        
        work_up = 0
        work_min = huge(1.0)
        do jj=1,dy
          do ii=1,dx
            xup = modulo(ii+x1-1, nx) + 1
            yup = modulo(jj+y1-1, ny) + 1
            associate(wup => work(:, :, xup, yup, zup))
              work_up = work_up + kernel(ii,jj)*wup
              if (ii>1 .and. ii<dx .and. jj>1 .and. jj<dy) then
                work_min = min(work_min, wup)
              endif
            end associate
          end do
        end do
        work_up = max(work_up, work_min)

        do f=self%fs, self%fe
          xx = work(1,f,x,y,z)
          ss = work(2,f,x,y,z)
          xxup = work_up(1,f)
          ssup = work_up(2,f)
          Inuup = work_up(3,f)

          dt = dsup05 * ( xxup + xx ) ! upwind optical pathlength
          if (dt < 80) then
            ex = exp(-dt)
            if ( dt > 1d-3 ) then
              u0 = 1 - ex
              Psi = 1 - u0 / dt
            else
              Psi = dt * (0.5 - 0.16666667*dt)
              u0 = dt * (1 - Psi)
            endif

            Inu = Inuup*ex + (u0 - Psi)*ssup + Psi*ss
          else
            !ex = 0
            u0 = 1 / dt
            Psi = 1 - u0
            Inu = u0*ssup + Psi*ss
          endif

          work(3,f,x,y,z) = Inu
          if (get_rates) self%Psi(f,x,y,z) = Psi / xx
        end do ! f-loop
      end do
    end do
  end do ! z-loop
  call trace%end(itimer)
END SUBROUTINE short_char

!===============================================================================
SUBROUTINE long_char (self, vec, get_rates)
  class(nlte_patch_t) :: self
  real(8) :: vec(3)
  logical :: get_rates
  real(4) :: ds, Inuup, xxup, ssup, dsup05, xx, ss, Inu, u0, Psi
  real(8) :: dt, ex ! important
  real(4), allocatable, dimension(:,:) :: kernel
  integer :: nx, ny, nz, dx, dy, x, y, z, f, zs, ze, xup, yup, zup, x1, y1, ii, jj
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('nlte_patch%long_char', itimer=itimer)
  nx = self%patch%n(1)
  ny = self%patch%n(2)
  nz = self%patch%n(3)

  ds = self%ds(3) * scaling%l
  dsup05 = 0.5*ds/abs(vec(3))

  zs = merge(1, nz, self%dir > 0)
  ze = merge(nz, 1, self%dir > 0)

  do z=zs, ze, self%dir
    zup = z - self%dir
    do y=1, ny
    do x=1, nx
      do f=self%fs, self%fe
        xxup = work(1,f,x,y,zup)
        ssup = work(2,f,x,y,zup)
        xx = work(1,f,x,y,z)
        ss = work(2,f,x,y,z)
        dt = dsup05 * ( xxup + xx ) ! upwind optical pathlength
        if (dt < 80) then
          Inuup = work(3,f,x,y,zup)
          ex = exp(-dt)
          if ( dt > 1d-3 ) then
            u0 = 1 - ex
            Psi = 1 - u0 / dt
          else
            Psi = dt * (0.5 - 0.16666667*dt)
            u0 = dt * (1 - Psi)
          endif

          Inu = Inuup*ex + (u0 - Psi)*ssup + Psi*ss
        else
          !ex = 0
          u0 = 1 / dt
          Psi = 1 - u0
          Inu = u0*ssup + Psi*ss
        endif

        work(3,f,x,y,z) = Inu
        if (get_rates) self%Psi(f,x,y,z) = Psi / xx
      end do ! f-loop
    end do
    end do
  end do ! z-loop
  call trace%end(itimer)
END SUBROUTINE long_char

!===============================================================================
SUBROUTINE get_vproj (self, vproj, vec)
  class(nlte_patch_t) :: self
  real(4), dimension(:,:,:) :: vproj
  real(8), dimension(:) :: vec
  real(4), dimension(:,:,:), pointer :: d, px, py, pz
  integer :: li(3), ui(3)
  !-----------------------------------------------------------------------------
  if (.not. use_vel) then
    vproj = 0
  else
    li = self%patch%li
    ui = self%patch%ui

    call self%check_access(self%patch)
    d  => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),self%patch%idx%d, self%it,1)
    px => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),self%patch%idx%px,self%it,1)
    py => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),self%patch%idx%py,self%it,1)
    pz => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),self%patch%idx%pz,self%it,1)

    if (self%out_ray) then
      vproj = (-vec(1)*px - vec(2)*py + vec(3)*pz) * scaling%u / (d*cc)
    else
      vproj = (vec(1)*px + vec(2)*py + vec(3)*pz) * scaling%u / (d*cc)
    endif
  endif
END SUBROUTINE get_vproj

  !===============================================================================
SUBROUTINE print_post_update_info (self)
  class(nlte_patch_t) :: self
  !-----------------------------------------------------------------------------
  call self%print_to_stdout(arr=memcheck(1:3))
END SUBROUTINE print_post_update_info

!===============================================================================
SUBROUTINE nlte_patch_nbor_relations(self, another, needed, needs_me, download)
  class(nlte_patch_t) :: self
  class(task_t), pointer :: another        ! another task in the list
  logical :: needed, needs_me, download    ! Is the other task: needed, 
                                           ! does it need me, should I download?
  class(task_t), pointer :: other
  logical :: default
  ! ----------------------------------------------------------------------------
  call trace%begin ('nlte_patch_t%nbor_relations')
  call self%ray_patch_t%nbor_relations(another, needed, needs_me, download)

  if (needed .or. needs_me) then
    select type(another)
    class is (nlte_patch_t)
      ! do nothing
    class is (m3d_patch_t)
      if (another%id == self%m3d_parent%id) then
        needs_me = .false.
      else
        select case(another%kind)
        case ('atom_patch')
          needs_me = .true.; needed = .true.
        end select
      endif
    class default
      needed = .false.; needs_me = .false.
    end select
  endif

  call trace%end()
END SUBROUTINE nlte_patch_nbor_relations

!===============================================================================
!> Fast method to generate nbor lists in cases with very long task lists
!===============================================================================
SUBROUTINE init_nbors_hash (self, nbors)
  class(nlte_patch_t):: self
  class(nbor_t), pointer:: nbors
  !.............................................................................
  class(task_t), pointer:: task
  class(nbor_t), pointer:: nbor
  integer, save:: itimer=0
  integer :: key(4)
  !-----------------------------------------------------------------------------
  call trace%begin('nlte_patch_t%init_nbors_hash', itimer=itimer)
  call self%ray_patch_t%init_nbors_hash(nbors)

  key = self%key
  key(2) = 0
  key(3) = 0
  key(4) = hash_string('atom_patch')
  call self%check_nbor(key, nbors)

  key(4) = hash_string('composition')
  call self%check_nbor(key, nbors)
  call trace%end (itimer)
END SUBROUTINE init_nbors_hash

!===============================================================================
SUBROUTINE allocate_mem(self, m_in, label)
  class(nlte_patch_t) :: self
  class(line_type), pointer :: line
  integer, intent(in), optional :: m_in(6)
  character(len=*), optional :: label
  integer :: i, m(6)
  ! ----------------------------------------------------------------------------
  if (self%is_set(bits%no_download)) then
    m = [3, self%nf, self%n(1), self%n(2), self%n_mu, self%n_set]
    call self%patch_t%allocate_mem(m)
  endif
END SUBROUTINE allocate_mem

!===============================================================================
SUBROUTINE save_intensity (self)
  class(nlte_patch_t) :: self
  real(4), dimension(:,:,:,:,:), allocatable :: ie, i3
  integer :: meta_unit, iostat, i, nmu
  character(len=128) :: meta_file
  !-----------------------------------------------------------------------------
  fs = self%fs
  fe = self%fe
  allocate(i3(fs:fe,nx,ny, self%mus:self%mue, self%n_set), source=0.)
  allocate(ie(fs:fe, 1, 1, self%mus:self%mue, self%n_set), source=0.)

  do f=fs, fe
    i3(f,:,:,:,:) = self%twohnu3_c2(f) * self%mem_remap(3,f,:,:,:,:)
  end do

  do y=1, ny
  do x=1, nx
    ie(:,1,1,:,:) = ie(:,1,1,:,:) + i3(:,x,y,:,:)
  end do
  end do
  ie = ie / real(nx*ny)

  call self%save(ie, 'ie')
  self%saved = .false.
  call self%save(i3, 'i3')

  meta_file = trim(io%run_dir) // 'angle_quadrature.txt'
  call io%lock%set('save_quadrature')
  if (.not. os%exists(meta_file)) then
    open(newunit=meta_unit, file=meta_file, form='formatted')
    nmu = size(long_wmu)
    write(meta_unit, '(A3, 4(A12))') '#', 'weight', 'mux', 'muy', 'muz'
    do i=1, nmu
      write(meta_unit, '(I3, F12.8, SP, 3(F12.8))') i, long_wmu(i), long_vec(:,i)
    end do
    close(meta_unit)
  endif
  call io%lock%unset('save_quadrature')
END SUBROUTINE save_intensity

!===============================================================================
SUBROUTINE mem_register(arr, label)
  real(4), dimension(..) :: arr
  character(len=*), optional :: label
  !-----------------------------------------------------------------------------
  call os%mem_register(storage_size(arr), product(shape(arr)), label)
END SUBROUTINE mem_register

!===============================================================================
SUBROUTINE mem_unregister(arr, label)
  real(4), dimension(..) :: arr
  character(len=*), optional :: label
  !-----------------------------------------------------------------------------
  call os%mem_register(storage_size(arr), -product(shape(arr)), label)
END SUBROUTINE mem_unregister

END MODULE nlte_patch_mod
