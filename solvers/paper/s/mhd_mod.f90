!> =============================================================================
!> Staggered mesh solver using entropy per unit volume as the energy variable.
!>
!> Flop counting is enable via annotation at the end of lines. To count flops,
!> use a Linux pipeline such as
!>
!> % grep ops: ../../../solvers/stagger2/mhd_mod.f90 | sed 's/.*ops://' | sumcol 2
!>
!> where 'sumcol' is any script that sums up values in column 2
!> =============================================================================
MODULE mhd_mod
  USE io_mod
  USE bits_mod
  USE omp_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  USE scalar_mod
  USE vector_mod
  USE vector_ops
  USE stagger_mod
  USE gpatch_mod
  USE force_mod
  USE index_mod
  USE timestep_mod
  USE planet_mod
  USE eos_mod
  implicit none
  private
  type, extends(gpatch_t), public:: mhd_t
    real, dimension(:,:,:),   pointer:: d, dddt, s, dsdt
    real, dimension(:,:,:,:), pointer:: p, dpdt, B, dBdt
    real:: nu(6)
    real:: csound
    real:: cdtd
    logical:: mhd=.true.
    real:: gamma_AD, rho_ions ! ambipolar diffusion parameters
    real(8):: dt_Ohm, dt_AD, niMHD_courant ! non-ideal MHD time-stepping variables
    logical:: first_time=.true., mhd_AD=.false., mhd_Ohm=.false.
    type(force_t):: force
  contains
    procedure:: init
    procedure:: pde
    procedure:: nonideal_MHD_courant_condition ! non-ideal MHD replacement for `courant_condition`
    procedure:: update
    procedure:: output
    procedure:: p2U
    procedure:: U2p
    procedure:: e2E_th
    procedure:: E_th2e
    procedure:: e2S
    procedure:: S2e
    procedure:: log_density => void
    procedure:: log_pressure => void
    procedure:: gas_pressure
    procedure:: gas_pressure_ngz
    procedure:: cell_gas_pressure
    procedure:: gas_temperature
    procedure:: gas_temperature_ngz
    procedure:: cell_gas_temperature
    procedure:: velocity_magnitude => void
    procedure:: magnetic_field_magnitude => void
    procedure:: grav_potential => void
    procedure:: apply_heating
  end type
CONTAINS

!> =============================================================================
!> The construction with a double test on first_time avoids entering a critical
!> region once the variable has been set.
!> =============================================================================
SUBROUTINE init (self)
  class(mhd_t):: self
  type(index_t):: idx
  integer:: i
  real, save:: nu(6)=[0.1,1.0,0.0,0.5,0.5,0.5], gamma=1.0, csound=1.0, cdtd=0.5
  real, save:: gamma_AD=75., rho_ions=1., niMHD_courant=0.1
  logical, save:: first_time=.true., mhd_AD=.false., mhd_Ohm=.false.
  character(len=16), save:: eos
  namelist /mhd_params/ nu, gamma, csound, cdtd, mhd_AD, mhd_Ohm, gamma_AD, &
    rho_ions, niMHD_courant, hardwire, eos
  !----------------------------------------------------------------------------
  call trace%begin('mhd_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    eos = self%eos
    rewind (io%input)
    read (io%input, mhd_params)
    if (io%master) write (*, mhd_params)
    call timestep%init
  end if
  !$omp end critical (input_cr)
  self%eos = eos
  self%nu = nu
  self%gamma = gamma
  self%csound = csound
  self%cdtd = cdtd
  if (self%mhd) then
    self%kind = 'stagger2_mhd_patch'
  else
    self%kind = 'stagger2_hd_patch'
  end if
  if (mhd_Ohm) then
    self%mhd_Ohm = mhd_Ohm
    self%kind = 'stagger2_mhd_Ohm_patch'
    self%niMHD_courant = niMHD_courant
  end if
  if (mhd_AD) then
    self%mhd_AD = mhd_AD
    self%gamma_AD  = gamma_AD
    self%rho_ions  = rho_ions
    self%kind = 'stagger2_mhd_AD_patch'
    self%niMHD_courant = niMHD_courant
  endif
  call self%gpatch_t%init
  self%i_d  = idx%d
  self%i_e  = idx%e
  self%i_s  = idx%s
  self%i_px = idx%px
  self%i_py = idx%py
  self%i_pz = idx%pz
  self%i_bx = idx%bx
  self%i_by = idx%by
  self%i_bz = idx%bz
  self%type = 'mhd_t'
  do i=1,3
    self%mesh(i)%h(2+i) = -0.5
    if (self%mhd) self%mesh(i)%h(5+i) = -0.5
  end do
  self%unsigned(idx%d)  = .true.
  self%pervolume(idx%s) = .true.
  call self%force%init(self%kind)
  call trace%end
END SUBROUTINE init

!> =============================================================================
SUBROUTINE pde (self)
  USE scalar_mod
  class(mhd_t):: self
  !.............................................................................
  real, dimension(:,:,:,:), pointer:: U, E, J, p, B, dpdt, dBdt
  real, dimension(:,:,:,:), pointer:: Ud, E_AD
  real, dimension(:,:,:,:), pointer:: fx, fy, fz, fd, ld, dd, fm
  real, dimension(:,:,:),   pointer:: fmx, fmy, fmz
  real, dimension(:,:,:),   pointer:: w1, w2, w3, w4, w5, w6
  real, dimension(:,:,:),   pointer:: d, s, dddt, dsdt
  real, dimension(:,:,:),   pointer:: du, lnd, ss, cs, pg, pa, pp, pb, Q, umx
  real, dimension(:,:,:),   pointer:: dpdtx, dpdty, dpdtz, dBdtx, dBdty, dBdtz
  real, dimension(:,:,:),   pointer:: Bx, By, Bz, px, py, pz
  real, dimension(:,:,:),   pointer:: Ux, Uy, Uz, Ex, Ey, Ez, Jx, Jy, Jz
  real, dimension(:,:,:),   pointer:: fdx, fdy, fdz, ldx, ldy, ldz, ddx, ddy, ddz
  real, dimension(:,:,:),   pointer:: fxx, fxy
  real, dimension(:,:,:),   pointer::      fyy, fyz
  real, dimension(:,:,:),   pointer:: fzx,      fzz
  real, dimension(:,:,:),   pointer:: Ud_x, Ud_y, Ud_z, eta_non, eta_th, eta_Ohm, ambi, H_AD
  real, dimension(:,:,:),   pointer:: Txx, Txy
  real, dimension(:,:,:),   pointer::      Tyy, Tyz
  real, dimension(:,:,:),   pointer:: Tzx,      Tzz
  integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
  logical, save:: first_time=.true., first_AD=.true.
  !.............................................................................
  integer:: ix, iy, iz, mx, my, mz
  real(8):: ds(3), dsmin, cx, cy, cz
  integer, save:: itimer=0
  type(index_t):: i
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%update',itimer=itimer)
  self%d => self%mem(:,:,:,i%d      ,self%it,1); self%dddt => self%mem(:,:,:,i%d      ,self%it,2)
  self%s => self%mem(:,:,:,i%s      ,self%it,1); self%dsdt => self%mem(:,:,:,i%s      ,self%it,2)
  self%p => self%mem(:,:,:,i%px:i%pz,self%it,1); self%dpdt => self%mem(:,:,:,i%px:i%pz,self%it,2)
  d=>self%d; dddt=>self%dddt; s=>self%s; dsdt=>self%dsdt;
  p=>self%p; dpdt=>self%dpdt
  px=>p(:,:,:,1); py=>p(:,:,:,2); pz=>p(:,:,:,3);
  dpdtx=>dpdt(:,:,:,1); dpdty=>dpdt(:,:,:,2); dpdtz=>dpdt(:,:,:,3);
  if (self%mhd) then
    self%B => self%mem(:,:,:,i%bx:i%bz,self%it,1); self%dBdt => self%mem(:,:,:,i%bx:i%bz,self%it,2)
    B=>self%B; dBdt=>self%dBdt
    Bx=>B(:,:,:,1); By=>B(:,:,:,2); Bz=>B(:,:,:,3);
    dBdtx=>dBdt(:,:,:,1); dBdty=>dBdt(:,:,:,2); dBdtz=>dBdt(:,:,:,3);
  end if
  call allocate_vectors (self%gn, ld, dd, E, J, fd, fx, fy, fz, fm)
  call allocate_scalars (self%gn, lnd, pa, du, ss, cs, pp, pb, Q)
  call allocate_scalars (self%gn, umx, w1, w2, w3, w4, w5, w6)
  if (self%mhd_Ohm) call allocate_scalars (self%gn, eta_Ohm)
  if (self%do_pebbles) then
    U => self%vgas(:,:,:,:,self%it)
    pg => self%pressure(:,:,:,self%it)
    self%density => self%mem(:,:,:,1,:,1)
  else
    call allocate_scalars (self%gn, pg)
    call allocate_vectors (self%gn, U)
  end if
  Ux => U(:,:,:,1); Uy => U(:,:,:,2); Uz => U(:,:,:,3)
  fmx=>fm(:,:,:,1); fmy=>fm(:,:,:,2); fmz=>fm(:,:,:,3)
  if (self%mhd_AD) then
    call allocate_vectors (self%gn, Ud, E_AD)
    Ud_x=>Ud(:,:,:,1)
    Ud_y=>Ud(:,:,:,2)
    Ud_z=>Ud(:,:,:,3)
    call allocate_scalars (self%gn, ambi, H_AD)
  end if
  Ex=>E(:,:,:,1); Ey=>E(:,:,:,2); Ez=>E(:,:,:,3)
  Jx=>J(:,:,:,1); Jy=>J(:,:,:,2); Jz=>J(:,:,:,3)
  fxx=>fx(:,:,:,1); fxy=>fx(:,:,:,2)
                    fyy=>fy(:,:,:,2); fyz=>fy(:,:,:,3)
  fzx=>fz(:,:,:,1);                   fzz=>fz(:,:,:,3)
  fdx=>fd(:,:,:,1); fdy=>fd(:,:,:,2); fdz=>fd(:,:,:,3)
  ldx=>ld(:,:,:,1); ldy=>ld(:,:,:,2); ldz=>ld(:,:,:,3)
  ddx=>dd(:,:,:,1); ddy=>dd(:,:,:,2); ddz=>dd(:,:,:,3)
  call allocate_scalars (self%gn, Txx, Tyy, Tzz)
  call allocate_scalars (self%gn, Txy, Tyz, Tzx)
  !-----------------------------------------------------------------------------
  ! Stagger operator macros
  !-----------------------------------------------------------------------------
  mx = self%mesh(1)%gn
  my = self%mesh(2)%gn
  mz = self%mesh(3)%gn
  ds = self%ds
  cx = 1./ds(1)
  cy = 1./ds(2)
  cz = 1./ds(3)
  dsmin = minval(ds)
#define SECOND
#ifdef SECOND
#define mXDN(f,ix,iy,iz)    (0.5*(f(ix,iy,iz)+f(max(1,ix-1),iy,iz)))
#define mYDN(f,ix,iy,iz)    (0.5*(f(ix,iy,iz)+f(ix,max(1,iy-1),iz)))
#define mZDN(f,ix,iy,iz)    (0.5*(f(ix,iy,iz)+f(ix,iy,max(1,iz-1))))
#define mXUP(f,ix,iy,iz)    (0.5*(f(ix,iy,iz)+f(min(mx,ix+1),iy,iz)))
#define mYUP(f,ix,iy,iz)    (0.5*(f(ix,iy,iz)+f(ix,min(my,iy+1),iz)))
#define mZUP(f,ix,iy,iz)    (0.5*(f(ix,iy,iz)+f(ix,iy,min(mz,iz+1))))
#define mDDXDN(f,ix,iy,iz)  (cx*(f(ix,iy,iz)-f(max(1,ix-1),iy,iz)))
#define mDDYDN(f,ix,iy,iz)  (cy*(f(ix,iy,iz)-f(ix,max(1,iy-1),iz)))
#define mDDZDN(f,ix,iy,iz)  (cz*(f(ix,iy,iz)-f(ix,iy,max(1,iz-1))))
#define mDDXUP(f,ix,iy,iz)  (cx*(f(min(mx,ix+1),iy,iz)-f(ix,iy,iz)))
#define mDDYUP(f,ix,iy,iz)  (cy*(f(ix,min(my,iy+1),iz)-f(ix,iy,iz)))
#define mDDZUP(f,ix,iy,iz)  (cz*(f(ix,iy,min(mz,iz+1))-f(ix,iy,iz)))
#endif
  !-----------------------------------------------------------------------------
  ! Compute velocity U=momentum/density, leaving down shifted density and log
  ! density for reuse below.  U is valid in (2:gn)
  !-----------------------------------------------------------------------------
  do iz=1,mz; do iy=1,my; do ix=1,mx
    lnd(ix,iy,iz)   =    log(d(ix,iy,iz))
    ld (ix,iy,iz,1) = mXDN(lnd,ix,iy,iz)
    ld (ix,iy,iz,2) = mYDN(lnd,ix,iy,iz)
    ld (ix,iy,iz,3) = mZDN(lnd,ix,iy,iz)
    dd (ix,iy,iz,1) =   exp(ld(ix,iy,iz,1))
    dd (ix,iy,iz,2) =   exp(ld(ix,iy,iz,2))
    dd (ix,iy,iz,3) =   exp(ld(ix,iy,iz,3))
    U  (ix,iy,iz,1) =       px(ix,iy,iz)/dd(ix,iy,iz,1)
    U  (ix,iy,iz,2) =       py(ix,iy,iz)/dd(ix,iy,iz,2)
    U  (ix,iy,iz,3) =       pz(ix,iy,iz)/dd(ix,iy,iz,3)
    ss (ix,iy,iz)   =        s(ix,iy,iz)/ d(ix,iy,iz)
  end do; end do; end do
  !-----------------------------------------------------------------------------
  ! Gas pressure. The 2nd table variable is entropy per unit mass, effectively
  ! a log variable, from the table point of view.
  !-----------------------------------------------------------------------------
  pg = self%gas_pressure(lnd,ss)
  !-----------------------------------------------------------------------------
  ! Divergence of velocity, and artificial pressure.  The divergence is valid
  ! in (2,gn-1)
  !-----------------------------------------------------------------------------
  self%u_max = 0.0
  do iz=1,mz; do iy=1,my; do ix=1,mx
    du(ix,iy,iz) = -dsmin*min(0.0, mDDXUP(Ux,ix,iy,iz) &
                                 + mDDYUP(Uy,ix,iy,iz) &
                                 + mDDZUP(Uz,ix,iy,iz) )
    !---------------------------------------------------------------------------
    ! Conservation of momentum -- first the diagonal part of the stress tensor
    ! The diagonal fxx is valid in (2:gn-1)
    !---------------------------------------------------------------------------
    Txx(ix,iy,iz) = mDDXUP(Ux,ix,iy,iz)
    Tyy(ix,iy,iz) = mDDYUP(Uy,ix,iy,iz)
    Tzz(ix,iy,iz) = mDDZUP(Uz,ix,iy,iz)
    ! Diagonal Reynolds stress
!    fxx = xup(Ux)**2
!    fyy = yup(Uy)**2
!    fzz = zup(Uz)**2
!    pp = (fxx+fyy+fzz)
    fxx(ix,iy,iz) = mXUP(Ux,ix,iy,iz)**2
    fyy(ix,iy,iz) = mYUP(Uy,ix,iy,iz)**2
    fzz(ix,iy,iz) = mZUP(Uz,ix,iy,iz)**2
    pp(ix,iy,iz) = (fxx(ix,iy,iz)+fyy(ix,iy,iz)+fzz(ix,iy,iz))
    ! Diagonal (negative) Maxwell stress (borrowing E for scratch)
!    Ex = xup(Bx)**2                                                             ! mhd: 0 1 0
!    Ey = yup(By)**2                                                             ! mhd: 0 1 0
!    Ez = zup(Bz)**2                                                             ! mhd: 0 1 0
    Ex(ix,iy,iz) = mXUP(Bx,ix,iy,iz)**2
    Ey(ix,iy,iz) = mYUP(By,ix,iy,iz)**2
    Ez(ix,iy,iz) = mZUP(Bz,ix,iy,iz)**2
!    ! Magnetic pressure
!    pb = 0.5*(Ex+Ey+Ez)                                                         ! mhd: 2 1 0
!    ! Fast mode speed + velocity
!    pa = self%nu(2)*d*du**2
!    cs = sqrt((self%gamma*pg+2.0*pb + pa)/d)                                    ! mhd: 1 1 0
    ! Magnetic pressure
    pb(ix,iy,iz) = 0.5*(Ex(ix,iy,iz)+Ey(ix,iy,iz)+Ez(ix,iy,iz))
    ! Fast mode speed + velocity
    pa(ix,iy,iz) = self%nu(2)*d(ix,iy,iz)*du(ix,iy,iz)**2
    cs(ix,iy,iz) = sqrt((self%gamma*pg(ix,iy,iz)+2.0*pb(ix,iy,iz) + pa(ix,iy,iz))/d(ix,iy,iz))
    umx(ix,iy,iz) = cs(ix,iy,iz) + sqrt(pp(ix,iy,iz)) + du(ix,iy,iz)
!    ! Isotropic pressure part
!    pp = pg + pa + pb                                                           ! mhd: 1 0 0
    ! Isotropic pressure part
    pp(ix,iy,iz) = pg(ix,iy,iz) + pa(ix,iy,iz) + pb(ix,iy,iz)
!    ! Total diagonal stress
!    du = (dsmin*self%nu(1))*cs
!    fxx = pp + d*(fxx - du*Txx) - Ex                                            ! mhd: 1 0 0
!    fyy = pp + d*(fyy - du*Tyy) - Ey                                            ! mhd: 1 0 0
!    fzz = pp + d*(fzz - du*Tzz) - Ez                                            ! mhd: 1 0 0
    ! Total diagonal stress
    du(ix,iy,iz) = (dsmin*self%nu(1))*cs(ix,iy,iz)
    fxx(ix,iy,iz) = pp(ix,iy,iz) + d(ix,iy,iz)*(fxx(ix,iy,iz) - du(ix,iy,iz)*Txx(ix,iy,iz)) - Ex(ix,iy,iz)
    fyy(ix,iy,iz) = pp(ix,iy,iz) + d(ix,iy,iz)*(fyy(ix,iy,iz) - du(ix,iy,iz)*Tyy(ix,iy,iz)) - Ey(ix,iy,iz)
    fzz(ix,iy,iz) = pp(ix,iy,iz) + d(ix,iy,iz)*(fzz(ix,iy,iz) - du(ix,iy,iz)*Tzz(ix,iy,iz)) - Ez(ix,iy,iz)
!  else
!    ! Diagonal Reynolds stress
!    fxx = xup(Ux)**2                                                            ! ops: 0 1 0
!    fyy = yup(Uy)**2                                                            ! ops: 0 1 0
!    fzz = zup(Uz)**2                                                            ! ops: 0 1 0
!    pp = (fxx+fyy+fzz)                                                          ! ops: 2 0 0
!    ! Isotropic pressure part
!    pa = self%nu(2)*d*du**2                          ! dimension: pressure      ! ops: 0 3 0
!    ! Sound speed + velocity
!    cs = sqrt((self%gamma*pg+pa)/d)                                             ! ops: 1 1 1 1 0
!    if (io%verbose>3) then
!      self%u_max = self%fmaxval(cs)
!      print*,'u_max',self%u_max
!      self%u_max = self%fmaxval(cs + sqrt(pp))
!      print*,'u_max',self%u_max
!    end if
!    self%u_max = self%fmaxval(cs + sqrt(pp) + du)                               ! ops: 2 0 0 1 0
!    if (io%verbose>3) print*,'u_max',self%u_max
!    pp = pg + pa                                                                ! ops: 1 0 0
!    ! Total diagondal stress, with density factor
!    du = (dsmin*self%nu(1))*cs                                                  ! ops: 0 1 0
!    !---------------------------------------------------------------------------
!    ! Incoming fxx is valid in (2:gn-1), as is du and ddxup(Ux)
!    !---------------------------------------------------------------------------
!    fxx = pp + d*(fxx - du*Txx)                                                 ! ops: 2 2 0
!    fyy = pp + d*(fyy - du*Tyy)                                                 ! ops: 2 2 0
!    fzz = pp + d*(fzz - du*Tzz)                                                 ! ops: 2 2 0
!  end if
  !self%u_max = self%fmaxval(cs + sqrt(pp) + du)
  !-----------------------------------------------------------------------------
  ! du has units area per unit time, which gives velocity squared times mass
  ! density per unit time = energy per unit volume and time
  !-----------------------------------------------------------------------------
!  Q = d*du*(Txx**2 + Tyy**2 + Tzz**2)                                           ! ops: 2 5 0
  Q(ix,iy,iz) = d(ix,iy,iz)*du(ix,iy,iz)*(Txx(ix,iy,iz)**2 + Tyy(ix,iy,iz)**2 + Tzz(ix,iy,iz)**2)
  !-----------------------------------------------------------------------------
  ! Conservation of mass
  !-----------------------------------------------------------------------------
!  fd = -(self%nu(4)*dd*down(du))*ddown(ds,lnd)          ! diffusive mass flux   ! ops: 0 9 0
  fd(ix,iy,iz,1) = -(self%nu(4)*dd(ix,iy,iz,1)*mXDN(du,ix,iy,iz))*mDDXDN(lnd,ix,iy,iz)
  fd(ix,iy,iz,2) = -(self%nu(4)*dd(ix,iy,iz,2)*mYDN(du,ix,iy,iz))*mDDYDN(lnd,ix,iy,iz)
  fd(ix,iy,iz,3) = -(self%nu(4)*dd(ix,iy,iz,3)*mZDN(du,ix,iy,iz))*mDDZDN(lnd,ix,iy,iz)
!  fm = p+fd                                             ! mass flux             ! ops: 3 0 0
  fm(ix,iy,iz,1) = p(ix,iy,iz,1) + fd(ix,iy,iz,1)
  fm(ix,iy,iz,2) = p(ix,iy,iz,2) + fd(ix,iy,iz,2)
  fm(ix,iy,iz,3) = p(ix,iy,iz,3) + fd(ix,iy,iz,3)
  end do; end do; end do
!  dddt = -div(ds,fm)                                    ! density derivative
  do iz=1,mz; do iy=1,my; do ix=1,mx
  dddt(ix,iy,iz) = - mDDXUP(fmx,ix,iy,iz) &
                   - mDDYUP(fmy,ix,iy,iz) &
                   - mDDZUP(fmz,ix,iy,iz)
  !-----------------------------------------------------------------------------
  ! Off-diagonal part of Reynolds stress, for now without the density factor.
  ! Note that, because of the symmmetry, we define only three components
  !-----------------------------------------------------------------------------
!  fxy = ydn(Ux)*xdn(Uy)                                                         ! ops: 0 1 0
!  fyz = zdn(Uy)*ydn(Uz)                                                         ! ops: 0 1 0
!  fzx = xdn(Uz)*zdn(Ux)                                                         ! ops: 0 1 0
  fxy(ix,iy,iz) = mYDN(Ux,ix,iy,iz)*mXDN(Uy,ix,iy,iz)
  fyz(ix,iy,iz) = mZDN(Uy,ix,iy,iz)*mYDN(Uz,ix,iy,iz)
  fzx(ix,iy,iz) = mXDN(Uz,ix,iy,iz)*mZDN(Ux,ix,iy,iz)
  !-----------------------------------------------------------------------------
  ! Add viscous part -- edge-centered
  !-----------------------------------------------------------------------------
!  Txy = ddxdn(ds,Uy)+ddydn(ds,Ux)                                               ! ops: 0 1 0
!  Tyz = ddydn(ds,Uz)+ddzdn(ds,Uy)                                               ! ops: 0 1 0
!  Tzx = ddzdn(ds,Ux)+ddxdn(ds,Uz)                                               ! ops: 0 1 0
  Txy(ix,iy,iz) = mDDXDN(Uy,ix,iy,iz)+mDDYDN(Ux,ix,iy,iz)
  Tyz(ix,iy,iz) = mDDYDN(Uz,ix,iy,iz)+mDDZDN(Uy,ix,iy,iz)
  Tzx(ix,iy,iz) = mDDZDN(Ux,ix,iy,iz)+mDDXDN(Uz,ix,iy,iz)
  end do; end do; end do
!  fxy = fxy - self%nu(3)*0.5*xdn1(ydn1(du))*Txy                                 ! ops: 1 3 0
!  fyz = fyz - self%nu(3)*0.5*ydn1(zdn1(du))*Tyz                                 ! ops: 1 3 0
!  fzx = fzx - self%nu(3)*0.5*zdn1(xdn1(du))*Tzx                                 ! ops: 1 3 0
  do iz=1,mz; do iy=1,my; do ix=1,mx
    w1(ix,iy,iz) = mYDN(du,ix,iy,iz)
    w2(ix,iy,iz) = mZDN(du,ix,iy,iz)
    w3(ix,iy,iz) = mXDN(du,ix,iy,iz)
  end do; end do; end do
  do iz=1,mz; do iy=1,my; do ix=1,mx
    fxy(ix,iy,iz) = fxy(ix,iy,iz) - self%nu(3)*0.5*mXDN(w1,ix,iy,iz)*Txy(ix,iy,iz)
    fyz(ix,iy,iz) = fyz(ix,iy,iz) - self%nu(3)*0.5*mYDN(w2,ix,iy,iz)*Tyz(ix,iy,iz)
    fzx(ix,iy,iz) = fzx(ix,iy,iz) - self%nu(3)*0.5*mZDN(w3,ix,iy,iz)*Tzx(ix,iy,iz)
  end do; end do; end do
!  Q = Q + self%nu(3)*0.5*d*du*(xup(yup(Txy))**2 + &                             ! ops: 2 5 0
!                               yup(zup(Tyz))**2 + &                             ! ops: 1 1 0
!                               zup(xup(Tzx))**2)                                ! ops: 1 0 0
  do iz=1,mz; do iy=1,my; do ix=1,mx
    w1(ix,iy,iz) = mYUP(Txy,ix,iy,iz)
    w2(ix,iy,iz) = mZUP(Tyz,ix,iy,iz)
    w3(ix,iy,iz) = mXUP(Tzx,ix,iy,iz)
  end do; end do; end do
  do iz=1,mz; do iy=1,my; do ix=1,mx
    Q(ix,iy,iz) = Q(ix,iy,iz) + self%nu(3)*0.5*d(ix,iy,iz)*du(ix,iy,iz)*( &
                                 mXUP(w1,ix,iy,iz)**2 + &
                                 mYUP(w2,ix,iy,iz)**2 + &
                                 mZUP(w3,ix,iy,iz)**2)
  end do; end do; end do
  !-----------------------------------------------------------------------------
  ! Multiply with edge-centered density
  !-----------------------------------------------------------------------------
!  fxy = exp(ydn(ldx))*fxy                                                       ! ops: 0 1 0 0 1
!  fyz = exp(zdn(ldy))*fyz                                                       ! ops: 0 1 0 0 1
!  fzx = exp(xdn(ldz))*fzx                                                       ! ops: 0 1 0 0 1
  do iz=1,mz; do iy=1,my; do ix=1,mx
    fxy(ix,iy,iz) = exp(mYDN(ldx,ix,iy,iz))*fxy(ix,iy,iz)
    fyz(ix,iy,iz) = exp(mZDN(ldy,ix,iy,iz))*fyz(ix,iy,iz)
    fzx(ix,iy,iz) = exp(mXDN(ldz,ix,iy,iz))*fzx(ix,iy,iz)
  end do; end do; end do
  !-----------------------------------------------------------------------------
  ! Off-diagonal Maxwell stress
  !-----------------------------------------------------------------------------
  if (self%mhd) then
!    fxy = fxy - ydn(Bx)*xdn(By)                                                 ! mhd: 1 1 0
!    fyz = fyz - zdn(By)*ydn(Bz)                                                 ! mhd: 1 1 0
!    fzx = fzx - xdn(Bz)*zdn(Bx)                                                 ! mhd: 1 1 0
    do iz=1,mz; do iy=1,my; do ix=1,mx
      fxy(ix,iy,iz) = fxy(ix,iy,iz) - mYDN(Bx,ix,iy,iz)*mXDN(By,ix,iy,iz)
      fyz(ix,iy,iz) = fyz(ix,iy,iz) - mZDN(By,ix,iy,iz)*mYDN(Bz,ix,iy,iz)
      fzx(ix,iy,iz) = fzx(ix,iy,iz) - mXDN(Bz,ix,iy,iz)*mZDN(Bx,ix,iy,iz)
    end do; end do; end do
  end if
  !-----------------------------------------------------------------------------
  ! Consistent addition of mass flux related momentum flux, symmetrized. The
  ! diagonal part is analogous, a centered version of \nu \rho U_i \grad\ln\rho,
  ! while the off-diagonal part needs to be properly symmetrized.
  !-----------------------------------------------------------------------------
  if (self%nu(4) > 0.0) then
!    pa = d*du*self%nu(4)                                                        ! ops: 0 2 0
  do iz=1,mz; do iy=1,my; do ix=1,mx
    pa(ix,iy,iz) = d(ix,iy,iz)*du(ix,iy,iz)*self%nu(4)
    !---------------------------------------------------------------------------
    ! fxx, pa, ddxup(ldx), and xup(Ux) are all valid in (2:gn-1)
    !---------------------------------------------------------------------------
!    fxx = fxx - pa*ddxup(ds,ldx)*xup(Ux)                                        ! ops: 1 2 0
!    fyy = fyy - pa*ddyup(ds,ldy)*yup(Uy)                                        ! ops: 1 2 0
!    fzz = fzz - pa*ddzup(ds,ldz)*zup(Uz)                                        ! ops: 1 2 0
    fxx(ix,iy,iz) = fxx(ix,iy,iz) - pa(ix,iy,iz)*mDDXUP(ldx,ix,iy,iz)*mXUP(Ux,ix,iy,iz)
    fyy(ix,iy,iz) = fyy(ix,iy,iz) - pa(ix,iy,iz)*mDDYUP(ldy,ix,iy,iz)*mYUP(Uy,ix,iy,iz)
    fzz(ix,iy,iz) = fzz(ix,iy,iz) - pa(ix,iy,iz)*mDDZUP(ldz,ix,iy,iz)*mZUP(Uz,ix,iy,iz)
!    fxy = fxy + (xdn1(fdy)*ydn(Ux)+ydn1(fdx)*xdn(Uy))                           ! ops: 2 2 0
!    fyz = fyz + (ydn1(fdz)*zdn(Uy)+zdn1(fdy)*ydn(Uz))                           ! ops: 2 2 0
!    fzx = fzx + (zdn1(fdx)*xdn(Uz)+xdn1(fdz)*zdn(Ux))                           ! ops: 2 2 0
    fxy(ix,iy,iz) = fxy(ix,iy,iz) + (mXDN(fdy,ix,iy,iz)*mYDN(Ux,ix,iy,iz)+mYDN(fdx,ix,iy,iz)*mXDN(Uy,ix,iy,iz))
    fyz(ix,iy,iz) = fyz(ix,iy,iz) + (mYDN(fdz,ix,iy,iz)*mZDN(Uy,ix,iy,iz)+mZDN(fdy,ix,iy,iz)*mYDN(Uz,ix,iy,iz))
    fzx(ix,iy,iz) = fzx(ix,iy,iz) + (mZDN(fdx,ix,iy,iz)*mXDN(Uz,ix,iy,iz)+mXDN(fdz,ix,iy,iz)*mZDN(Ux,ix,iy,iz))
  end do; end do; end do
  end if
  !-----------------------------------------------------------------------------
  ! Momentum equations.  The ddxdn(fxx) part is valid in (3:gn-1), and hence px
  ! only needs guard zones in (1:2) and (gn:gn)
  !-----------------------------------------------------------------------------
!  dpdtx = - ddxdn1(ds,fxx) - ddyup1(ds,fxy) - ddzup1(ds,fzx)                    ! ops: 2 0 0
!  dpdty = - ddxup1(ds,fxy) - ddydn1(ds,fyy) - ddzup1(ds,fyz)                    ! ops: 2 0 0
!  dpdtz = - ddxup1(ds,fzx) - ddyup1(ds,fyz) - ddzdn1(ds,fzz)                    ! ops: 2 0 0
  do iz=1,mz; do iy=1,my; do ix=1,mx
    dpdtx(ix,iy,iz) = - mDDXDN(fxx,ix,iy,iz) - mDDYUP(fxy,ix,iy,iz) - mDDZUP(fzx,ix,iy,iz)
    dpdty(ix,iy,iz) = - mDDXUP(fxy,ix,iy,iz) - mDDYDN(fyy,ix,iy,iz) - mDDZUP(fyz,ix,iy,iz)
    dpdtz(ix,iy,iz) = - mDDXUP(fzx,ix,iy,iz) - mDDYUP(fyz,ix,iy,iz) - mDDZDN(fzz,ix,iy,iz)
  end do; end do; end do
  !-----------------------------------------------------------------------------
  ! External force; e.g. point source gravity.  Note that the density factor is
  ! applied here, and should not be included in the function definition
  !-----------------------------------------------------------------------------
  if (associated(self%force%selected)) &
   dpdt = dpdt + dd*self%force%selected (self%time, d, p, Ux, Uy, Uz, self%mesh)! ops: 1 0 0
  !-----------------------------------------------------------------------------
  ! Conservation of entropy.  The formulation that works best so far is to make
  ! the diffusive entropy flux directly proportional to the diffusive mass flux,
  ! plus a contribution proportional to the gradient of specific entropy ss.
  !-----------------------------------------------------------------------------
  if (self%gamma==1.0) then
    dsdt = 0.0
  else
!    fm = fm*down(ss) - (self%nu(5)*dd*down(du))*ddown(ds,ss)                    ! ops: 1 4 0
    do iz=1,mz; do iy=1,my; do ix=1,mx
      fmx(ix,iy,iz) = fmx(ix,iy,iz)*mXDN(ss,ix,iy,iz) &
         - (self%nu(5)*dd(ix,iy,iz,1)*mXDN(du,ix,iy,iz))*mDDXDN(ss,ix,iy,iz)
      fmy(ix,iy,iz) = fmy(ix,iy,iz)*mYDN(ss,ix,iy,iz) &
         - (self%nu(5)*dd(ix,iy,iz,1)*mYDN(du,ix,iy,iz))*mDDYDN(ss,ix,iy,iz)
      fmz(ix,iy,iz) = fmz(ix,iy,iz)*mZDN(ss,ix,iy,iz) &
         - (self%nu(5)*dd(ix,iy,iz,1)*mZDN(du,ix,iy,iz))*mDDZDN(ss,ix,iy,iz)
    end do; end do; end do
    !---------------------------------------------------------------------------
    ! Q is heating per unit volume (energy per unit volume and time), which when
    ! divided by presssure (energy per unit volume) becomes entropy (unitless)
    ! per unit mass.  Multiplied by mass density it becomes mass-entropy per
    ! unit volume, which is our units for s.
    !---------------------------------------------------------------------------
!    dsdt = -div(ds,fm) + Q*d/pg                                                 ! ops: 1 1 1
    do iz=1,mz; do iy=1,my; do ix=1,mx
      dsdt(ix,iy,iz) = - (mDDXUP(fmx,ix,iy,iz) + &
                          mDDYUP(fmy,ix,iy,iz) + &
                          mDDZUP(fmz,ix,iy,iz)) &
                       + Q(ix,iy,iz)*d(ix,iy,iz)/pg(ix,iy,iz)
    end do; end do; end do
  end if
  !-----------------------------------------------------------------------------
  ! Induction equation
  !-----------------------------------------------------------------------------
  if (self%mhd) then
    !---------------------------------------------------------------------------
    ! Electric current, resistive electric field and magnetic dissipation
    !---------------------------------------------------------------------------
!    J = curl_dn(ds,B)
    do iz=1,mz; do iy=1,my; do ix=1,mx
      Jx(ix,iy,iz) = mDDYDN(Bz,ix,iy,iz) - mDDZDN(By,ix,iy,iz)
      Jy(ix,iy,iz) = mDDZDN(Bx,ix,iy,iz) - mDDXDN(Bz,ix,iy,iz)
      Jz(ix,iy,iz) = mDDXDN(By,ix,iy,iz) - mDDYDN(Bx,ix,iy,iz)
    end do; end do; end do
    if (self%mhd_Ohm) then
      eta_Ohm = self%nu(6)
      self%dt_Ohm = self%niMHD_courant*minval(dsmin**2/eta_Ohm)
      E = eta_Ohm*J
    else
      self%dt_Ohm = 100.
!      E = self%nu(6)*edge(du)*J                                                 ! mhd: 0 6 0
      do iz=1,mz; do iy=1,my; do ix=1,mx
        w1(ix,iy,iz) = mYDN(du,ix,iy,iz)
        w2(ix,iy,iz) = mZDN(du,ix,iy,iz)
        w3(ix,iy,iz) = mXDN(du,ix,iy,iz)
      end do; end do; end do
      do iz=1,mz; do iy=1,my; do ix=1,mx
        Ex(ix,iy,iz) = self%nu(6)*mZDN(w1,ix,iy,iz)*Jx(ix,iy,iz)
        Ey(ix,iy,iz) = self%nu(6)*mXDN(w2,ix,iy,iz)*Jy(ix,iy,iz)
        Ez(ix,iy,iz) = self%nu(6)*mYDN(w3,ix,iy,iz)*Jz(ix,iy,iz)
      end do; end do; end do
    end if
    if (self%gamma/=1.0) then
!      Q = yup(zup(Jx*Ex)) + zup(xup(Jy*Ey)) + xup(yup(Jz*Ez))                   ! mhd: 2 3 0
!      dsdt = dsdt + Q*d/pg                                                      ! mhd: 1 1 1
      do iz=mz,1,-1; do iy=my,1,-1; do ix=mx,1,-1
        w1(ix,iy,iz) = Jx(ix,iy,iz)*Ex(ix,iy,iz)
        w2(ix,iy,iz) = Jy(ix,iy,iz)*Ey(ix,iy,iz)
        w3(ix,iy,iz) = Jz(ix,iy,iz)*Ez(ix,iy,iz)
        w4(ix,iy,iz) = mZUP(w1,ix,iy,iz)
        w5(ix,iy,iz) = mXUP(w2,ix,iy,iz)
        w6(ix,iy,iz) = mYUP(w3,ix,iy,iz)
        Q(ix,iy,iz)  = mYUP(w4,ix,iy,iz) &
                     + mZUP(w5,ix,iy,iz) &
                     + mXUP(w6,ix,iy,iz)
        dsdt(ix,iy,iz) = dsdt(ix,iy,iz) + Q(ix,iy,iz)*d(ix,iy,iz)/pg(ix,iy,iz)
      end do; end do; end do
    end if
    !---------------------------------------------------------------------------
    ! Advective electric field on cell edges
    !---------------------------------------------------------------------------
!    Ex = Ex + ydn(Uz)*zdn(By) - zdn(Uy)*ydn(Bz)                                 ! mhd: 2 2 0
!    Ey = Ey + zdn(Ux)*xdn(Bz) - xdn(Uz)*zdn(Bx)                                 ! mhd: 2 2 0
!    Ez = Ez + xdn(Uy)*ydn(Bx) - ydn(Ux)*xdn(By)                                 ! mhd: 2 2 0
    do iz=1,mz; do iy=1,my; do ix=1,mx
      Ex(ix,iy,iz) = Ex(ix,iy,iz) + mYDN(Uz,ix,iy,iz)*mZDN(By,ix,iy,iz) &
                                  - mZDN(Uy,ix,iy,iz)*mYDN(Bz,ix,iy,iz)
      Ey(ix,iy,iz) = Ey(ix,iy,iz) + mZDN(Ux,ix,iy,iz)*mXDN(Bz,ix,iy,iz) &
                                  - mXDN(Uz,ix,iy,iz)*mZDN(Bx,ix,iy,iz)
      Ez(ix,iy,iz) = Ez(ix,iy,iz) + mXDN(Uy,ix,iy,iz)*mYDN(Bx,ix,iy,iz) &
                                  - mYDN(Ux,ix,iy,iz)*mXDN(By,ix,iy,iz)
    end do; end do; end do
!    dBdt = -curl_up(ds,E)
    do iz=1,mz; do iy=1,my; do ix=1,mx
      dBdtx(ix,iy,iz) = - mDDYUP(Ez,ix,iy,iz) + mDDZUP(Ey,ix,iy,iz)
      dBdty(ix,iy,iz) = - mDDZUP(Ex,ix,iy,iz) + mDDXUP(Ez,ix,iy,iz)
      dBdtz(ix,iy,iz) = - mDDXUP(Ey,ix,iy,iz) + mDDYUP(Ex,ix,iy,iz)
    end do; end do; end do
    if (self%mhd_AD) then
      !---------------------------------------------------------------------------
      ! Below is the actual ambipolar drift parameter as given in equation 20 in
      ! Padoan, Zweibel and Nordlund (2000).
      ! To only write ambi, the definition in Padoan et al. is already multiplied
      ! by density**(-3/2) as it appears in the calculation of the drift velocity.
      ! The term d**1.5 is a cosnequence of the single fluid approximation and the
      ! relation n_i = K(n_n/10^5cm^-3)^k + K' (n_n/10^3 cm^-3)^-2
      ! K = 3*10^-3cm^-3, k=0.5 and K'=4.64*10^-4cm^-3.
      ! based on assuming that dust grains do NOT affect ionization balance significantly
      ! => second term falls off quickly => n_i~n_n^0.5
      ! In order to compute this parameter correctly, we need to know TEMPERATURE,
      ! RESOLUTION (not directly applicable here because we use AMR), BOX LENGTH,
      ! AVERAGE DENSITY, COSMIC-RAY IONISATION (K_i parameter:
      ! K_i = n_i*n_H**(-0.5) = x_i*n_H**(-0.5) ) and AVERAGE VOLUMETRIC FLOW RATE av_sigv:
      !---------------------------------------------------------------------------
      ! ambi = 0.28 * (0.1*T)**0.5 * (N/128d0) * (L/10d0)**(-1) * (av_n/200d0)**(-0.5)
      !      * (K_i/1e-5)**(-1) * (av_sigv/2e-9)**(-1) * d**(-1.5)
      !
      ! Alternative definition in Duffin and Pudritz (2008) assumes constant parameter
      ! for gamma_AD
      ! gamma_AD = av_sigv / (m_i + m_n) = 3.28e13 g^-1cm^3s^-1
      ! ambi = 1.4 / (gamma_AD*d**1.5)
      ! rho_i = 1.
      ambi = 1./(self%gamma_AD*d*self%rho_ions)
      Ud_x = xdn(ambi) * (zup(Jy)*xdn(zup(Bz)) - yup(Jz)*xdn(yup(By)))
      Ud_y = ydn(ambi) * (xup(Jz)*ydn(xup(Bx)) - zup(Jx)*ydn(zup(Bz)))
      Ud_z = zdn(ambi) * (yup(Jx)*zdn(yup(By)) - xup(Jy)*zdn(xup(Bx)))
      E_AD(:,:,:,1) = (Ud_y*Bz - Ud_z*By)
      E_AD(:,:,:,2) = (Ud_z*Bx - Ud_x*Bz)
      E_AD(:,:,:,3) = (Ud_x*By - Ud_y*Bx)
      dBdt(:,:,:,1) = dBdt(:,:,:,1) + (ddyup(ds, (Ud_x*By - Ud_y*Bx)) - ddzup(ds, (Ud_z*Bx - Ud_x*Bz)) )
      dBdt(:,:,:,2) = dBdt(:,:,:,2) + (ddzup(ds, (Ud_y*Bz - Ud_z*By)) - ddxup(ds, (Ud_x*By - Ud_y*Bx)) )
      dBdt(:,:,:,3) = dBdt(:,:,:,3) + (ddxup(ds, (Ud_z*Bx - Ud_x*Bz)) - ddyup(ds, (Ud_y*Bz - Ud_z*By)) )
      self%dt_AD = self%niMHD_courant*minval(self%gamma_AD*self%rho_ions*d/(Bx**2+By**2+Bz**2)*dsmin**2)
      if (self%gamma==1.0) then
        dsdt = 0.
      else
        ! ambipolar heating H_AD = rho_i*rho_n*nu_in*v_D^2. nu_in is the collisional frequency between ions and neutrals
        ! and so far this is all hidden in the constant gamma_AD.
        H_AD = self%rho_ions * d * self%gamma_AD * (xup(Ud_x)**2 + yup(Ud_y)**2 + zup(Ud_z)**2)
        dsdt = dsdt + H_AD*d/pg
      end if
    end if
  end if
  !-----------------------------------------------------------------------------
  self%u_max = self%fmaxval(umx)
  self%u_max = max(self%u_max,self%cdtd*dsmin*self%fmaxval(abs(dddt/d)))        ! ops: 0 0 1
  if (io%verbose>3) print*,'u_max',self%u_max
  !-----------------------------------------------------------------------------
  call deallocate_scalars (Txx, Tyy, Tzz)
  call deallocate_scalars (Txy, Tyz, Tzx)
  call deallocate_scalars (w1, w2, w3, w4, w5 ,w6)
  call deallocate_vectors (ld, dd, E, J, fd, fx, fy, fz, fm)
  call deallocate_scalars (lnd, pa, du, ss, cs, pp, pb, Q, umx)
  if (self%mhd_Ohm) then
    call deallocate_scalars (eta_Ohm)
  end if
  if (self%mhd_AD) then
    call deallocate_vectors (Ud, E_AD)
    call deallocate_scalars (ambi,H_AD)
  end if
  if (.not.self%do_pebbles) then
    call deallocate_vectors (U)
    call deallocate_scalars (pg)
  end if
  call trace%end (itimer)
END SUBROUTINE pde

!> =============================================================================
!> Courant condition, one way or another
!> =============================================================================
SUBROUTINE nonideal_MHD_courant_condition (self,dt_Ohm_inp,dt_AD_inp)
  class(mhd_t):: self
  real(8):: dt_AD_inp
  real(8):: dt_Ohm_inp
  !.............................................................................
  call trace%begin('mhd_t%nonideal_MHD_courant_condition')
  call self%lock%set ('nonideal_MHD_courant')
  if (io%verbose>2) &
    print*, 'nonideal_MHD_courant_condition: u_max, ds', self%u_max, minval(self%ds)
  if (self%dt_fixed > 0.0) then
    self%dtime = self%dt_fixed
    self%courant = self%dtime*self%u_max/minval(self%ds)
  else if (self%u_fixed > 0.0) then
    self%dtime = self%courant*minval(self%ds)/self%u_fixed
  else
    self%dtime = self%courant*minval(self%ds)/self%u_max
  endif
  !-----------------------------------------------------------------------------
  ! Non-ideal time step considerations; these are calculated in `pde`.
  !-----------------------------------------------------------------------------
  if (self%mhd_AD ) self%dtime = min (self%dtime, dt_AD_inp )
  if (self%mhd_Ohm) self%dtime = min (self%dtime, dt_Ohm_inp)
  !-----------------------------------------------------------------------------
  ! Prevent large (cheap) patches from taking time steps longer than the out_time
  !-----------------------------------------------------------------------------
  if (io%out_time > 0.0) &
    self%dtime = min (self%dtime, io%out_time)

  self%dt(self%it) = self%dtime

  if (io%verbose>3) &
    print'(a,i7,f12.6,f8.4,1p,e10.2)', 'nonideal_MHD_courant_condition:', &
      self%id, self%dtime, self%courant, self%u_max
  call self%lock%unset ('nonideal_MHD_courant')

  call trace%end
END SUBROUTINE nonideal_MHD_courant_condition

!> =============================================================================
!> Take a complete update step.  This can be overloaded in higher levels, e.g.
!> in experiment_mod.
!> =============================================================================
SUBROUTINE update (self)
  class(mhd_t):: self
  !-----------------------------------------------------------------------------
  ! The call to output below catches the state after both guard zone downloads
  ! and boundary conditions.  The call here is a fall back; if another call is
  ! placed at an earlier point in the task handling the call here will have no
  ! effect, since out_next has already been updated.
  !-----------------------------------------------------------------------------
  call self%output
  call self%pde
  call self%courant_condition
  if (io%verbose>2) &
    call self%check_density('before timestep')
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call self%check_density(' after timestep')
END SUBROUTINE update

!> =============================================================================
!> Add extra, MHD-related information to the patch .txt output file
!> =============================================================================
SUBROUTINE output (self, experiment_name)
  class(mhd_t):: self
  character(len=64):: experiment_name, filename
  optional:: experiment_name
  logical:: exists
  !-----------------------------------------------------------------------------
  if (self%time < self%out_next .or. .not.io%do_output .or. .not.io%do_legacy) return
  call self%gpatch_t%output (experiment_name)
  !$omp critical (out_cr)
  write (filename,'(a,i4.4,"_",i4.4,".txt")') trim(io%outputname), &
    self%id, self%iout
  inquire(file=filename, exist=exists)
  if (exists) then
    open (io%data_unit, file=trim(filename), form='formatted', status='old', &
      position='append')
    write (io%data_unit,'(8f10.4)') self%nu
    close (io%data_unit)
  end if
  !$omp end critical (out_cr)
END SUBROUTINE output

!> =====================================================================
!> Get velocities from momenta
!> =====================================================================
SUBROUTINE p2U(self, U, it)
  class(mhd_t) :: self
  real, dimension(:,:,:,:), pointer:: U
  integer:: it
  !.....................................................................
  real, dimension(:,:,:,:), pointer:: dd
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%i_d,           it,1), &
            p => self%mem(:,:,:,self%i_px:self%i_pz,it,1))
  allocate(dd(size(d,1),size(d,2),size(d,3),3))
  dd=down(log(d))
  U = p/exp(dd)
  deallocate(dd)
  end associate
END SUBROUTINE p2U

!> =====================================================================
!> Put velocities to momenta
!> =====================================================================
SUBROUTINE U2p(self, U, it)
  class(mhd_t) :: self
  real, dimension(:,:,:,:), pointer:: U
  integer:: it
  !.....................................................................
  real, dimension(:,:,:,:), pointer:: dd
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%i_d           ,it,1), &
            p => self%mem(:,:,:,self%i_px:self%i_pz,it,1))
  allocate(dd(size(d,1),size(d,2),size(d,3),3))
  dd=down(log(d))
  p = U*exp(dd)
  deallocate(dd)
  end associate
END SUBROUTINE U2p

!> =====================================================================
!> Get thermal energy per unit mass from entropy
!> =====================================================================
SUBROUTINE e2E_th(self, E_th, it)
  class(mhd_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:), pointer:: E_th
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%i_d,it,1), &
            s => self%mem(:,:,:,self%i_s,it,1))
  g1 = self%gamma-1.0
  E_th = d**g1*exp(s*g1/d)/g1
  end associate
END SUBROUTINE e2E_th

!> =====================================================================
!> Put thermal energy to entropy per unit volume
!> =====================================================================
SUBROUTINE E_th2e(self, E_th, it)
  class(mhd_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:), pointer:: E_th
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%i_d,it,1), &
            s => self%mem(:,:,:,self%i_s,it,1))
  g1 = self%gamma-1.0
  s = d*log(E_th*g1/d**g1)/g1
  end associate
END SUBROUTINE E_th2e

!> =====================================================================
!> Get thermal energy per unit mass from entropy
!> =====================================================================
SUBROUTINE e2S(self, s, it)
  class(mhd_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:), pointer:: s
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%i_d,it,1), &
            e => self%mem(:,:,:,self%i_s,it,1))
  g1 = self%gamma-1.0
  s = log(e*g1/d**self%gamma)/g1
  end associate
END SUBROUTINE e2S

!> =====================================================================
!> Get thermal energy per unit mass from entropy
!> =====================================================================
SUBROUTINE S2e(self, s, it)
  class(mhd_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:), pointer:: s
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%i_d,it,1), &
            e => self%mem(:,:,:,self%i_s,it,1))
  g1 = self%gamma-1.0
  e = d**self%gamma*exp(s*g1)/g1
  end associate
END SUBROUTINE S2e

!> =============================================================================
FUNCTION gas_pressure (self, lnd, ss) RESULT (pg)
  class(mhd_t):: self
  real, dimension(:,:,:), pointer:: lnd, ss, ee
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: pg
  !-----------------------------------------------------------------------------
  if (trim(self%eos)=='ideal') then
    self%d => self%mem(:,:,:,self%i_d,self%it,1)
    if (self%gamma==1.0) then
      pg = self%d*self%csound**2
    else
      self%s => self%mem(:,:,:,self%i_s,self%it,1)
      pg = self%d**self%gamma*exp(self%s/self%d*(self%gamma-1))
    end if
  else
    allocate (ee(self%gn(1),self%gn(2),self%gn(3)))
    call self%e2E_th(ee, self%it)
    call eos%lookup (shape(lnd), lnx=lnd, y=ee, pg=pg)
    deallocate(ee)
  end if
END FUNCTION gas_pressure

FUNCTION gas_pressure_ngz (self) RESULT (pg)
  class(mhd_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: pg
  integer :: l(3), u(3)
  !-----------------------------------------------------------------------------
  l = self%li; u = self%ui
  associate (d => self%mem(l(1):u(1),l(2):u(2),l(3):u(3), 1,self%it,1), &
             s => self%mem(l(1):u(1),l(2):u(2),l(3):u(3), 2,self%it,1))
  if (self%gamma==1.0) then
    pg = d
  else
    pg = d**self%gamma*exp(s/d*(self%gamma-1))
  end if
  end associate
END FUNCTION gas_pressure_ngz

FUNCTION cell_gas_pressure(self, c) RESULT (pg)
  class(mhd_t) :: self
  real :: pg
  integer :: c(3)
  !---------------------------------------------------------------------
  associate (d => self%mem(c(1),c(2),c(3), 1,self%it,1), &
             s => self%mem(c(1),c(2),c(3), 2,self%it,1))
  if (self%gamma==1.0) then
    pg = d
  else
    pg = d**self%gamma*exp(s/d*(self%gamma-1))
  end if
  end associate
END FUNCTION cell_gas_pressure

!> =============================================================================
!> temperature from ideal gas law, taking into account that the temperature is
!> normalised by k_B and mu (so that specific gas constant = 1)
!> =============================================================================
FUNCTION gas_temperature (self, lnd, ss) RESULT (tmp)
  class(mhd_t):: self
  real, dimension(:,:,:), pointer:: lnd, ss
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: tmp
  !-----------------------------------------------------------------------------
  associate (d => self%mem(:,:,:,  1,self%it,1))
  tmp = self%gas_pressure(lnd,ss)/d
  if (io%verbose > 2) then
    if (any(tmp<0)) then
      error stop "mhd_t%gas_temperature: T<0"
    end if
  end if
  end associate
END FUNCTION gas_temperature

FUNCTION gas_temperature_ngz (self) RESULT (tmp)
  class(mhd_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: tmp
  integer :: l(3), u(3)
  !---------------------------------------------------------------------
  l = self%li; u = self%ui
  associate (d => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),  1,self%it,1))
  tmp = self%gas_pressure_ngz()/d

  if (io%verbose > 2) then
    if (any(tmp<0)) then
      error stop "mhd_t%gas_temperature: T<0"
    end if
  end if
  end associate
END FUNCTION gas_temperature_ngz

FUNCTION cell_gas_temperature (self, c) RESULT (tmp)
  class(mhd_t):: self
  real :: tmp
  integer :: c(3)
  !-----------------------------------------------------------------------------
  associate (d => self%mem(c(1),c(2),c(3),  1,self%it,1))
  tmp = cell_gas_pressure(self,c)/d
  end associate
END FUNCTION cell_gas_temperature

!> =============================================================================
SUBROUTINE log_density (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = log(self%d)
END SUBROUTINE

!> =============================================================================
SUBROUTINE log_pressure (self, lnd, ss, v)
  class(mhd_t):: self
  real, dimension(:,:,:), pointer:: lnd, ss, v
  !-----------------------------------------------------------------------------
  v = log(self%gas_pressure(lnd,ss))
END SUBROUTINE

!> =============================================================================
SUBROUTINE velocity_magnitude (self, v)
  class (mhd_t):: self
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = norm(self%p/exp(down(log(self%d))))
END SUBROUTINE

!> =============================================================================
SUBROUTINE magnetic_field_magnitude (self, v)
  class (mhd_t):: self
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = norm(self%B)
END SUBROUTINE

!> =============================================================================
SUBROUTINE void (self, v)
  class (mhd_t):: self
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = 0.0
END SUBROUTINE

SUBROUTINE apply_heating (self, q)
  class (mhd_t):: self
  real, dimension(:,:,:):: q
END SUBROUTINE

END MODULE mhd_mod
