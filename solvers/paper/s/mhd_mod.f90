!===============================================================================
!> Staggered mesh solver using entropy per unit volume as the energy variable.
!>
!> Flop counting is enable via annotation at the end of lines. To count flops,
!> use a Linux pipeline such as
!>
!> % grep ops: ../../../solvers/stagger2/mhd_mod.f90 | sed 's/.*ops://' | sumcol 2
!>
!> where 'sumcol' is any script that sums up values in column 2
!===============================================================================
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

!===============================================================================
!> The construction with a double test on first_time avoids entering a critical
!> region once the variable has been set.
!===============================================================================
SUBROUTINE init (self)
  class(mhd_t):: self
  integer:: i
  real, save:: nu(6)=[0.1,1.0,0.0,0.5,0.5,0.5], gamma=1.0, csound=1.0, cdtd=0.5
  real, save:: gamma_AD=75., rho_ions=1., niMHD_courant=0.1
  logical, save:: first_time=.true., mhd_AD=.false., mhd_Ohm=.false.
  character(len=16), save:: eos
  namelist /mhd_params/ nu, gamma, csound, cdtd, mhd_AD, mhd_Ohm, gamma_AD, &
    rho_ions, niMHD_courant, hardwire, eos
  !----------------------------------------------------------------------------
  call trace%begin('mhd_t%init')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    eos = self%eos
    rewind (io%input)
    read (io%input, mhd_params)
    if (io%master) write (*, mhd_params)
    call timestep%init
  end if
  call io%lock%unset
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
  call self%gpatch_t%init
  self%type = 'mhd_t'
  do i=1,3
    self%mesh(i)%h(2+i) = -0.5
    if (self%mhd) self%mesh(i)%h(5+i) = -0.5
  end do
  self%unsigned(self%idx%d)  = .true.
  self%pervolume(self%idx%s) = .true.
  call self%force%init(self%kind)
  call trace%end()
END SUBROUTINE init

!===============================================================================
SUBROUTINE pde (self)
  USE scalar_mod
  class(mhd_t):: self
  !.............................................................................
  real, dimension(:,:,:,:), pointer:: U, p, B, dpdt, dBdt
  real, dimension(:,:,:,:), pointer:: Ud, E_AD
  real, dimension(:,:,:),   pointer:: fmx, fmy, fmz
  real, dimension(:,:,:),   pointer:: d, s, dddt, dsdt
  real, dimension(:,:,:),   pointer:: dpdtx, dpdty, dpdtz, dBdtx, dBdty, dBdtz
  real, dimension(:,:,:),   pointer:: Bx, By, Bz, px, py, pz
  !dir$ attributes align: 64::        fmx, fmy, fmz
  real, dimension(:,:,:),   pointer:: w1, w2, w3, w4, w5, w6
  !dir$ attributes align: 64::        w1, w2, w3, w4, w5, w6
  real, dimension(:,:,:),   pointer:: du, lnd, ss, cs, pg, pa, pp, pb, Q, umx
  !dir$ attributes align: 64::        du, lnd, ss, cs, pg, pa, pp, pb, Q, umx
  real, dimension(:,:,:),   pointer:: Ux, Uy, Uz, Ex, Ey, Ez, Jx, Jy, Jz
  !dir$ attributes align: 64::                    Ex, Ey, Ez, Jx, Jy, Jz
  real, dimension(:,:,:),   pointer:: fdx, fdy, fdz, ldx, ldy, ldz, ddx, ddy, ddz
  !dir$ attributes align: 64::        fdx, fdy, fdz, ldx, ldy, ldz, ddx, ddy, ddz
  real, dimension(:,:,:),   pointer:: fxx, fyy, fzz, fxy, fyz, fzx
  !dir$ attributes align: 64::        fxx, fyy, fzz, fxy, fyz, fzx
  real, dimension(:,:,:),   pointer:: Txx, Tyy, Tzz, Txy, Tyz, Tzx
  !dir$ attributes align: 64::        Txx, Tyy, Tzz, Txy, Tyz, Tzx
  real, dimension(:,:,:),   pointer:: Ud_x, Ud_y, Ud_z, eta_non, eta_th, eta_Ohm, ambi, H_AD
  integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
  logical, save:: first_time=.true., first_AD=.true.
  !.............................................................................
  integer:: ix, iy, iz, mx, my, mz
  real(8):: ds(3)
  real:: dsmin, cx, cy, cz
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
  mx = self%mesh(1)%gn
  my = self%mesh(2)%gn
  mz = self%mesh(3)%gn
  allocate (fxx(mx,my,mz), fyy(mx,my,mz), fzz(mx,my,mz))
  allocate (fxy(mx,my,mz), fyz(mx,my,mz), fzx(mx,my,mz))
  allocate (fdx(mx,my,mz), fdy(mx,my,mz), fdz(mx,my,mz))
  allocate (ldx(mx,my,mz), ldy(mx,my,mz), ldz(mx,my,mz))
  allocate (fmx(mx,my,mz), fmy(mx,my,mz), fmz(mx,my,mz))
  allocate (ddx(mx,my,mz), ddy(mx,my,mz), ddz(mx,my,mz))
  allocate (lnd(mx,my,mz), pa(mx,my,mz), du(mx,my,mz), ss(mx,my,mz))
  allocate ( cs(mx,my,mz), pp(mx,my,mz), pb(mx,my,mz),  Q(mx,my,mz))
  allocate ( w1(mx,my,mz), w2(mx,my,mz), w3(mx,my,mz))
  allocate ( w4(mx,my,mz), w5(mx,my,mz), w6(mx,my,mz))
  allocate (umx(mx,my,mz))
  allocate (Txx(mx,my,mz), Tyy(mx,my,mz), Tzz(mx,my,mz))
  allocate (Txy(mx,my,mz), Tyz(mx,my,mz), Tzx(mx,my,mz))
  if (self%mhd) then
    allocate (Ex(mx,my,mz),  Ey(mx,my,mz),  Ez(mx,my,mz))
    allocate (Jx(mx,my,mz),  Jy(mx,my,mz),  Jz(mx,my,mz))
  end if
  if (self%do_pebbles) then
    U => self%vgas(:,:,:,:,self%it)
    pg => self%pressure(:,:,:,self%it)
    self%density => self%mem(:,:,:,1,:,1)
  else
    allocate (pg(mx,my,mz), U(mx,my,mz,3))
  end if
  Ux => U(:,:,:,1); Uy => U(:,:,:,2); Uz => U(:,:,:,3)
  !-----------------------------------------------------------------------------
  ! Stagger operator macros
  !-----------------------------------------------------------------------------
  ds = self%ds
  cx = 1./ds(1)
  cy = 1./ds(2)
  cz = 1./ds(3)
  dsmin = minval(ds)
#define SECOND
#ifdef SECOND
#define mXDN(f,ix,iy,iz)    (0.5*(f(ix,iy,iz)+f(ix,iy,iz)))
#define mYDN(f,ix,iy,iz)    (0.5*(f(ix,iy,iz)+f(ix,iy,iz)))
#define mZDN(f,ix,iy,iz)    (0.5*(f(ix,iy,iz)+f(ix,iy,iz)))
#define mXUP(f,ix,iy,iz)    (0.5*(f(ix,iy,iz)+f(ix,iy,iz)))
#define mYUP(f,ix,iy,iz)    (0.5*(f(ix,iy,iz)+f(ix,iy,iz)))
#define mZUP(f,ix,iy,iz)    (0.5*(f(ix,iy,iz)+f(ix,iy,iz)))
#define mDDXDN(f,ix,iy,iz)  (cx*(f(ix,iy,iz)-f(ix,iy,iz)))
#define mDDYDN(f,ix,iy,iz)  (cy*(f(ix,iy,iz)-f(ix,iy,iz)))
#define mDDZDN(f,ix,iy,iz)  (cz*(f(ix,iy,iz)-f(ix,iy,iz)))
#define mDDXUP(f,ix,iy,iz)  (cx*(f(ix,iy,iz)-f(ix,iy,iz)))
#define mDDYUP(f,ix,iy,iz)  (cy*(f(ix,iy,iz)-f(ix,iy,iz)))
#define mDDZUP(f,ix,iy,iz)  (cz*(f(ix,iy,iz)-f(ix,iy,iz)))
#endif
  !-----------------------------------------------------------------------------
  ! Compute velocity U=momentum/density, leaving down shifted density and log
  ! density for reuse below.  U is valid in (2:gn)
  !-----------------------------------------------------------------------------
  do iz=2,mz; do iy=2,my
  !$omp simd linear(ix)
  do ix=2,mx
    lnd(ix,iy,iz) =    log(d(ix,iy,iz))
    ldx(ix,iy,iz) = mXDN(lnd,ix,iy,iz)
    ldy(ix,iy,iz) = mYDN(lnd,ix,iy,iz)
    ldz(ix,iy,iz) = mZDN(lnd,ix,iy,iz)
    ddx(ix,iy,iz) =   exp(ldx(ix,iy,iz))
    ddy(ix,iy,iz) =   exp(ldy(ix,iy,iz))
    ddz(ix,iy,iz) =   exp(ldz(ix,iy,iz))
    Ux (ix,iy,iz) =       px(ix,iy,iz)/ddx(ix,iy,iz)
    Uy (ix,iy,iz) =       py(ix,iy,iz)/ddy(ix,iy,iz)
    Uz (ix,iy,iz) =       pz(ix,iy,iz)/ddz(ix,iy,iz)
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
  do iz=1,mz-1; do iy=1,my-1
    !$omp simd
    do ix=1,mx-1
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
!      fxx = xup(Ux)**2
!      fyy = yup(Uy)**2
!      fzz = zup(Uz)**2
!      pp = (fxx+fyy+fzz)
      fxx(ix,iy,iz) = mXUP(Ux,ix,iy,iz)**2
      fyy(ix,iy,iz) = mYUP(Uy,ix,iy,iz)**2
      fzz(ix,iy,iz) = mZUP(Uz,ix,iy,iz)**2
      pp(ix,iy,iz) = (fxx(ix,iy,iz)+fyy(ix,iy,iz)+fzz(ix,iy,iz))
    end do
    if (self%mhd) then
      !$omp simd
      do ix=1,mx-1
        ! Diagonal (negative) Maxwell stress (borrowing E for scratch)
!        Ex = xup(Bx)**2                                                             ! mhd: 0 1 0
!        Ey = yup(By)**2                                                             ! mhd: 0 1 0
!        Ez = zup(Bz)**2                                                             ! mhd: 0 1 0
        Ex(ix,iy,iz) = mXUP(Bx,ix,iy,iz)**2
        Ey(ix,iy,iz) = mYUP(By,ix,iy,iz)**2
        Ez(ix,iy,iz) = mZUP(Bz,ix,iy,iz)**2
!        ! Magnetic pressure
!        pb = 0.5*(Ex+Ey+Ez)                                                         ! mhd: 2 1 0
!        ! Fast mode speed + velocity
!        pa = self%nu(2)*d*du**2
!        cs = sqrt((self%gamma*pg+2.0*pb + pa)/d)                                    ! mhd: 1 1 0
        ! Magnetic pressure
        pb(ix,iy,iz) = 0.5*(Ex(ix,iy,iz)+Ey(ix,iy,iz)+Ez(ix,iy,iz))
        ! Fast mode speed + velocity
        pa(ix,iy,iz) = self%nu(2)*d(ix,iy,iz)*du(ix,iy,iz)**2
        cs(ix,iy,iz) = sqrt((self%gamma*pg(ix,iy,iz)+2.0*pb(ix,iy,iz) + pa(ix,iy,iz))/d(ix,iy,iz))
      end do
    else
      !$omp simd
      do ix=1,mx-1
        pa(ix,iy,iz) = self%nu(2)*d(ix,iy,iz)*du(ix,iy,iz)**2
        cs(ix,iy,iz) = sqrt((self%gamma*pg(ix,iy,iz) + pa(ix,iy,iz))/d(ix,iy,iz))
      end do
    end if
    !$omp simd
    do ix=1,mx-1
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
  end do; end do; end do
  !-----------------------------------------------------------------------------
  ! Conservation of mass
  !-----------------------------------------------------------------------------
  do iz=2,mz; do iy=2,my
  !$omp simd
  do ix=2,mx
!  fd = -(self%nu(4)*dd*down(du))*ddown(ds,lnd)          ! diffusive mass flux   ! ops: 0 9 0
  fdx(ix,iy,iz) = -(self%nu(4)*ddx(ix,iy,iz)*mXDN(du,ix,iy,iz))*mDDXDN(lnd,ix,iy,iz)
  fdy(ix,iy,iz) = -(self%nu(4)*ddy(ix,iy,iz)*mYDN(du,ix,iy,iz))*mDDYDN(lnd,ix,iy,iz)
  fdz(ix,iy,iz) = -(self%nu(4)*ddz(ix,iy,iz)*mZDN(du,ix,iy,iz))*mDDZDN(lnd,ix,iy,iz)
!  fm = p+fd                                             ! mass flux             ! ops: 3 0 0
  fmx(ix,iy,iz) = p(ix,iy,iz,1) + fdx(ix,iy,iz)
  fmy(ix,iy,iz) = p(ix,iy,iz,2) + fdy(ix,iy,iz)
  fmz(ix,iy,iz) = p(ix,iy,iz,3) + fdz(ix,iy,iz)
  end do; end do; end do
!===============================================================================
  do iz=2,mz; do iy=2,my
  !$omp simd
  do ix=2,mx
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
!  fxy = fxy - self%nu(3)*0.5*xdn1(ydn1(du))*Txy                                 ! ops: 1 3 0
!  fyz = fyz - self%nu(3)*0.5*ydn1(zdn1(du))*Tyz                                 ! ops: 1 3 0
!  fzx = fzx - self%nu(3)*0.5*zdn1(xdn1(du))*Tzx                                 ! ops: 1 3 0
    w1(ix,iy,iz) = mYDN(du,ix,iy,iz)
    w2(ix,iy,iz) = mZDN(du,ix,iy,iz)
    w3(ix,iy,iz) = mXDN(du,ix,iy,iz)
  end do; end do; end do
!===============================================================================
  do iz=2,mz; do iy=2,my
  !$omp simd
  do ix=2,mx
    fxy(ix,iy,iz) = fxy(ix,iy,iz) - self%nu(3)*0.5*mXDN(w1,ix,iy,iz)*Txy(ix,iy,iz)
    fyz(ix,iy,iz) = fyz(ix,iy,iz) - self%nu(3)*0.5*mYDN(w2,ix,iy,iz)*Tyz(ix,iy,iz)
    fzx(ix,iy,iz) = fzx(ix,iy,iz) - self%nu(3)*0.5*mZDN(w3,ix,iy,iz)*Tzx(ix,iy,iz)
  end do; end do; end do
!===============================================================================
  do iz=1,mz-1; do iy=1,my-1
  !$omp simd
  do ix=1,mx-1
    w1(ix,iy,iz) = mYUP(Txy,ix,iy,iz)
    w2(ix,iy,iz) = mZUP(Tyz,ix,iy,iz)
    w3(ix,iy,iz) = mXUP(Tzx,ix,iy,iz)
  end do; end do; end do
!===============================================================================
  do iz=1,mz-1; do iy=1,my-1
  !$omp simd
  do ix=1,mx-1
!  Q = Q + self%nu(3)*0.5*d*du*(xup(yup(Txy))**2 + &                             ! ops: 2 5 0
!                               yup(zup(Tyz))**2 + &                             ! ops: 1 1 0
!                               zup(xup(Tzx))**2)                                ! ops: 1 0 0
    Q(ix,iy,iz) = Q(ix,iy,iz) + self%nu(3)*0.5*d(ix,iy,iz)*du(ix,iy,iz)*( &
                                 mXUP(w1,ix,iy,iz)**2 + &
                                 mYUP(w2,ix,iy,iz)**2 + &
                                 mZUP(w3,ix,iy,iz)**2)
  end do; end do; end do
  !-----------------------------------------------------------------------------
  ! Multiply with edge-centered density
  !-----------------------------------------------------------------------------
!===============================================================================
  do iz=2,mz; do iy=2,my
  !$omp simd
  do ix=2,mx
!  fxy = exp(ydn(ldx))*fxy                                                       ! ops: 0 1 0 0 1
!  fyz = exp(zdn(ldy))*fyz                                                       ! ops: 0 1 0 0 1
!  fzx = exp(xdn(ldz))*fzx                                                       ! ops: 0 1 0 0 1
    fxy(ix,iy,iz) = exp(mYDN(ldx,ix,iy,iz))*fxy(ix,iy,iz)
    fyz(ix,iy,iz) = exp(mZDN(ldy,ix,iy,iz))*fyz(ix,iy,iz)
    fzx(ix,iy,iz) = exp(mXDN(ldz,ix,iy,iz))*fzx(ix,iy,iz)
  end do; end do; end do
  !-----------------------------------------------------------------------------
  ! Off-diagonal Maxwell stress
  !-----------------------------------------------------------------------------
  if (self%mhd) then
!===============================================================================
    do iz=2,mz; do iy=2,my
    !$omp simd
    do ix=2,mx
!    fxy = fxy - ydn(Bx)*xdn(By)                                                 ! mhd: 1 1 0
!    fyz = fyz - zdn(By)*ydn(Bz)                                                 ! mhd: 1 1 0
!    fzx = fzx - xdn(Bz)*zdn(Bx)                                                 ! mhd: 1 1 0
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
!===============================================================================
  do iz=2,mz-1; do iy=2,my-1
  !$omp simd
  do ix=2,mx-1
!    pa = d*du*self%nu(4)                                                        ! ops: 0 2 0
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
!===============================================================================
  do iz=2,mz-1; do iy=2,my-1
  !$omp simd
  do ix=2,mx-1
!  dpdtx = - ddxdn1(ds,fxx) - ddyup1(ds,fxy) - ddzup1(ds,fzx)                    ! ops: 2 0 0
!  dpdty = - ddxup1(ds,fxy) - ddydn1(ds,fyy) - ddzup1(ds,fyz)                    ! ops: 2 0 0
!  dpdtz = - ddxup1(ds,fzx) - ddyup1(ds,fyz) - ddzdn1(ds,fzz)                    ! ops: 2 0 0
    dpdtx(ix,iy,iz) = - mDDXDN(fxx,ix,iy,iz) - mDDYUP(fxy,ix,iy,iz) - mDDZUP(fzx,ix,iy,iz)
    dpdty(ix,iy,iz) = - mDDXUP(fxy,ix,iy,iz) - mDDYDN(fyy,ix,iy,iz) - mDDZUP(fyz,ix,iy,iz)
    dpdtz(ix,iy,iz) = - mDDXUP(fzx,ix,iy,iz) - mDDYUP(fyz,ix,iy,iz) - mDDZDN(fzz,ix,iy,iz)
  end do; end do; end do
  !-----------------------------------------------------------------------------
  ! External force; e.g. point source gravity.  Note that the density factor is
  ! applied here, and should not be included in the function definition
  !-----------------------------------------------------------------------------
  !if (associated(self%force%selected)) &
  ! dpdt = dpdt + dd*self%force%selected (self%time, d, p, Ux, Uy, Uz, self%mesh)! ops: 1 0 0
  !-----------------------------------------------------------------------------
  ! Conservation of entropy.  The formulation that works best so far is to make
  ! the diffusive entropy flux directly proportional to the diffusive mass flux,
  ! plus a contribution proportional to the gradient of specific entropy ss.
  !-----------------------------------------------------------------------------
  if (self%gamma==1.0) then
    dsdt = 0.0
  else
!===============================================================================
    do iz=2,mz; do iy=2,my
    !$omp simd
    do ix=2,mx
!    fm = fm*down(ss) - (self%nu(5)*dd*down(du))*ddown(ds,ss)                    ! ops: 1 4 0
      fmx(ix,iy,iz) = fmx(ix,iy,iz)*mXDN(ss,ix,iy,iz) &
         - (self%nu(5)*ddx(ix,iy,iz)*mXDN(du,ix,iy,iz))*mDDXDN(ss,ix,iy,iz)
      fmy(ix,iy,iz) = fmy(ix,iy,iz)*mYDN(ss,ix,iy,iz) &
         - (self%nu(5)*ddy(ix,iy,iz)*mYDN(du,ix,iy,iz))*mDDYDN(ss,ix,iy,iz)
      fmz(ix,iy,iz) = fmz(ix,iy,iz)*mZDN(ss,ix,iy,iz) &
         - (self%nu(5)*ddz(ix,iy,iz)*mZDN(du,ix,iy,iz))*mDDZDN(ss,ix,iy,iz)
    end do; end do; end do
    !---------------------------------------------------------------------------
    ! Q is heating per unit volume (energy per unit volume and time), which when
    ! divided by presssure (energy per unit volume) becomes entropy (unitless)
    ! per unit mass.  Multiplied by mass density it becomes mass-entropy per
    ! unit volume, which is our units for s.
    !---------------------------------------------------------------------------
!===============================================================================
    do iz=2,mz-1; do iy=2,my-1
    !$omp simd
    do ix=2,mx-1
!    dsdt = -div(ds,fm) + Q*d/pg                                                 ! ops: 1 1 1
      dsdt(ix,iy,iz) = - mDDXUP(fmx,ix,iy,iz) &
                       - mDDYUP(fmy,ix,iy,iz) &
                       - mDDZUP(fmz,ix,iy,iz) &
                       + Q(ix,iy,iz)*d(ix,iy,iz)/pg(ix,iy,iz)
!    dddt = -div(ds,fm)                                    ! density derivative
      dddt(ix,iy,iz) = - mDDXUP(fmx,ix,iy,iz) &
                       - mDDYUP(fmy,ix,iy,iz) &
                       - mDDZUP(fmz,ix,iy,iz)
    end do; end do; end do
  end if
  !-----------------------------------------------------------------------------
  ! Induction equation
  !-----------------------------------------------------------------------------
  if (self%mhd) then
    !---------------------------------------------------------------------------
    ! Electric current, resistive electric field and magnetic dissipation
    !---------------------------------------------------------------------------
!===============================================================================
    self%dt_Ohm = 100.
    do iz=2,mz; do iy=2,my
    !$omp simd
    do ix=2,mx
!    J = curl_dn(ds,B)
      Jx(ix,iy,iz) = mDDYDN(Bz,ix,iy,iz) - mDDZDN(By,ix,iy,iz)
      Jy(ix,iy,iz) = mDDZDN(Bx,ix,iy,iz) - mDDXDN(Bz,ix,iy,iz)
      Jz(ix,iy,iz) = mDDXDN(By,ix,iy,iz) - mDDYDN(Bx,ix,iy,iz)
      w1(ix,iy,iz) = mYDN(du,ix,iy,iz)
      w2(ix,iy,iz) = mZDN(du,ix,iy,iz)
      w3(ix,iy,iz) = mXDN(du,ix,iy,iz)
    end do; end do; end do
!===============================================================================
    do iz=2,mz; do iy=2,my
    !$omp simd
    do ix=2,mx
!    E = self%nu(6)*edge(du)*J                                                 ! mhd: 0 6 0
      Ex(ix,iy,iz) = self%nu(6)*mZDN(w1,ix,iy,iz)*Jx(ix,iy,iz)
      Ey(ix,iy,iz) = self%nu(6)*mXDN(w2,ix,iy,iz)*Jy(ix,iy,iz)
      Ez(ix,iy,iz) = self%nu(6)*mYDN(w3,ix,iy,iz)*Jz(ix,iy,iz)
    end do; end do; end do
    if (self%gamma/=1.0) then
!===============================================================================
      do iz=1,mz-1; do iy=1,my-1
      !$omp simd
      do ix=1,mx-1
!      Q = yup(zup(Jx*Ex)) + zup(xup(Jy*Ey)) + xup(yup(Jz*Ez))                   ! mhd: 2 3 0
        w1(ix,iy,iz) = Jx(ix,iy,iz)*Ex(ix,iy,iz)
        w2(ix,iy,iz) = Jy(ix,iy,iz)*Ey(ix,iy,iz)
        w3(ix,iy,iz) = Jz(ix,iy,iz)*Ez(ix,iy,iz)
        w4(ix,iy,iz) = mZUP(w1,ix,iy,iz)
        w5(ix,iy,iz) = mXUP(w2,ix,iy,iz)
        w6(ix,iy,iz) = mYUP(w3,ix,iy,iz)
        Q(ix,iy,iz)  = mYUP(w4,ix,iy,iz) &
                     + mZUP(w5,ix,iy,iz) &
                     + mXUP(w6,ix,iy,iz)
!      dsdt = dsdt + Q*d/pg                                                      ! mhd: 1 1 1
        dsdt(ix,iy,iz) = dsdt(ix,iy,iz) + Q(ix,iy,iz)*d(ix,iy,iz)/pg(ix,iy,iz)
      end do; end do; end do
    end if
    !---------------------------------------------------------------------------
    ! Advective electric field on cell edges
    !---------------------------------------------------------------------------
!===============================================================================
    do iz=2,mz; do iy=2,my
    !$omp simd
    do ix=2,mx
!    Ex = Ex + ydn(Uz)*zdn(By) - zdn(Uy)*ydn(Bz)                                 ! mhd: 2 2 0
!    Ey = Ey + zdn(Ux)*xdn(Bz) - xdn(Uz)*zdn(Bx)                                 ! mhd: 2 2 0
!    Ez = Ez + xdn(Uy)*ydn(Bx) - ydn(Ux)*xdn(By)                                 ! mhd: 2 2 0
      Ex(ix,iy,iz) = Ex(ix,iy,iz) + mYDN(Uz,ix,iy,iz)*mZDN(By,ix,iy,iz) &
                                  - mZDN(Uy,ix,iy,iz)*mYDN(Bz,ix,iy,iz)
      Ey(ix,iy,iz) = Ey(ix,iy,iz) + mZDN(Ux,ix,iy,iz)*mXDN(Bz,ix,iy,iz) &
                                  - mXDN(Uz,ix,iy,iz)*mZDN(Bx,ix,iy,iz)
      Ez(ix,iy,iz) = Ez(ix,iy,iz) + mXDN(Uy,ix,iy,iz)*mYDN(Bx,ix,iy,iz) &
                                  - mYDN(Ux,ix,iy,iz)*mXDN(By,ix,iy,iz)
    end do; end do; end do
!===============================================================================
    do iz=2,mz-1; do iy=2,my-1
    !$omp simd
    do ix=2,mx-1
!    dBdt = -curl_up(ds,E)
      dBdtx(ix,iy,iz) = - mDDYUP(Ez,ix,iy,iz) + mDDZUP(Ey,ix,iy,iz)
      dBdty(ix,iy,iz) = - mDDZUP(Ex,ix,iy,iz) + mDDXUP(Ez,ix,iy,iz)
      dBdtz(ix,iy,iz) = - mDDXUP(Ey,ix,iy,iz) + mDDYUP(Ex,ix,iy,iz)
    end do; end do; end do
  end if
  !-----------------------------------------------------------------------------
  self%u_max = self%fmaxval(umx)
  self%u_max = max(self%u_max,self%cdtd*dsmin*self%fmaxval(abs(dddt/d)))        ! ops: 0 0 1
  if (io%verbose>3) print*,'u_max',self%u_max
  !-----------------------------------------------------------------------------
  deallocate (fxx, fyy, fzz)
  deallocate (fxy, fyz, fzx)
  deallocate (fdx, fdy, fdz)
  deallocate (ldx, ldy, ldz)
  deallocate (fmx, fmy, fmz)
  deallocate (ddx, ddy, ddz)
  deallocate (lnd, pa, du, ss)
  deallocate ( cs, pp, pb,  Q)
  deallocate ( w1, w2, w3)
  deallocate ( w4, w5, w6)
  deallocate (umx)
  deallocate (Txx, Tyy, Tzz)
  deallocate (Txy, Tyz, Tzx)
  if (self%mhd) then
    deallocate ( Ex,  Ey,  Ez)
    deallocate ( Jx,  Jy,  Jz)
  end if
  if (.not.self%do_pebbles) then
    deallocate (U)
    deallocate (pg)
  end if
  call trace%end (itimer)
END SUBROUTINE pde

!===============================================================================
!> Courant condition, one way or another
!===============================================================================
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

  call trace%end()
END SUBROUTINE nonideal_MHD_courant_condition

!===============================================================================
!> Take a complete update step.  This can be overloaded in higher levels, e.g.
!> in experiment_mod.
!===============================================================================
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

!===============================================================================
!> Add extra, MHD-related information to the patch .txt output file
!===============================================================================
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

!=======================================================================
!> Get velocities from momenta
!=======================================================================
SUBROUTINE p2U(self, U, it)
  class(mhd_t) :: self
  real, dimension(:,:,:,:), pointer:: U
  integer:: it
  !.....................................................................
  real, dimension(:,:,:,:), pointer:: dd
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,           it,1), &
            p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1))
  allocate(dd(size(d,1),size(d,2),size(d,3),3))
  dd=down(log(d))
  U = p/exp(dd)
  deallocate(dd)
  end associate
END SUBROUTINE p2U

!=======================================================================
!> Put velocities to momenta
!=======================================================================
SUBROUTINE U2p(self, U, it)
  class(mhd_t) :: self
  real, dimension(:,:,:,:), pointer:: U
  integer:: it
  !.....................................................................
  real, dimension(:,:,:,:), pointer:: dd
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d           ,it,1), &
            p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1))
  allocate(dd(size(d,1),size(d,2),size(d,3),3))
  dd=down(log(d))
  p = U*exp(dd)
  deallocate(dd)
  end associate
END SUBROUTINE U2p

!=======================================================================
!> Get thermal energy per unit mass from entropy
!=======================================================================
SUBROUTINE e2E_th(self, E_th, it)
  class(mhd_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:), pointer:: E_th
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,it,1), &
            s => self%mem(:,:,:,self%idx%s,it,1))
  g1 = self%gamma-1.0
  E_th = d**g1*exp(s*g1/d)/g1
  end associate
END SUBROUTINE e2E_th

!=======================================================================
!> Put thermal energy to entropy per unit volume
!=======================================================================
SUBROUTINE E_th2e(self, E_th, it)
  class(mhd_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:), pointer:: E_th
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,it,1), &
            s => self%mem(:,:,:,self%idx%s,it,1))
  g1 = self%gamma-1.0
  s = d*log(E_th*g1/d**g1)/g1
  end associate
END SUBROUTINE E_th2e

!=======================================================================
!> Get thermal energy per unit mass from entropy
!=======================================================================
SUBROUTINE e2S(self, s, it)
  class(mhd_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:), pointer:: s
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,it,1), &
            e => self%mem(:,:,:,self%idx%s,it,1))
  g1 = self%gamma-1.0
  s = log(e*g1/d**self%gamma)/g1
  end associate
END SUBROUTINE e2S

!=======================================================================
!> Get thermal energy per unit mass from entropy
!=======================================================================
SUBROUTINE S2e(self, s, it)
  class(mhd_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:), pointer:: s
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,it,1), &
            e => self%mem(:,:,:,self%idx%s,it,1))
  g1 = self%gamma-1.0
  e = d**self%gamma*exp(s*g1)/g1
  end associate
END SUBROUTINE S2e

!===============================================================================
FUNCTION gas_pressure (self, lnd, ss) RESULT (pg)
  class(mhd_t):: self
  real, dimension(:,:,:), pointer:: lnd, ss, ee
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: pg
  !-----------------------------------------------------------------------------
  if (trim(self%eos)=='ideal') then
    self%d => self%mem(:,:,:,self%idx%d,self%it,1)
    if (self%gamma==1.0) then
      pg = self%d*self%csound**2
    else
      self%s => self%mem(:,:,:,self%idx%s,self%it,1)
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

!===============================================================================
!> temperature from ideal gas law, taking into account that the temperature is
!> normalised by k_B and mu (so that specific gas constant = 1)
!===============================================================================
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

!===============================================================================
SUBROUTINE log_density (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = log(self%d)
END SUBROUTINE

!===============================================================================
SUBROUTINE log_pressure (self, lnd, ss, v)
  class(mhd_t):: self
  real, dimension(:,:,:), pointer:: lnd, ss, v
  !-----------------------------------------------------------------------------
  v = log(self%gas_pressure(lnd,ss))
END SUBROUTINE

!===============================================================================
SUBROUTINE velocity_magnitude (self, v)
  class (mhd_t):: self
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = norm(self%p/exp(down(log(self%d))))
END SUBROUTINE

!===============================================================================
SUBROUTINE magnetic_field_magnitude (self, v)
  class (mhd_t):: self
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = norm(self%B)
END SUBROUTINE

!===============================================================================
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
