!===============================================================================
!===============================================================================
MODULE mhd_mod
  USE io_mod
  USE bits_mod
  USE omp_mod
  USE math_mod
  USE omp_timer_mod
  USE mpi_mod
  USE trace_mod
  USE kinds_mod
  USE units_mod
  USE scaling_mod
  USE scalar_mod
  USE vector_mod
  USE vector_ops
  USE stagger_mod
  USE extras_mod
  USE index_mod
  USE timestep_mod
  USE eos_mod
  USE non_ideal_mod
  USE initial_mod
  USE divb_clean_mod
  USE scratch_mod
  implicit none
  private
  type, extends(extras_t), public:: mhd_t
    real, dimension(:,:,:),   pointer:: d, dddt, e, dedt
    real, dimension(:,:,:,:), pointer:: p, dpdt, B, dBdt
    real, dimension(:,:,:),   pointer:: lnd, ee, pg
    real, dimension(:,:,:), allocatable:: qrad
    real:: nu(6)
    real:: cdtd
    real:: maxJ=0.0
    logical:: mhd=.false.
    logical:: first_time=.true.
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: pde
    procedure:: update
    procedure:: output
    procedure:: gas_pressure
    procedure:: divB_clean => divb_clean2
  end type
  logical, save:: unsigned=.true.
  logical, save:: do_smooth=.false.
  logical, save:: do_maxwell = .true.
  integer, save:: verbose=0
  integer, save:: divb_cleaner=2
CONTAINS

!===============================================================================
!> Read the namelist parameters, store them, determine the type of the solver
!> (HD or MHD) and set the indices.
!===============================================================================
SUBROUTINE pre_init (self)
  class(mhd_t) :: self
  type(index_t):: idx
  integer:: i, iv
  real, save:: nu(6)=[0.1,1.0,0.0,0.5,0.5,0.5], csound=1.0, cdtd=0.5, courant=0.2
  real(8), save:: gamma=1.4_8
  logical, save:: first_time=.true.
  logical, save:: mhd=.false.
  character(len=16), save:: eos
  namelist /stagger_params/ nu, courant, gamma, csound, cdtd, hardwire, eos, &
    do_smooth, unsigned, do_maxwell, verbose, divb_cleaner
  !----------------------------------------------------------------------------
  call trace%begin('mhd_t%pre_init')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    eos = self%eos
    mhd = self%mhd
    rewind (io%input)
    read (io%input, stagger_params)
    if (io%master) write (*, stagger_params)
    call non_ideal%init
    call timestep%init
    ! fool-proofing
    if (timestep%time_order <= 0) call mpi%abort("The Stagger solvers require `time_order` > 0. Abort!")
  end if
  call io%lock%unset
  self%eos = eos
  self%mhd = mhd
  self%nu = nu
  self%gamma = gamma
  self%csound = csound
  self%courant = courant
  self%cdtd = cdtd
  !-----------------------------------------------------------------------------
  ! Read IC parameters, which determine if this is HD or MHD.  The IC init is
  ! called again, from extras_t%init, when gamma is known and can be passed on.
  !-----------------------------------------------------------------------------
  self%initial%mhd = self%mhd
  call self%initial%init (self%kind, real(self%gamma))
  self%mhd = self%initial%mhd
  self%kind = 'stagger2e_pic'
  if (self%nv == 0) then
    if (self%mhd) then
      self%nv = 8
    else
      self%nv = 5
    end if
  else
    self%nv = self%nv + 5
    if (self%mhd) self%nv = self%nv + 3
  end if
  call self%idx%init (5, self%mhd)
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> The construction with a double test on first_time avoids entering a critical
!> region once the variable has been set.
!===============================================================================
SUBROUTINE init (self)
  class(mhd_t):: self
  integer :: i
  ! ----------------------------------------------------------------------------
  call trace%begin('mhd_t%init')
  !-----------------------------------------------------------------------------
  ! Allocate memory and mesh, etc
  !-----------------------------------------------------------------------------
  call self%gpatch_t%init
  self%type = 'mhd_t'
  !----------------------------------------------------------------------------
  ! Stagger offsets for (up-) face centered mass flux components
  !----------------------------------------------------------------------------
  self%mesh(1)%h(self%idx%px) = +0.5
  self%mesh(2)%h(self%idx%py) = +0.5
  self%mesh(3)%h(self%idx%pz) = +0.5
  !----------------------------------------------------------------------------
  ! Stagger offsets for (up-up) edge centered magnetic field components
  !----------------------------------------------------------------------------
  if (self%mhd) then
    self%mesh(2)%h(self%idx%bx) = +0.5
    self%mesh(3)%h(self%idx%bx) = +0.5
    self%mesh(1)%h(self%idx%by) = +0.5
    self%mesh(3)%h(self%idx%by) = +0.5
    self%mesh(1)%h(self%idx%bz) = +0.5
    self%mesh(2)%h(self%idx%bz) = +0.5
  end if
  if (verbose>1 .and. self%track) then
    do i=1,3
      print '("h:",8f5.2)', self%mesh(i)%h
    end do
  end if
  self%unsigned(self%idx%d)  = unsigned
  self%pervolume(self%idx%e) = unsigned
 call trace%end()
END SUBROUTINE init

!===============================================================================
SUBROUTINE pde (self)
  USE scalar_mod
  class(mhd_t):: self
  !.............................................................................
  real, dimension(:,:,:),   pointer:: d, e, dddt, dedt
  real, dimension(:,:,:),   pointer:: dpdtx, dpdty, dpdtz, dBdtx, dBdty, dBdtz
  real, dimension(:,:,:),   pointer:: Bx, By, Bz, px, py, pz, phi, qr
  real, dimension(:,:,:),   pointer:: Ux, Uy, Uz, pg
  real, dimension(:,:,:,:), pointer:: U, p, B, dpdt, dBdt
  real, dimension(:,:,:),   pointer:: Ex, Ey, Ez, Jx, Jy, Jz
  real, dimension(:,:,:),   pointer:: fdx, fdy, fdz, ldx, ldy, ldz, ddx, ddy, ddz
  real, dimension(:,:,:),   pointer:: fxx, fxy
  real, dimension(:,:,:),   pointer::      fyy, fyz
  real, dimension(:,:,:),   pointer:: fzx,      fzz
  real, dimension(:,:,:),   allocatable:: du, ss, cs, pa, pp, pb, Q, uu, divU
  real, dimension(:,:,:),   allocatable, target:: ee, lnd
  real, dimension(:,:,:),   allocatable:: Txx, Txy
  real, dimension(:,:,:),   allocatable::      Tyy, Tyz
  real, dimension(:,:,:),   allocatable:: Tzx,      Tzz
  real, dimension(:,:,:,:), allocatable, target:: emf, J
  real, dimension(:,:,:,:), allocatable, target:: fx, fy, fz, fd, ld, dd, fm
  integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
  logical, save:: first_time=.true.
  !.............................................................................
  real(8):: ds(3), dsmax
  real:: u_max
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%pde',itimer=itimer)
  self%d    => self%mem(:,:,:,self%idx%d,self%it,1)
  self%dddt => self%mem(:,:,:,self%idx%d,self%it,2)
  self%e    => self%mem(:,:,:,self%idx%e,self%it,1)
  self%dedt => self%mem(:,:,:,self%idx%e,self%it,2)
  self%p    => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  self%dpdt => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,2)
  d=>self%d; dddt=>self%dddt
  e=>self%e; dedt=>self%dedt
  p=>self%p; dpdt=>self%dpdt
  px=>p(:,:,:,1); py=>p(:,:,:,2); pz=>p(:,:,:,3);
  dpdtx=>dpdt(:,:,:,1); dpdty=>dpdt(:,:,:,2); dpdtz=>dpdt(:,:,:,3);
  if (self%mhd) then
    self%B    => self%mem(:,:,:,self%idx%bx:self%idx%bz,self%it,1)
    self%dBdt => self%mem(:,:,:,self%idx%bx:self%idx%bz,self%it,2)
    B=>self%B; dBdt=>self%dBdt
    Bx   =>   B(:,:,:,1); By   =>   B(:,:,:,2); Bz   =>   B(:,:,:,3);
    dBdtx=>dBdt(:,:,:,1); dBdty=>dBdt(:,:,:,2); dBdtz=>dBdt(:,:,:,3);
  end if
  call allocate_vectors_a (self%gn, ld, dd, emf, J, fd, fx, fy, fz, fm)
  call allocate_scalars_a (self%gn, lnd, pa, du, ee, cs, pp, pb, Q, uu, divU)
  if (self%do_pebbles) then
    U => self%vgas(:,:,:,:,self%it)
    pg => self%pressure(:,:,:,self%it)
    self%density => self%mem(:,:,:,1,:,1)
  else
    call allocate_scalars (self%gn, pg)
    call allocate_vectors (self%gn, U)
  end if
  Ux=>U(:,:,:,1); Uy=>U(:,:,:,2); Uz=>U(:,:,:,3)
  Ex=>emf(:,:,:,1); Ey=>emf(:,:,:,2); Ez=>emf(:,:,:,3)
  Jx=>J(:,:,:,1); Jy=>J(:,:,:,2); Jz=>J(:,:,:,3)
  fxx=>fx(:,:,:,1); fxy=>fx(:,:,:,2)
                    fyy=>fy(:,:,:,2); fyz=>fy(:,:,:,3)
  fzx=>fz(:,:,:,1);                   fzz=>fz(:,:,:,3)
  fdx=>fd(:,:,:,1); fdy=>fd(:,:,:,2); fdz=>fd(:,:,:,3)
  ldx=>ld(:,:,:,1); ldy=>ld(:,:,:,2); ldz=>ld(:,:,:,3)
  ddx=>dd(:,:,:,1); ddy=>dd(:,:,:,2); ddz=>dd(:,:,:,3)
  !-----------------------------------------------------------------------------
  ! Compute velocity U=momentum/density, leaving up shifted density and log
  ! density for reuse below.  U is valid in (2:gn)
  !-----------------------------------------------------------------------------
  lnd = log(d)
  ld = up(lnd)
  dd = exp(ld)
  U = p/dd
  !-----------------------------------------------------------------------------
  ! Gas pressure
  !-----------------------------------------------------------------------------
  ee = e/d                                                                      ! ops: 0 0 1
  self%lnd => lnd
  self%ee  => ee
  pg = self%gas_pressure ()
  self%pg => pg
  !-----------------------------------------------------------------------------
  ! Conservation of momentum -- first the diagonal part of the stress tensor
  ! The diagonal fxx is valid in (2:gn-1)
  !-----------------------------------------------------------------------------
  call allocate_scalars_a (self%gn, Txx, Tyy, Tzz)
  ds = self%ds
  Txx = ddxdn(ds,Ux)
  Tyy = ddydn(ds,Uy)
  Tzz = ddzdn(ds,Uz)
  divU = Txx+Tyy+Tzz                                    ! convergence rate      ! ops: 0 2 0
  !-----------------------------------------------------------------------------
  ! Divergence of velocity, and artificial pressure.  The divergence is valid
  ! in (2,gn-1)
  !-----------------------------------------------------------------------------
  dsmax = maxval(ds,self%n > 1)
  du = -dsmax*min(divU,0.0)                             ! positive part
  if (do_smooth) then
    pa = self%nu(2)*d*xyzsm(du**2)
  else
    pa = self%nu(2)*d*du**2                                                     ! ops: 0 3 0
  end if
  !-----------------------------------------------------------------------------
  ! Diagonal Reynolds stress
  !-----------------------------------------------------------------------------
  fxx = xdn(Ux)**2
  fyy = ydn(Uy)**2
  fzz = zdn(Uz)**2
  uu = sqrt(fxx+fyy+fzz)
  self%u_max = 0.0
  if (self%mhd) then
    ! Diagonal (negative) Maxwell stress (borrowing emf for scratch)
    Ex = ydn(zdn(Bx))**2
    Ey = zdn(xdn(By))**2
    Ez = xdn(ydn(Bz))**2
    ! Magnetic pressure
    pb = 0.5*(Ex+Ey+Ez)                                                         ! mhd: 2 1 0
    ! Fast mode speed + velocity
    cs = sqrt((self%gamma*pg+2.0*pb + pa)/d)                                    ! mhd: 1 1 0
    if (verbose>3) then
      u_max = self%fmaxval(cs,outer=.true.)
      self%u_max = max(self%u_max,u_max)
      print 1, 'u_max(ca)  ', u_max, self%u_max
    1 format(1x,a,1p,2g12.3)
    end if
    ! Isotropic pressure part
    pp = pg + pa + pb                                                           ! mhd: 1 0 0
    ! Total diagonal stress
    du = (dsmax*self%nu(1))*(cs + 2.*uu)                                        ! mhd: 0 2 0
    !either use Maxwell stress tensor or Lorentz force
    if (do_maxwell) then
      fxx = pp + d*(fxx - du*Txx) - Ex                                          ! mhd: 1 0 0
      fyy = pp + d*(fyy - du*Tyy) - Ey                                          ! mhd: 1 0 0
      fzz = pp + d*(fzz - du*Tzz) - Ez                                          ! mhd: 1 0 0
    else
      fxx = pp + d*(fxx - du*Txx)                                               ! mhd: 1 0 0
      fyy = pp + d*(fyy - du*Tyy)                                               ! mhd: 1 0 0
      fzz = pp + d*(fzz - du*Tzz)                                               ! mhd: 1 0 0
    end if
  else
    ! Sound speed + velocity
    cs = sqrt((self%gamma*pg+pa)/d)                                             ! ops: 1 1 1 1 0
    if (verbose>3) then
      u_max = self%fmaxval(cs,outer=.true.)
      self%u_max = max(self%u_max,u_max)
      print 1, 'u_max(cs)  ', u_max, self%u_max
    end if
    pp = pg + pa                                                                ! ops: 1 0 0
    du = (dsmax*self%nu(1))*(cs + 2.*uu)                                        ! ops: 0 2 0
    !---------------------------------------------------------------------------
    ! Total diagonal stress, with density factor
    ! Incoming fxx is valid in (2:gn-1), as is du and ddxup(Ux)
    !---------------------------------------------------------------------------
    fxx = pp + d*(fxx - du*Txx)                                                 ! ops: 2 2 0
    fyy = pp + d*(fyy - du*Tyy)                                                 ! ops: 2 2 0
    fzz = pp + d*(fzz - du*Tzz)                                                 ! ops: 2 2 0
  end if
  u_max = self%fmaxval(cs+uu,outer=.true.)
  self%u_max = max(self%u_max,u_max)
  if (verbose>2) print 1, 'u_max(cs+u)', u_max, self%u_max
  !-----------------------------------------------------------------------------
  ! du has units area per unit time, which gives velocity squared times mass
  ! density per unit time = energy per unit volume and time
  !-----------------------------------------------------------------------------
  Q = d*du*(Txx**2 + Tyy**2 + Tzz**2)                                           ! ops: 2 5 0
  !if (verbose==1) print *, 'average(Q):', self%faver(Q), &
    !self%faver(e-d**self%gamma/(self%gamma-1.0)), self%faver(fm(:,:,:,1))
  !-----------------------------------------------------------------------------
  ! Conservation of mass
  !-----------------------------------------------------------------------------
  fd = -(self%nu(4)*dd*up(du))*dup(ds,lnd)              ! diffusive mass flux   ! ops: 0 9 0
  fm = p+fd                                             ! mass flux
  dddt = -div_dn(ds,fm)                                 ! density derivative
  u_max = self%cdtd*dsmax*self%fmaxval(abs(dddt/d),outer=.true.)                ! ops: 0 0 1
  self%u_max = max(self%u_max,u_max)
  if (verbose>2) print 1, 'u_max(dddt)', u_max, self%u_max
  !-----------------------------------------------------------------------------
  ! Off-diagonal part of Reynolds stress, for now without the density factor.
  ! Note that, because of the symmmetry, we define only three components
  !-----------------------------------------------------------------------------
  fxy = yUP(Ux)*xUP(Uy)
  fyz = zUP(Uy)*yUP(Uz)
  fzx = xUP(Uz)*zUP(Ux)
  if (self%id==io%id_debug) print *, 'DPDTZ1', fzx(19:20,19,5)
  !-----------------------------------------------------------------------------
  ! Add viscous part -- edge-centered
  !-----------------------------------------------------------------------------
  call allocate_scalars_a (self%gn, Txy, Tyz, Tzx)
  Txy = ddxUP(ds,Uy)+ddyUP(ds,Ux)                                               ! ops: 0 1 0
  Tyz = ddyUP(ds,Uz)+ddzUP(ds,Uy)                                               ! ops: 0 1 0
  Tzx = ddzUP(ds,Ux)+ddxUP(ds,Uz)                                               ! ops: 0 1 0
  fxy = fxy - self%nu(3)*0.5*xUP1(yUP1(du))*Txy
  fyz = fyz - self%nu(3)*0.5*yUP1(zUP1(du))*Tyz
  fzx = fzx - self%nu(3)*0.5*zUP1(xUP1(du))*Tzx
  Q = Q + self%nu(3)*0.5*d*du*(xdn(ydn(Txy))**2 + &                             ! ops: 2 5 0
                               ydn(zdn(Tyz))**2 + &                             ! ops: 1 1 0
                               zdn(xdn(Tzx))**2)                                ! ops: 1 0 0
  call deallocate_scalars_a (Txy, Tyz, Tzx)
  !-----------------------------------------------------------------------------
  ! Multiply with edge-centered density
  !-----------------------------------------------------------------------------
  fxy = exp(yup(ldx))*fxy                                                       ! ops: 0 1 0 0 1
  fyz = exp(zup(ldy))*fyz                                                       ! ops: 0 1 0 0 1
  fzx = exp(xup(ldz))*fzx                                                       ! ops: 0 1 0 0 1
  !-----------------------------------------------------------------------------
  ! Off-diagonal Maxwell stress
  !----------------------------------------------------------------------------
  if (self%mhd.and.do_maxwell) then
    fxy = fxy - xUP(zDN(Bx))*yUP(zDN(By))                                       ! mhd: 1 1 0
    fyz = fyz - yUP(xDN(By))*zUP(xDN(Bz))                                       ! mhd: 1 1 0
    fzx = fzx - zUP(yDN(Bz))*xUP(yDN(Bx))                                       ! mhd: 1 1 0
  end if
  !-----------------------------------------------------------------------------
  ! Consistent addition of mass flux related momentum flux, symmetrized. The
  ! diagonal part is analogous, a centered version of \nu \rho U_i \grad\ln\rho,
  ! while the off-diagonal part needs to be properly symmetrized.
  !-----------------------------------------------------------------------------
  if (self%nu(4) > 0.0) then
    pa = d*du*self%nu(4)                                                        ! ops: 0 2 0
    !---------------------------------------------------------------------------
    ! fxx, pa, ddxup(ldx), and xup(Ux) are all valid in (2:gn-1)
    !---------------------------------------------------------------------------
    fxx = fxx - pa*ddxDN(ds,ldx)*xDN(Ux)                                        ! ops: 1 2 0
    fyy = fyy - pa*ddyDN(ds,ldy)*yDN(Uy)                                        ! ops: 1 2 0
    fzz = fzz - pa*ddzDN(ds,ldz)*zDN(Uz)                                        ! ops: 1 2 0
    fxy = fxy + (xUP1(fdy)*yUP(Ux)+yUP1(fdx)*xUP(Uy))                           ! ops: 2 2 0
    fyz = fyz + (yUP1(fdz)*zUP(Uy)+zUP1(fdy)*yUP(Uz))                           ! ops: 2 2 0
    fzx = fzx + (zUP1(fdx)*xUP(Uz)+xUP1(fdz)*zUP(Ux))                           ! ops: 2 2 0
    if (self%id==io%id_debug) print *, 'DPDTZ4', fzx(19:20,19,5)
  end if
  !-----------------------------------------------------------------------------
  ! Momentum equations.  The ddxdn(fxx) part is valid in (3:gn-1), and hence px
  ! only needs guard zones in (1:2) and (gn:gn)
  !-----------------------------------------------------------------------------
  dpdtx = - ddxUP1(ds,fxx) - ddyDN1(ds,fxy) - ddzDN1(ds,fzx)                   ! ops: 2 0 0
  dpdty = - ddxDN1(ds,fxy) - ddyUP1(ds,fyy) - ddzDN1(ds,fyz)                   ! ops: 2 0 0
  dpdtz = - ddxDN1(ds,fzx) - ddyDN1(ds,fyz) - ddzUP1(ds,fzz)                   ! ops: 2 0 0
  !-----------------------------------------------------------------------------
  ! External force; e.g. point source gravity.  Note that the density factor is
  ! applied here, and should not be included in the function definition
  !-----------------------------------------------------------------------------
  if (associated(scratch%force_per_unit_mass)) then
    dpdt = dpdt + dd*scratch%force_per_unit_mass                                  ! ops: 2 2 0
  end if
  if (associated(scratch%force_per_unit_volume)) then
    dpdt = dpdt + scratch%force_per_unit_volume
  end if
  !-----------------------------------------------------------------------------
  ! Conservation of energy.  The formulation that works best so far is to make
  ! the diffusive energy flux directly proportional to the diffusive mass flux,
  ! plus a contribution proportional to the gradient of specific energy ee.
  !-----------------------------------------------------------------------------
  if (self%gamma == 1d0) then                                                   ! unless isothermal
    dedt = dddt*ee                                                              ! ops: 0 1 0
  else
    block
    real:: afm
    afm = self%faver(fm(:,:,:,1))
    fm = fm*up(ee) - (self%nu(5)*dd*up(du))*dup(ds,ee)                        ! ops: 1 4 0
    dedt = -pg*divU - div_DN(ds,fm) + Q
    if (verbose==1) print '(a,10f12.6)', 'average(Q):', &
      self%time, self%time*self%faver(Q), &
      self%faver(e-d**self%gamma/(self%gamma-1.0)), &
      self%faver(- pg*divU), &
      afm, self%faver(fm(:,:,:,1))
    end block
  end if
  call deallocate_scalars_a (Txx, Tyy, Tzz)
  !-----------------------------------------------------------------------------
  ! Induction equation
  !-----------------------------------------------------------------------------
  if (self%mhd) then
    !---------------------------------------------------------------------------
    ! Electric current, resistive electric field and magnetic dissipation
    ! emf = {Ex,Ey,Ez} is the vector with components Ex, Ey, Ez
    !---------------------------------------------------------------------------
    J = curl_dn(ds,B)
    self%maxJ=maxval(abs(J))
    emf = self%nu(6)*up(du)*J                                                 ! mhd: 0 6 0  
    Q = xDN(Jx*Ex) + yDN(Jy*Ey) + zDN(Jz*Ez)                                    ! mhd: 2 3 0
    if (self%gamma /= 1.0) dedt = dedt + Q                                      ! mhd: 1 0 0
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Lorentz force
    ! dpdtz(i,j,k+½) = zup(xdn(jx(i+½,j,k)))*xdn(By(i+½,j,k+½)) - jy(i,j+½,k)*Bx(i,j+½,k+½)
    !---------------------------------------------------------------------------
    if (.not.do_maxwell) then
    !borrowing fxyz for scratch to store cell-centered J values
      ! jx -> fxx, jy -> fyy, jz -> fzz
      fxx = xdn(jx)
      fyy = ydn(jy)
      fzz = zdn(jz)
      dpdtx = dpdtx + xup(fyy)*(ydn(Bz)) - xup(fzz)*(zdn(By))
      dpdty = dpdty + yup(fzz)*(zdn(Bx)) - yup(fxx)*(xdn(Bz))
      dpdtz = dpdtz + zup(fxx)*(xdn(By)) - zup(fyy)*(ydn(Bx))
    end if
    !---------------------------------------------------------------------------
    ! Advective electric field on cell FACES
    ! Ex(i+½,j,k) = xup(zdn(Uz(i,j+½,k)))*zdn(By(i+½,j,k+½)) - ...
    !---------------------------------------------------------------------------
    !Ex = Ex + xUP(yDN(Uz))*zDN(By) - xUP(zDN(Uy))*yDN(Bz)                       ! mhd: 2 2 0
    !Ey = Ey + yUP(zDN(Ux))*xDN(Bz) - yUP(xDN(Uz))*zDN(Bx)                       ! mhd: 2 2 0
    !Ez = Ez + zUP(xDN(Uy))*yDN(Bx) - zUP(yDN(Ux))*xDN(By)                       ! mhd: 2 2 0
    Ex = Ex + xUP(zDN(Uz))*zDN(By) - xUP(yDN(Uy))*yDN(Bz)                       ! mhd: 2 2 0
    Ey = Ey + yUP(xDN(Ux))*xDN(Bz) - yUP(zDN(Uz))*zDN(Bx)                       ! mhd: 2 2 0
    Ez = Ez + zUP(yDN(Uy))*yDN(Bx) - zUP(xDN(Ux))*xDN(By)                       ! mhd: 2 2 0
    if (non_ideal%is_used) call io%abort("non ideal MHD has to be checked for staggering for PIC")
    call non_ideal%update (self%gn, minval(ds), d, Jx, Jy, Jz, Bx, By, Bz, Q, emf, self%u_max)
    dBdt = -curl_up(ds,emf)
  end if
  !-----------------------------------------------------------------------------
  ! External heating / cooling, if any; e.g. radiative or conductive
  !-----------------------------------------------------------------------------
  if (associated(scratch%heating_per_unit_mass)) then
    dedt = dedt + d*scratch%heating_per_unit_mass                                  ! ops: 2 2 0
  end if
  if (associated(scratch%heating_per_unit_volume)) then
    dedt = dedt + scratch%heating_per_unit_volume
  end if
  !-----------------------------------------------------------------------------
  ! Time step limitation for rate of change dE/dt = (de/dt-(e/d)*dddt)/d
  !-----------------------------------------------------------------------------
  if (self%gamma /= 1d0) then
    u_max = self%cdtd*dsmax*self%fmaxval(abs((dedt-ee*dddt)/d),outer=.true.)    ! ops: 1 1 1
    self%u_max = max(self%u_max,u_max)
    if (verbose>2) print 1, 'u_max(dedt)', u_max, self%u_max
  end if
  !-----------------------------------------------------------------------------
  call deallocate_vectors_a (ld, dd, emf, J, fd, fx, fy, fz, fm)
  call deallocate_scalars_a (lnd, pa, du, ee, cs, pp, pb, Q, uu, divU)
  if (.not.self%do_pebbles) then
    call deallocate_vectors (U)
    call deallocate_scalars (pg)
  end if
  call trace%end (itimer)
  if (mpi%master .and. first_time) then
    print *, stagger%count, ' stagger calls', stagger%flops/product(real(self%n)),' flops per interior point'
    stagger%count=0
    stagger%flops=0
    first_time = .false.
  end if
END SUBROUTINE pde

!===============================================================================
!> Take a complete update step.  This can be overloaded in higher levels, e.g.
!> in experiment_mod, where boundary conditions may be inserted
!===============================================================================
SUBROUTINE update (self)
  class(mhd_t):: self
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%update',itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Output the state after guard zone downloads, selfgravity, and boundary
  ! conditions. If another call is placed at an earlier point in the task handling
  ! the call here will have no effect, since out_next has already been updated.
  !-----------------------------------------------------------------------------
  call self%output
  if (self%gamma==1d0) then
    self%mem(:,:,:,self%idx%e,self%it,1) = &
    self%mem(:,:,:,self%idx%d,self%it,1)*self%csound**2
  end if
  if (self%mhd) call self%divB_clean
  call self%pde
  call self%courant_condition
  if (verbose>2) &
    call self%check_density('before timestep')
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call self%check_density(' after timestep')
  !-----------------------------------------------------------------------------
  call self%counter_update
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Perform a single iteration of a div(B) clean.  The limit of stability is
!> eps=1./6.  This is intended to prevent build-up of layers of large div(B)
!> between patches. An exact solution would be to solve an internal Laplace
!> equation, with boundary conditiosn = the boundary value of div(B), and then
!> subtract the gradient of the resulting potential.
!>
!> The problem can also be greatly reduced by using a higher order time
!> interpoation scheme, thus making the electric field computed from the
!> guard zone values of B and U a better match to the neighbors E-fields.
!===============================================================================
SUBROUTINE divb_clean1(self)
  class(mhd_t):: self
  !.............................................................................
  real, dimension(:,:,:,:), pointer:: b
  real, dimension(:,:,:)  , pointer:: bx, by, bz, divb
  real, dimension(:,:,:,:)  , allocatable:: dndivb
  real, parameter:: eps=0.11
  !-----------------------------------------------------------------------------
  b  => self%mem(:,:,:,self%idx%bx:self%idx%bz,self%it,1)
  bx => b(:,:,:,1)
  by => b(:,:,:,2)
  bz => b(:,:,:,3)
  call allocate_scalars (self%gn, divb)
  call allocate_vectors_a (self%gn, dndivb)
  divb = ddup(xdn(ydn(zdn(Bx))),1) + &
         ddup(ydn(zdn(xdn(By))),2) + &
         ddup(zdn(xdn(ydn(Bz))),3)
  bx = bx + zup(yup(xup(eps*dddn(divb,1))))
  by = by + xup(zup(yup(eps*dddn(divb,2))))
  bz = bz + yup(xup(zup(eps*dddn(divb,3))))
  call deallocate_vectors_a (dndivb)
  call deallocate_scalars(divb)
    return
contains
  !-----------------------------------------------------------------------------
  function ddup(f,i) result (g)
  real:: f(:,:,:), g(size(f,1),size(f,2),size(f,3))
  integer:: i
  !.............................................................................
  g = cshift(f,1,i) - f
  end function
  !-----------------------------------------------------------------------------
  function dddn(f,i) result (g)
  real:: f(:,:,:), g(size(f,1),size(f,2),size(f,3))
  integer:: i
  !.............................................................................
  g = f - cshift(f,-1,i)
  end function
END SUBROUTINE divb_clean1

!===============================================================================
!> Enforce div(B)=0 in the guard zones, by succesively imposing the condition,
!> layer-by-layer, on the component perpendicular to each of the patch faces.
!===============================================================================
SUBROUTINE divb_clean2 (self)
  class(mhd_t):: self
  !.............................................................................
  real, dimension(:,:,:), pointer:: bx, by, bz
  real, dimension(:,:,:), allocatable:: divb
  real, parameter:: eps=0.11
  integer:: l(3), u(3), ix, iy, iz
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%divb_clean2', itimer=itimer)
  bx => self%mem(:,:,:,self%idx%bx,self%it,1)
  by => self%mem(:,:,:,self%idx%by,self%it,1)
  bz => self%mem(:,:,:,self%idx%bz,self%it,1)
  l = self%mesh%lb
  u = self%mesh%ub-1
  !--- x-direction -------------------------------------------------------------
  do ix=self%mesh(1)%lo,self%mesh(1)%lb,-1
    bx(ix  ,l(2):u(2),l(3):u(3)) = &
    bx(ix+1,l(2):u(2),l(3):u(3)) + ( &
      by(ix,l(2)+1:u(2)+1,l(3)  :u(3)  ) - by(ix,l(2):u(2),l(3):u(3)) + &
      bz(ix,l(2)  :u(2)  ,l(3)+1:u(3)+1) - bz(ix,l(2):u(2),l(3):u(3)))
  end do
  do ix=self%mesh(1)%ui,self%mesh(1)%ub-1
    bx(ix+1,l(2):u(2),l(3):u(3)) = &
    bx(ix  ,l(2):u(2),l(3):u(3)) - ( &
      by(ix,l(2)+1:u(2)+1,l(3)  :u(3)  ) - by(ix,l(2):u(2),l(3):u(3)) + &
      bz(ix,l(2)  :u(2)  ,l(3)+1:u(3)+1) - bz(ix,l(2):u(2),l(3):u(3)))
  end do
  !--- y-direction -------------------------------------------------------------
  do iy=self%mesh(2)%lo,self%mesh(2)%lb,-1
    by(l(1):u(1),iy  ,l(3):u(3)) = &
    by(l(1):u(1),iy+1,l(3):u(3)) + ( &
      bx(l(1)+1:u(1)+1,iy,l(3)  :u(3)  ) - bx(l(1):u(1),iy,l(3):u(3)) + &
      bz(l(1)  :u(1)  ,iy,l(3)+1:u(3)+1) - bz(l(1):u(1),iy,l(3):u(3)))
  end do
  do iy=self%mesh(2)%ui,self%mesh(2)%ub-1
    by(l(1):u(1),iy+1,l(3):u(3)) = &
    by(l(1):u(1),iy  ,l(3):u(3)) - ( &
      bx(l(1)+1:u(1)+1,iy,l(3)  :u(3)  ) - bx(l(1):u(1),iy,l(3):u(3)) + &
      bz(l(1)  :u(1)  ,iy,l(3)+1:u(3)+1) - bz(l(1):u(1),iy,l(3):u(3)))
  end do
  !--- z-direction -------------------------------------------------------------
  do iz=self%mesh(3)%lo,self%mesh(3)%lb,-1
    bz(l(1):u(1),l(2):u(2),iz  ) = &
    bz(l(1):u(1),l(2):u(2),iz+1) + ( &
      bx(l(1)+1:u(1)+1,l(2)  :u(2)  ,iz) - bx(l(1):u(1),l(2):u(2),iz) + &
      by(l(1)  :u(1)  ,l(2)+1:u(2)+1,iz) - by(l(1):u(1),l(2):u(2),iz))
  end do
  do iz=self%mesh(3)%ui,self%mesh(3)%ub-1
    bz(l(1):u(1),l(2):u(2),iz+1) = &
    bz(l(1):u(1),l(2):u(2),iz  ) - ( &
      bx(l(1)+1:u(1)+1,l(2)  :u(2)  ,iz) - bx(l(1):u(1),l(2):u(2),iz) + &
      by(l(1)  :u(1)  ,l(2)+1:u(2)+1,iz) - by(l(1):u(1),l(2):u(2),iz))
  end do
  call trace%end (itimer)
END SUBROUTINE divb_clean2

!===============================================================================
!> Add extra, MHD-related information to the patch .txt output file
!===============================================================================
SUBROUTINE output (self)
  class(mhd_t):: self
  character(len=64):: filename
  real(8):: out_next
  logical:: exists
  !-----------------------------------------------------------------------------
  call trace%begin ('mhd_t%output')
  out_next = self%out_next
  call self%extras_t%output
  if (self%time >= out_next .and. io%do_output .and. io%do_legacy) then
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
  end if
  call trace%end()
END SUBROUTINE output

!===============================================================================
FUNCTION gas_pressure (self) RESULT (pg)
  class(mhd_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: pg
  real, dimension(:,:,:), pointer:: d, e, tt
  integer, save:: nprint=1
  !-----------------------------------------------------------------------------
  call trace%begin ('mhd_t%gas_pressure')
  if (trim(self%eos)=='ideal') then
    d => self%mem(:,:,:,self%idx%d,self%it,1)
    e => self%mem(:,:,:,self%idx%e,self%it,1)
    if (self%gamma==1d0) then
      pg = d*self%csound**2
    else
      pg = e*(self%gamma-1)
    end if
  else
    tt => self%mem(:,:,:,self%idx%tt,self%it,1)
    call eos%lookup_table (shape(self%lnd), lnd=self%lnd, ee=self%ee, pg=pg, tt=tt)
    if (io%master .and. (verbose>1 .or. nprint>0)) then
      nprint = nprint-1
      write(io%output,*) self%id, 'lookup_table: lnd =', minval(self%lnd), maxval(self%lnd)
      write(io%output,*) self%id, 'lookup_table:  ee =', minval(self%ee ), maxval(self%ee )
      write(io%output,*) self%id, 'lookup_table:  pg =', minval(pg),       maxval(pg)
      write(io%output,*) self%id, 'lookup_table:  tt =', minval(tt),       maxval(tt)
    end if
  end if
  call trace%end()
END FUNCTION gas_pressure

END MODULE mhd_mod
