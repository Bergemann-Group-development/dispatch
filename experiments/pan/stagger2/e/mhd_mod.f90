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
  USE force_mod
  USE index_mod
  USE timestep_mod
  USE eos_mod
  USE non_ideal_mod
  USE initial_mod
  USE scratch_mod
  USE eos_mod
  implicit none
  private
  type, extends(extras_t), public:: mhd_t
    real, dimension(:,:,:),   pointer, contiguous:: d, dddt, e, dedt
    real, dimension(:,:,:,:), pointer:: p, dpdt, B, dBdt
    real, dimension(:,:,:),   pointer, contiguous:: lnd, ee, pg
    real, dimension(:,:,:),   pointer, contiguous:: qrad
    procedure(divB_clean1),   pointer    :: divB_clean => null()
    real:: nu(6)
    real:: cdtd
    logical:: mhd=.false.
    logical:: first_time=.true.
    type(force_t):: force
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: pde
    procedure:: update
    procedure:: output
    procedure:: gas_pressure
  end type
  real, save:: d_floor=1e-4
  logical, save:: test_p=.false.
  logical, save:: unsigned=.false.
  logical, save:: do_smooth=.false.
  logical, save:: do_temperature = .false.
  integer, save:: verbose=0
  integer, save:: divb_cleaner=2
CONTAINS

!===============================================================================
!> Read the namelist parameters, store them, determine the type of the solver
!> (HD or MHD) and set the indices.
!===============================================================================
SUBROUTINE pre_init (self)
  class(mhd_t), target:: self
  type(index_t):: idx
  integer:: i, iv, iostat
  real, save:: nu(6)=[0.1,1.0,0.0,0.5,0.5,0.5], csound=1.0, cdtd=0.5, courant=0.2
  logical, save:: first_time=.true.
  logical, save:: mhd=.false.
  character(len=16), save:: eos
  namelist /stagger_params/ nu, courant, csound, cdtd, hardwire, eos, &
    do_smooth, unsigned, do_temperature, verbose, divb_cleaner
  !----------------------------------------------------------------------------
  call trace%begin ('mhd_t%pre_init')
  !-----------------------------------------------------------------------------
  ! Read IC parameters, which determine if this is HD or MHD.  The IC init is
  ! called again, from extras_t%init, when gamma is known and can be passed on.
  !-----------------------------------------------------------------------------
  self%initial%mhd = self%mhd
  call self%initial%init (self%kind)
  self%mhd = self%initial%mhd
  if (self%mhd) then
    self%kind = 'stagger2e_mhd_patch'
    self%nv = 8
  else
    self%kind = 'stagger2e_hd_patch'
    self%nv = 5
  end if
  self%nw = 2
  !-----------------------------------------------------------------------------
  ! Call extras pre-init
  !-----------------------------------------------------------------------------
  call self%extras_t%pre_init
  call io%lock%set ('mhd_t%pre_init')
  if (first_time) then
    first_time = .false.
    eos = self%eos
    mhd = self%mhd
    rewind (io%input)
    read (io%input, stagger_params, iostat=iostat)
    call io%namelist_warning ('stagger_params', iostat)
    write (io%namelist, stagger_params)
    call non_ideal%init
    call timestep%init
    ! fool-proofing
    if (timestep%time_order <= 0) call mpi%abort("The Stagger solvers require `time_order` > 0. Abort!")
  end if
  call io%lock%unset ('mhd_t%pre_init')
  self%eos = eos
  self%mhd = mhd
  self%nu = nu
  self%csound = csound
  self%courant = courant
  self%cdtd = cdtd
  !-----------------------------------------------------------------------------
  ! Choose div(B) cleaner
  !-----------------------------------------------------------------------------
  select case(divb_cleaner)
  case(1)
    self%divb_clean => divb_clean1
  case(2)
    self%divb_clean => divb_clean2
  case default
    call io%abort('mhd_mod:: divb cleaner method not understood')
  end select
  call self%idx%init (ienergy=5, mhd=self%mhd)
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
  call self%extras_t%init
  self%type = 'mhd_t'
  self%unsigned(self%idx%d)  = unsigned
  self%pervolume(self%idx%e) = unsigned
  do i=1,3
    self%mesh(i)%h(self%idx%px+i-1) = -0.5
    if (self%mhd) self%mesh(i)%h(self%idx%bx+i-1) = -0.5
  end do
  if (verbose > 1 .and. self%track) then
    do i=1,3
      print '("h:",8f5.2)', self%mesh(i)%h
    end do
  end if
  call trace%end()
END SUBROUTINE init

!===============================================================================
SUBROUTINE pde (self)
  USE scalar_mod
  class(mhd_t):: self
  !.............................................................................
  real, dimension(:,:,:),   pointer, contiguous:: d, e, dddt, dedt
  real, dimension(:,:,:),   pointer:: dpdtx, dpdty, dpdtz, dBdtx, dBdty, dBdtz
  real, dimension(:,:,:),   pointer:: Bx, By, Bz, px, py, pz, phi
  real, dimension(:,:,:),   pointer:: Ux, Uy, Uz
  real, dimension(:,:,:),   pointer, contiguous:: pg
  real, dimension(:,:,:,:), pointer:: U, p, B, dpdt, dBdt
  real, dimension(:,:,:),   pointer:: Ex, Ey, Ez, Jx, Jy, Jz
  real, dimension(:,:,:),   pointer:: fdx, fdy, fdz, ldx, ldy, ldz, ddx, ddy, ddz
  real, dimension(:,:,:),   pointer:: fxx, fxy
  real, dimension(:,:,:),   pointer::      fyy, fyz
  real, dimension(:,:,:),   pointer:: fzx,      fzz
  real, dimension(:,:,:),   allocatable:: du, ss, cs, pa, pp, pb, Q, uu, divU
  real, dimension(:,:,:),   allocatable, target:: ee, lnd
  real, dimension(:,:,:),   allocatable:: dxy, dyz, dzx
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
  associate (i => self%idx)
  self%d => self%mem(:,:,:,i%d      ,self%new,1); self%dddt => self%mem(:,:,:,i%d      ,self%new,2)
  self%e => self%mem(:,:,:,i%e      ,self%new,1); self%dedt => self%mem(:,:,:,i%e      ,self%new,2)
  self%p => self%mem(:,:,:,i%px:i%pz,self%new,1); self%dpdt => self%mem(:,:,:,i%px:i%pz,self%new,2)
  d=>self%d; dddt=>self%dddt; e=>self%e; dedt=>self%dedt;
  p=>self%p; dpdt=>self%dpdt
  px=>p(:,:,:,1); py=>p(:,:,:,2); pz=>p(:,:,:,3);
  dpdtx=>dpdt(:,:,:,1); dpdty=>dpdt(:,:,:,2); dpdtz=>dpdt(:,:,:,3);
  if (self%mhd) then
    self%B => self%mem(:,:,:,i%bx:i%bz,self%new,1); self%dBdt => self%mem(:,:,:,i%bx:i%bz,self%new,2)
    B=>self%B; dBdt=>self%dBdt
    Bx=>B(:,:,:,1); By=>B(:,:,:,2); Bz=>B(:,:,:,3);
    dBdtx=>dBdt(:,:,:,1); dBdty=>dBdt(:,:,:,2); dBdtz=>dBdt(:,:,:,3);
  end if
  end associate
  call allocate_vectors_a (self%gn, ld, dd, emf, J, fd, fx, fy, fz, fm)
  call allocate_scalars_a (self%gn, lnd, pa, du, ee, cs, pp, pb, Q, uu, divU)
  if (self%do_pebbles) then
    U => self%vgas(:,:,:,:,self%new)
    pg => self%pressure(:,:,:,self%new)
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
  ! Compute velocity U=momentum/density, leaving down shifted density and log
  ! density for reuse below.  U is valid in (2:gn)
  !-----------------------------------------------------------------------------
  lnd = log(d)                                                                  ! ops: 0 0 0 0 3
  ld = down(lnd)
  dd = exp(ld)                                                                  ! ops: 0 0 0 0 3
  U = p/dd                                                                      ! ops: 0 0 3
  !-----------------------------------------------------------------------------
  ! Gas pressure. The 2nd table variable is entropy per unit mass, effectively
  ! a log variable, from the table point of view.
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
  Txx = ddxup(ds,Ux)
  Tyy = ddyup(ds,Uy)
  Tzz = ddzup(ds,Uz)
  divU = Txx+Tyy+Tzz                                    ! convergence rate      ! ops: 0 2 0
  !-----------------------------------------------------------------------------
  ! Divergence of velocity, and artificial pressure.  The divergence is valid
  ! in (2,gn-1)
  !-----------------------------------------------------------------------------
  dsmax = maxval(ds,self%n > 2)
  du = -dsmax*min(divU,0.0)                             ! positive part
  if (do_smooth) then
    pa = self%nu(2)*d*xyzsm(du**2)
  else
    pa = self%nu(2)*d*du**2                                                     ! ops: 0 3 0
  end if
  !-----------------------------------------------------------------------------
  ! Diagonal Reynolds stress
  !-----------------------------------------------------------------------------
  fxx = xup(Ux)**2
  fyy = yup(Uy)**2
  fzz = zup(Uz)**2
  uu = sqrt(fxx+fyy+fzz)
  self%u_max = 0.0
  if (self%mhd) then
    ! Diagonal Reynolds stress
    ! Diagonal (negative) Maxwell stress (borrowing emf for scratch)
    Ex = xup(Bx)**2                                                             ! mhd: 0 1 0
    Ey = yup(By)**2                                                             ! mhd: 0 1 0
    Ez = zup(Bz)**2                                                             ! mhd: 0 1 0
    ! Magnetic pressure
    pb = 0.5*(Ex+Ey+Ez)                                                         ! mhd: 2 1 0
    cs = sqrt((eos%gamma*pg+2.0*pb + pa)/d)                                    ! mhd: 1 1 0
    if (verbose>3) then
      u_max = self%fmaxval(cs)
      self%u_max = max(self%u_max,u_max)
      print 1, 'u_max(ca)  ', u_max, self%u_max
    1 format(1x,a,1p,2g12.3)
    end if
    ! Isotropic pressure part
    pp = pg + pa + pb                                                           ! mhd: 1 0 0
    ! Total diagonal stress
    du = (dsmax*self%nu(1))*cs                                                  ! mhd: 0 2 0
    fxx = pp + d*(fxx - du*Txx) - Ex                                            ! mhd: 1 0 0
    fyy = pp + d*(fyy - du*Tyy) - Ey                                            ! mhd: 1 0 0
    fzz = pp + d*(fzz - du*Tzz) - Ez                                            ! mhd: 1 0 0
  else
    ! Sound speed + velocity
    cs = sqrt((eos%gamma*pg+pa)/d)                                             ! ops: 1 1 1 1 0
    if (verbose>3) then
      u_max = self%fmaxval(cs)
      self%u_max = max(self%u_max,u_max)
      print 1, 'u_max(cs)  ', u_max, self%u_max
    end if
    pp = pg + pa                                                                ! ops: 1 0 0
    du = (dsmax*self%nu(1))*cs                                                  ! ops: 0 2 0
    !---------------------------------------------------------------------------
    ! Total diagonal stress, with density factor
    ! Incoming fxx is valid in (2:gn-1), as is du and ddxup(Ux)
    !---------------------------------------------------------------------------
    fxx = pp + d*(fxx - du*Txx)                                                 ! ops: 2 2 0
    fyy = pp + d*(fyy - du*Tyy)                                                 ! ops: 2 2 0
    fzz = pp + d*(fzz - du*Tzz)                                                 ! ops: 2 2 0
  end if
  u_max = self%fmaxval(cs+uu) + self%fmaxval(cs)*self%nu(1)*3.0
  self%u_max = max(self%u_max,u_max)
  if (verbose>2) print 1, 'u_max(cs+u)', u_max, self%u_max
  !-----------------------------------------------------------------------------
  ! du has units area per unit time, which gives velocity squared times mass
  ! density per unit time = energy per unit volume and time
  !-----------------------------------------------------------------------------
  Q = d*du*(Txx**2 + Tyy**2 + Tzz**2)                                           ! ops: 2 5 0
  !if (verbose==1) print *, 'average(Q):', self%faver(Q), &
    !self%faver(e-d**eos%gamma/(eos%gamma-1.0)), self%faver(fm(:,:,:,1))
  !-----------------------------------------------------------------------------
  ! Conservation of mass
  !-----------------------------------------------------------------------------
  fd = -(self%nu(4)*dd*down(du))*ddown(ds,lnd)          ! diffusive mass flux   ! ops: 0 9 0
  fm = p+fd                                             ! mass flux             ! ops: 3 0 0
  dddt = -div(ds,fm)                                    ! density derivative
  u_max = self%cdtd*dsmax*self%fmaxval(abs(dddt/d))                             ! ops: 0 0 1
  self%u_max = max(self%u_max,u_max)
  if (verbose>2) print 1, 'u_max(dddt)', u_max, self%u_max
  !-----------------------------------------------------------------------------
  ! Test, explictly equal to the IDL re-inaction
  !-----------------------------------------------------------------------------
  call allocate_scalars_a (self%gn, Txy, Tyz, Tzx)
  if (test_p) then
    if (verbose>2) print *,self%id,' special case'
    pg=d*self%csound**2
    ux=px/exp(xdn(lnd))
    uy=py/exp(ydn(lnd))
    uz=pz/exp(zdn(lnd))
    !
    Txx=ddxup(ds,ux)
    Tyy=ddyup(ds,uy)
    Tzz=ddzup(ds,uz)
    !
    Txy=(ddxdn(ds,uy)+ddydn(ds,ux))*0.5
    Tyz=(ddydn(ds,uz)+ddzdn(ds,uy))*0.5
    Tzx=(ddzdn(ds,ux)+ddxdn(ds,uz))*0.5
    !
    call allocate_scalars_a (self%gn, dxy, dyz, dzx)
    dxy=exp(xdn(ydn(lnd)))
    dyz=exp(ydn(zdn(lnd)))
    dzx=exp(zdn(xdn(lnd)))
    !
    fxx=d  *(xup(ux)**2     -ds(1)*self%nu(1)*Txx)+pg
    fyy=d*  (yup(uy)**2     -ds(1)*self%nu(1)*Tyy)+pg
    fzz=d*  (zup(uz)**2     -ds(1)*self%nu(1)*Tzz)+pg
    fxy=dxy*(xdn(uy)*ydn(ux)-ds(1)*self%nu(1)*Txy)
    fyz=dyz*(ydn(uz)*zdn(uy)-ds(1)*self%nu(1)*Tyz)
    fzx=dzx*(zdn(ux)*xdn(uz)-ds(1)*self%nu(1)*Tzx)
    !
    call deallocate_scalars_a (dxy, dyz, dzx)
    !
    dpdtx=-ddxdn(ds,fxx)-ddyup(ds,fxy)-ddzup(ds,fzx)
    dpdty=-ddxup(ds,fxy)-ddydn(ds,fyy)-ddzup(ds,fyz)
    dpdtz=-ddxup(ds,fzx)-ddyup(ds,fyz)-ddzdn(ds,fzz)
  !-----------------------------------------------------------------------------
  ! Original part
  !-----------------------------------------------------------------------------
  else
    !-----------------------------------------------------------------------------
    ! Off-diagonal part of Reynolds stress, for now without the density factor.
    ! Note that, because of the symmmetry, we define only three components
    !-----------------------------------------------------------------------------
    fxy = ydn(Ux)*xdn(Uy)                                                         ! ops: 0 1 0
    fyz = zdn(Uy)*ydn(Uz)                                                         ! ops: 0 1 0
    fzx = xdn(Uz)*zdn(Ux)                                                         ! ops: 0 1 0
    !-----------------------------------------------------------------------------
    ! Add viscous part -- edge-centered
    !-----------------------------------------------------------------------------
    Txy = ddxdn(ds,Uy)+ddydn(ds,Ux)                                               ! ops: 0 1 0
    Tyz = ddydn(ds,Uz)+ddzdn(ds,Uy)                                               ! ops: 0 1 0
    Tzx = ddzdn(ds,Ux)+ddxdn(ds,Uz)                                               ! ops: 0 1 0
    fxy = fxy - self%nu(3)*0.5*xdn1(ydn1(du))*Txy                                 ! ops: 1 3 0
    fyz = fyz - self%nu(3)*0.5*ydn1(zdn1(du))*Tyz                                 ! ops: 1 3 0
    fzx = fzx - self%nu(3)*0.5*zdn1(xdn1(du))*Tzx                                 ! ops: 1 3 0
    Q = Q + self%nu(3)*0.5*d*du*(xup(yup(Txy))**2 + &                             ! ops: 2 5 0
                                 yup(zup(Tyz))**2 + &                             ! ops: 1 1 0
                                 zup(xup(Tzx))**2)                                ! ops: 1 0 0
    call deallocate_scalars_a (Txy, Tyz, Tzx)
    !-----------------------------------------------------------------------------
    ! Multiply with edge-centered density
    !-----------------------------------------------------------------------------
    fxy = exp(ydn(ldx))*fxy                                                       ! ops: 0 1 0 0 1
    fyz = exp(zdn(ldy))*fyz                                                       ! ops: 0 1 0 0 1
    fzx = exp(xdn(ldz))*fzx                                                       ! ops: 0 1 0 0 1
    !-----------------------------------------------------------------------------
    ! Off-diagonal Maxwell stress
    !-----------------------------------------------------------------------------
    if (self%mhd) then
      fxy = fxy - ydn(Bx)*xdn(By)                                                 ! mhd: 1 1 0
      fyz = fyz - zdn(By)*ydn(Bz)                                                 ! mhd: 1 1 0
      fzx = fzx - xdn(Bz)*zdn(Bx)                                                 ! mhd: 1 1 0
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
      fxx = fxx - pa*ddxup(ds,ldx)*xup(Ux)                                        ! ops: 1 2 0
      fyy = fyy - pa*ddyup(ds,ldy)*yup(Uy)                                        ! ops: 1 2 0
      fzz = fzz - pa*ddzup(ds,ldz)*zup(Uz)                                        ! ops: 1 2 0
      fxy = fxy + (xdn1(fdy)*ydn(Ux)+ydn1(fdx)*xdn(Uy))                           ! ops: 2 2 0
      fyz = fyz + (ydn1(fdz)*zdn(Uy)+zdn1(fdy)*ydn(Uz))                           ! ops: 2 2 0
      fzx = fzx + (zdn1(fdx)*xdn(Uz)+xdn1(fdz)*zdn(Ux))                           ! ops: 2 2 0
    end if
    !-----------------------------------------------------------------------------
    ! Momentum equations.  The ddxdn(fxx) part is valid in (3:gn-1), and hence px
    ! only needs guard zones in (1:2) and (gn:gn)
    !-----------------------------------------------------------------------------
    dpdtx = - ddxdn1(ds,fxx) - ddyup1(ds,fxy) - ddzup1(ds,fzx)                    ! ops: 2 0 0
    dpdty = - ddxup1(ds,fxy) - ddydn1(ds,fyy) - ddzup1(ds,fyz)                    ! ops: 2 0 0
    dpdtz = - ddxup1(ds,fzx) - ddyup1(ds,fyz) - ddzdn1(ds,fzz)                    ! ops: 2 0 0
  end if
  !-----------------------------------------------------------------------------
  ! External force
  !-----------------------------------------------------------------------------
  if (associated(self%force%selected)) then
    if (abs(self%time-self%out_next) > 4.0*self%dtime) then
      dpdt = dpdt + dd*self%force%selected (self%time, d, p, Ux, Uy, Uz, self%mesh)! ops: 1 0 0
    end if
  end if
  if (associated(scratch%force_per_unit_mass)) then
   dpdt = dpdt + dd*scratch%force_per_unit_mass                                    ! ops: 3 3 0
  end if
  if (associated(scratch%force_per_unit_volume)) then
   dpdt = dpdt + scratch%force_per_unit_volume                                     ! ops: 3 3 0
  end if
  !-----------------------------------------------------------------------------
  ! Conservation of energy; reuse the mass flux for the energy flux
  !-----------------------------------------------------------------------------
  if (eos%gamma == 1d0) then                                                   ! unless isothermal
    dedt = dddt*ee                                                              ! ops: 0 1 0
  else
    block
    real:: afm
    afm = self%faver(fm(:,:,:,1))
    fm = fm*down(ee) - (self%nu(5)*dd*down(du))*ddown(ds,ee)                    ! ops: 1 4 0
    dedt = -div(ds,fm) - pg*divU + Q
    if (verbose==1) print '(a,10f12.6)', 'average(Q):', &
      self%time, self%time*self%faver(Q), &
      self%faver(e-d**eos%gamma/(eos%gamma-1.0)), &
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
    !---------------------------------------------------------------------------
    J = curl_dn(ds,B)
    emf = self%nu(6)*edge(du)*J                                                 ! mhd: 0 6 0
    Q = yup(zup(Jx*Ex)) + zup(xup(Jy*Ey)) + xup(yup(Jz*Ez))                     ! mhd: 2 3 0
    if (eos%gamma/=1.0) &
      dedt = dedt + Q                                                           ! mhd: 1 0 0
    !---------------------------------------------------------------------------
    ! Advective electric field on cell edges
    !---------------------------------------------------------------------------
    Ex = Ex + ydn(Uz)*zdn(By) - zdn(Uy)*ydn(Bz)                                 ! mhd: 2 2 0
    Ey = Ey + zdn(Ux)*xdn(Bz) - xdn(Uz)*zdn(Bx)                                 ! mhd: 2 2 0
    Ez = Ez + xdn(Uy)*ydn(Bx) - ydn(Ux)*xdn(By)                                 ! mhd: 2 2 0
    call non_ideal%update (self%gn, minval(ds), d, Jx, Jy, Jz, Bx, By, Bz, Q, emf, self%u_max)
    dBdt = -curl_up(ds,emf)
  end if
  !-----------------------------------------------------------------------------
  ! Radiative heating / cooling, if any
  !-----------------------------------------------------------------------------
  if (self%idx%qr > 0) then
    dedt = dedt + self%mem(:,:,:,self%idx%qr,self%new,1)
  end if
  !-----------------------------------------------------------------------------
  ! Time step limitation for rate of change dE/dt = (de/dt-(e/d)*dddt)/d
  !-----------------------------------------------------------------------------
  if (eos%gamma /= 1d0) then
    u_max = self%cdtd*dsmax*self%fmaxval(abs((dedt-ee*dddt)/d)) ! ops: 1 1 1
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
  integer:: i, it
  real, dimension(:,:,:), pointer:: d
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%update',itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Output the state after guard zone downloads, selfgravity, and boundary
  ! conditions. If another call is placed at an earlier point in the task handling
  ! the call here will have no effect, since out_next has already been updated.
  !-----------------------------------------------------------------------------
!print *, self%id, self%fmaxval(self%mem(:,:,:,1,1,1)), 'mhd_t%update 1'
!print *, self%id, self%fmaxval(self%mem(:,:,:,2,1,1)), 'mhd_t%update 2'
!print *, self%id, self%fmaxval(self%mem(:,:,:,3,1,1)), 'mhd_t%update 3'
!print *, self%id, self%fmaxval(self%mem(:,:,:,4,1,1)), 'mhd_t%update 4'
  call self%output
  if (eos%gamma==1d0) then
    self%mem(:,:,:,self%idx%e,self%new,1) = &
    self%mem(:,:,:,self%idx%d,self%new,1)*self%csound**2
  end if
  call self%pde
  call self%courant_condition
  if (verbose>2) &
    call self%check_density('before timestep')
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh, self%lock)
  d => self%mem(:,:,:,self%idx%d,self%new,1)
  d = max(d,d_floor)
  if (self%id==io%id_debug) then
    do i=1,self%nt
      it = self%iit(i)
      print'(a,i7,2x,4i3,2x,1p,e14.5,a,4e12.4)', &
        'id,it,new,i,iit(i),t(iit(i)):', &
        self%id, self%new, self%new, i, it, self%t(it), 'minmax(d):', &
        self%fminval(self%mem(:,:,:,1,it,2)), self%fmaxval(self%mem(:,:,:,1,it,2)), &
        self%fminval(self%mem(:,:,:,1,it,1)), self%fmaxval(self%mem(:,:,:,1,it,1))
    end do
  end if
  call self%check_density(' after timestep')
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
  divb = ddup(bx,1) + ddup(by,2) + ddup(bz,3)
  bx = bx + eps*dddn(divb,1)
  by = by + eps*dddn(divb,2)
  bz = bz + eps*dddn(divb,3)
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
  integer:: l(3), u(3), ix, iy, iz, i1, j1, k1
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%divb_clean2', itimer=itimer)
  bx => self%mem(:,:,:,self%idx%bx,self%it,1)
  by => self%mem(:,:,:,self%idx%by,self%it,1)
  bz => self%mem(:,:,:,self%idx%bz,self%it,1)
  l = self%mesh%lb
  u = self%mesh%ub-1
  i1 = merge(1, 0, self%n(1) > 1)
  j1 = merge(1, 0, self%n(2) > 1)
  k1 = merge(1, 0, self%n(3) > 1)
  !--- x-direction -------------------------------------------------------------
  do ix=self%mesh(1)%lo,self%mesh(1)%lb,-1
    bx(ix   ,l(2):u(2),l(3):u(3)) = &
    bx(ix+i1,l(2):u(2),l(3):u(3)) + ( &
      by(ix,l(2)+j1:u(2)+j1,l(3)   :u(3)   ) - by(ix,l(2):u(2),l(3):u(3)) + &
      bz(ix,l(2)   :u(2)   ,l(3)+k1:u(3)+k1) - bz(ix,l(2):u(2),l(3):u(3)))
  end do
  do ix=self%mesh(1)%ui,self%mesh(1)%ub-1
    bx(ix+i1,l(2):u(2),l(3):u(3)) = &
    bx(ix   ,l(2):u(2),l(3):u(3)) - ( &
      by(ix,l(2)+j1:u(2)+j1,l(3)   :u(3)   ) - by(ix,l(2):u(2),l(3):u(3)) + &
      bz(ix,l(2)   :u(2)   ,l(3)+k1:u(3)+k1) - bz(ix,l(2):u(2),l(3):u(3)))
  end do
  !--- y-direction -------------------------------------------------------------
  do iy=self%mesh(2)%lo,self%mesh(2)%lb,-1
    by(l(1):u(1),iy   ,l(3):u(3)) = &
    by(l(1):u(1),iy+j1,l(3):u(3)) + ( &
      bx(l(1)+i1:u(1)+i1,iy,l(3)   :u(3)   ) - bx(l(1):u(1),iy,l(3):u(3)) + &
      bz(l(1)   :u(1)   ,iy,l(3)+k1:u(3)+k1) - bz(l(1):u(1),iy,l(3):u(3)))
  end do
  do iy=self%mesh(2)%ui,self%mesh(2)%ub-1
    by(l(1):u(1),iy+j1,l(3):u(3)) = &
    by(l(1):u(1),iy   ,l(3):u(3)) - ( &
      bx(l(1)+i1:u(1)+i1,iy,l(3)   :u(3)   ) - bx(l(1):u(1),iy,l(3):u(3)) + &
      bz(l(1)   :u(1)   ,iy,l(3)+k1:u(3)+k1) - bz(l(1):u(1),iy,l(3):u(3)))
  end do
  !--- z-direction -------------------------------------------------------------
  do iz=self%mesh(3)%lo,self%mesh(3)%lb,-1
    bz(l(1):u(1),l(2):u(2),iz  ) = &
    bz(l(1):u(1),l(2):u(2),iz+k1) + ( &
      bx(l(1)+i1:u(1)+i1,l(2)   :u(2)   ,iz) - bx(l(1):u(1),l(2):u(2),iz) + &
      by(l(1)   :u(1)   ,l(2)+j1:u(2)+j1,iz) - by(l(1):u(1),l(2):u(2),iz))
  end do
  do iz=self%mesh(3)%ui,self%mesh(3)%ub-1
    bz(l(1):u(1),l(2):u(2),iz+k1) = &
    bz(l(1):u(1),l(2):u(2),iz   ) - ( &
      bx(l(1)+i1:u(1)+i1,l(2)   :u(2)   ,iz) - bx(l(1):u(1),l(2):u(2),iz) + &
      by(l(1)   :u(1)   ,l(2)+j1:u(2)+j1,iz) - by(l(1):u(1),l(2):u(2),iz))
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
  real, dimension(:,:,:), pointer:: d, e
  integer, save:: nprint=3
  !-----------------------------------------------------------------------------
  if (trim(self%eos)=='ideal') then
    d => self%mem(:,:,:,self%idx%d,self%new,1)
    e => self%mem(:,:,:,self%idx%e,self%new,1)
    if (eos%gamma==1d0) then
      pg = d*self%csound**2
    else
      pg = e*(eos%gamma-1)
    end if
  else
    call eos%lookup (shape(self%lnd), lnd=self%lnd, ee=self%ee, pg=pg)
    if (verbose>1) then
      print *, self%id, 'lookup_table: lnd =', minval(self%lnd), maxval(self%lnd)
      print *, self%id, 'lookup_table:  ee =', minval(self%ee ), maxval(self%ee )
      print *, self%id, 'lookup_table:  pg =', minval(pg),       maxval(pg)
    end if
  end if
  if (io%master .and. nprint > 0) then
    nprint = nprint-1
    print *, 'mhd_t%gas_pressure: min, max =', self%fminval(pg), self%fmaxval(pg)
  end if
END FUNCTION gas_pressure

END MODULE mhd_mod
