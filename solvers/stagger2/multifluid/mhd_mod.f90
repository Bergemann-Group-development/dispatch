!===============================================================================
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
  USE stagger_mod
  USE vector_ops
  USE gpatch_mod
  USE force_mod
  USE viscosity_mod
  USE index_mod
  USE timestep_mod
  USE eos_mod
  USE non_ideal_mod
  implicit none
  private
  type, extends(gpatch_t), public:: mhd_t
    real, dimension(:,:,:),   pointer:: d, dddt, s, dsdt, pg, dtr, ddtrdt           ! ADDED TRACER DENSITY AND TIME DERIVATIVE
    real, dimension(:,:,:,:), pointer:: p, dpdt, B, dBdt
    real, allocatable:: nu(:)
    real:: csound
    real:: cdtd
    integer:: nc=1                                                                  ! ADD NUMBER OF COMPONENTS nc
    integer, dimension(:), allocatable:: no_pressure                                ! ADDED VECTOR TO SPECIFY WHICH COMPONENTS ARE PRESSURELESS DUST
    logical:: mhd=.false., do_force=.true., &
              force_coupled=.false., do_viscosity=.true.                            ! ADDED force_coupled LOGICAL TO ENABLE A DIFFERENT INPUT TO FORCE MOD.
                                                                                    ! THE INPUT IS THEN, self%time, ALL OF self%mem, AND THE COMPONTENT CURRENTLY UNDER CONSIDERATION
    logical:: first_time=.true., tracer=.false.                                     ! ADDED TRACER
    type(force_t):: force
    type(viscosity_t):: viscosity
  contains
    procedure:: init
    procedure:: pde
    procedure:: update
    procedure:: output
    procedure:: gas_pressure
  end type
CONTAINS

!===============================================================================
!> The construction with a double test on first_time avoids entering a critical
!> region once the variable has been set.
!===============================================================================
SUBROUTINE init (self)
  class(mhd_t):: self
  integer:: i, k, ic                                                                    ! ADDED INTEGER FOR LOOPS OVER ALL COMPONENTS
  real, allocatable, save:: nu(:)
  real, save:: gamma=1.0, csound=1.0, cdtd=0.5, courant=0.2
  logical, save:: first_time=.true., tracer=.false., force_coupled=.false.              ! ADDED tracer = false AS DEFAULT
  character(len=16), save:: eos
  namelist /stagger_params/ nu, gamma, csound, cdtd, eos, tracer, & 
                        force_coupled, courant                                          ! ADDED tracer AND alpha TO THE NAMELIST
  !----------------------------------------------------------------------------
  call trace%begin('mhd_t%init')
  call io%lock%set
  if (first_time) then
    allocate(nu(6*self%nc))
    first_time = .false.
    eos = self%eos
    rewind (io%input)
    read (io%input, stagger_params)
    if (io%master) then
      write (*, stagger_params)
      write (io_unit%nml, stagger_params)
    end if
    call non_ideal%init
    call timestep%init
  end if
  call io%lock%unset
  self%eos           = eos
  self%gamma         = gamma
  self%csound        = csound
  self%courant       = courant
  self%cdtd          = cdtd
  self%force_coupled = force_coupled
  allocate(self%nu(6*self%nc))
  self%nu = nu
  allocate(self%no_pressure(1))
  self%no_pressure = 0
  call self%idx%init(5, self%mhd)
  if (self%mhd) then
    self%kind = 'stagger2_mhd_patch'
  else
    self%kind = 'stagger2_hd_patch'
  end if
  if (tracer) then
    self%tracer  = tracer                         ! ADDED TRACER TRUTH STATEMENT TO self
    self%nv      = self%nv+1                      ! INCREASE NUMBER OF VARIABLES TO MAKE ROOM FOR PASSIVE TRACER
    self%idx%dtr = self%nv                        ! ADDED INDEX FOR TRACER DENSITY. CURRENTLY, ONLY THE FIRST 
  end if                                          ! FLUID MAY BE ASSIGNED A TRACER, AND THUS ADDITIONAL FLUIDS
                                                  ! ALWAYS HAVE 5 COMPONENTS ONLY
  if (force_coupled) self%force_coupled = force_coupled
  call self%gpatch_t%init
  self%type = 'mhd_t'
  do k=1,self%nc                                  ! LOOP OVER ALL COMPONENTS
    ic = 5*(k-1)                                  ! ADDED INDEX TO CYCLE TO NEXT FLUID
    do i=1,3
        self%mesh(i)%h(self%idx%px-1+i+ic) = -0.5
        if (self%mhd) self%mesh(i)%h(self%idx%bx-1+i) = -0.5
    end do
    self%unsigned(idx%d+ic)  = .true.
    self%pervolume(idx%s+ic) = .true.
  end do
  call self%force%init(self%kind, self%id, self%mesh)
  call self%viscosity%init(self%kind)
  call trace%end()
END SUBROUTINE init

!===============================================================================
SUBROUTINE pde (self)
  USE scalar_mod
  class(mhd_t):: self
  !.............................................................................
  real, dimension(:,:,:,:,:), pointer:: T                                                       ! ADEDD TENSOR
  real, dimension(:,:,:,:), pointer:: U, E, J, p, B, dpdt, dBdt, Urot, plin                     ! ADDED FRAME VELOCITY WITH 3 DIFFERENT CENTERINGS
  real, dimension(:,:,:,:), pointer:: Ud, E_AD                                                  ! ADDED POINTER TO ALL DENSITIES AND ALL MOMENTA
  real, dimension(:,:,:,:), pointer:: fx, fy, fz, fd, ld, dd, fm, fdtr, ldtr, ddtr, fmtr        ! ADDED TEMPORARY TRACER POINTERS
  real, dimension(:,:,:),   pointer:: d, s, dddt, dsdt, dtr, ddtrdt                             ! ADDED NEW POINTERS FOR TRACER dt AND ddtdt
  real, dimension(:,:,:),   pointer:: du, lnd, ss, cs, pg, pa, pp, pb, Q, lndtr, dnu            ! ADDED TEMPORARY TRACER POINTER AND VISCOSITY PARAMETER
  real, dimension(:,:,:),   pointer:: dpdtx, dpdty, dpdtz, dBdtx, dBdty, dBdtz
  real, dimension(:,:,:),   pointer:: Bx, By, Bz, px, py, pz, phi
  real, dimension(:,:,:),   pointer:: Ux, Uy, Uz, Ex, Ey, Ez, Jx, Jy, Jz
  real, dimension(:,:,:),   pointer:: Urotx, Uroty, Urotz                                       ! ADDED POINTERS TO DIFFERENT CENTERINGS OF Urot
  real, dimension(:,:,:),   pointer:: cx1, cx2, cy1, cy2, cy3, cy4, cz1, cz2, cz3               ! ADDED ADDITIONAL POINTERS FOR USE IN COROTATIONAL FRAME
  real, dimension(:,:,:),   pointer:: fdx, fdy, fdz, ldx, ldy, ldz, ddx, ddy, ddz
  real, dimension(:,:,:),   pointer:: divU
  real, dimension(:,:,:),   pointer:: fxx, fxy
  real, dimension(:,:,:),   pointer::      fyy, fyz
  real, dimension(:,:,:),   pointer:: fzx,      fzz
  real, dimension(:,:,:),   pointer:: fyyx, fzzx, fzzy                                          ! ADDED SEPERATE TERMS FOR CENTRIFUGAL FORCES
  real, dimension(:,:,:),   pointer:: Txx, Txy
  real, dimension(:,:,:),   pointer::      Tyy, Tyz
  real, dimension(:,:,:),   pointer:: Tzx,      Tzz
  integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
  logical, save:: first_time=.true., first_AD=.true.
  !.............................................................................
  real(8):: ds(3), ds_min(3), dsmax, dsmin
  integer, save:: itimer=0
  real:: u_max
  integer:: k, ic, i2                                                                           ! ADDED INTEGER FOR LOOPS OVER ALL COMPONENTS
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%pde',itimer=itimer)
  self%u_max = 0.0                                                                              ! INITIALIZE u_max
  if (self%mhd) then
    self%mem(:,:,:,self%idx%bx:self%idx%bz,self%it,2) = 0.0                                                   ! INITIALIZE dBdt
  end if 
  do k=1,self%nc                                                                                ! ADDED LOOP OVER ALL COMPONENTS
  ic = 5*(k-1)                                                                                  ! ADDED INDEX TO CYCLE TO NEXT FLUID

  self%d => self%mem(:,:,:,self%idx%d+ic         ,self%it,1); self%dddt => self%mem(:,:,:,self%idx%d+ic         ,self%it,2)
  self%s => self%mem(:,:,:,self%idx%s+ic         ,self%it,1); self%dsdt => self%mem(:,:,:,self%idx%s+ic         ,self%it,2)
  self%p => self%mem(:,:,:,self%idx%px+ic:self%idx%pz+ic,self%it,1); self%dpdt => self%mem(:,:,:,self%idx%px+ic:self%idx%pz+ic,self%it,2)
  d=>self%d; dddt=>self%dddt; s=>self%s; dsdt=>self%dsdt;
  p=>self%p; dpdt=>self%dpdt
  px=>p(:,:,:,1); py=>p(:,:,:,2); pz=>p(:,:,:,3);
  dpdtx=>dpdt(:,:,:,1); dpdty=>dpdt(:,:,:,2); dpdtz=>dpdt(:,:,:,3);
  if (self%mhd) then
    self%B => self%mem(:,:,:,self%idx%bx:self%idx%bz,self%it,1); self%dBdt => self%mem(:,:,:,self%idx%bx:self%idx%bz,self%it,2)
    B=>self%B; dBdt=>self%dBdt
    Bx=>B(:,:,:,1); By=>B(:,:,:,2); Bz=>B(:,:,:,3);
    dBdtx=>dBdt(:,:,:,1); dBdty=>dBdt(:,:,:,2); dBdtz=>dBdt(:,:,:,3);
  end if
  if (self%tracer .and. k==1) then                                                              ! ADDED POINTER TO TRACER IF tracer IS TRUE
    self%dtr => self%mem(:,:,:,self%idx%dtr,self%it,1); self%ddtrdt => self%mem(:,:,:,self%idx%dtr,self%it,2)
    dtr => self%dtr; ddtrdt => self%ddtrdt
    call allocate_vectors (self%gn, ldtr, ddtr, fdtr, fmtr)
    call allocate_scalars (self%gn, lndtr)
  end if
  call allocate_vectors (self%gn, ld, dd, E, J, fd, fx, fy, fz, fm)
  call allocate_vectors (self%gn, Urot, plin)
  call allocate_scalars (self%gn, lnd, pa, du, ss, cs, pp, pb, Q)
  call allocate_scalars (self%gn, fyyx, fzzx, fzzy, dnu, divU)
  if (self%do_pebbles) then
    U => self%vgas(:,:,:,:,self%it)
    pg => self%pressure(:,:,:,self%it)
    self%density => self%mem(:,:,:,1,:,1)
  else
    call allocate_scalars (self%gn, pg)
    call allocate_vectors (self%gn, U)
  end if
  Ux=>U(:,:,:,1); Uy=>U(:,:,:,2); Uz=>U(:,:,:,3)
  Urotx=>Urot(:,:,:,1); Uroty=>Urot(:,:,:,2); Urotz=>Urot(:,:,:,3)                              ! Urot is always in the 3 direction. The difference is the centering

  Ex=>E(:,:,:,1); Ey=>E(:,:,:,2); Ez=>E(:,:,:,3)
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
  lnd = log(d)
  ld  = down(lnd)
  dd  = exp(ld)
  plin = hinv_vector(self%mesh,p)
  U   = plin/dd
  !-----------------------------------------------------------------------------
  ! Gradient of velocity vector. The resulting rank 2 tensor contains all of the
  ! relevant velocity derivatives necessary to compute the viscous stress tensor
  ! Must be computed before background velocity is subtracted (if co-rotating)
  !-----------------------------------------------------------------------------
  allocate(T(self%gn(1),self%gn(2),self%gn(3),3,3))
  T  = grad_vector(self%mesh,U)
  !-----------------------------------------------------------------------------
  ! Compute and subtract background velocity to convert to co-rotating frame
  !-----------------------------------------------------------------------------
  Urot  = 1.0
  Urotx = self%grid_vel(3)*h31f(self%mesh,Urotx)*h32c(self%mesh,Urotx)                          ! Compute Urot centered on x-faces.
  Uroty = self%grid_vel(3)*h31c(self%mesh,Uroty)*h32f(self%mesh,Uroty)                          ! Compute Urot centered on y-faces.
  Urotz = self%grid_vel(3)*h3c(self%mesh,Urotz)
  Uz  = Uz - Urotz
  !-----------------------------------------------------------------------------
  ! Gas pressure. The 2nd table variable is entropy per unit mass, effectively
  ! a log variable, from the table point of view.
  !-----------------------------------------------------------------------------
  ss = s/d
  if (any(self%no_pressure.EQ.k)) then                   ! ADDED no_pressure CONDITION
    pg(:,:,:) = 0d0
  else
    pg = self%gas_pressure(d,s)
  end if
  self%pg => pg
  !-----------------------------------------------------------------------------
  ! Divergence of velocity, and artificial pressure.  The divergence is valid
  ! in (2,gn-1)
  !-----------------------------------------------------------------------------
  ds = self%ds
  ds(2) = ds(2)*self%mesh(1)%h2c(self%mesh(1)%ub)
  ds(3) = ds(3)*self%mesh(1)%h31c(self%mesh(1)%ub)*abs(self%mesh(2)%h32c(self%mesh(2)%lb))
  ds_min(1) = ds(1)
  ds_min(2) = ds(2)*self%mesh(1)%h2c(self%mesh(1)%lb)
  ds_min(3) = ds(3)*self%mesh(1)%h31c(self%mesh(1)%lb)*abs(self%mesh(2)%h32c(self%mesh(2)%lb))
  if (trim(self%mesh(3)%type) == 'spherical_mesh') then
    do i2=self%mesh(2)%lb,self%mesh(2)%ub
      ds(3)     = max(ds(3)    ,self%ds(3)*self%mesh(1)%h31c(self%mesh(1)%ub)*abs(self%mesh(2)%h32c(i2)))
      ds_min(3) = min(ds_min(3),self%ds(3)*self%mesh(1)%h31c(self%mesh(1)%lb)*abs(self%mesh(2)%h32c(i2)))
    end do
  end if
  dsmin = minval(ds,self%n > 1)                         ! smallest {dx,dy,dz}
  dsmax = maxval(ds,self%n > 1)                         ! largest {dx,dy,dz}
  divU = div(self%mesh,U)
  du = -dsmax*divU                                      ! convergence rate
  du = max(du,0.0)                                      ! positive part
  !-----------------------------------------------------------------------------
  ! Conservation of momentum -- first the diagonal part of the stress tensor
  ! The diagonal fxx is valid in (2:gn-1)
  !-----------------------------------------------------------------------------
  call allocate_scalars (self%gn, Txx, Tyy, Tzz)
  Txx = T(:,:,:,1,1)
  Tyy = T(:,:,:,2,2)
  Tzz = T(:,:,:,3,3)
  if (self%mhd) then
    ! Diagonal Reynolds stress
    fxx = xup(Ux)**2
    fyy = yup(Uy)**2
    fzz = zup(Uz)**2
    ! Face centered centrifugal terms
    fyyx = ddx*xdn(yup(Uy      ))**2 - xdn(yup(By))**2
    fzzx = ddx*xdn(zup(Uz+Urotz))**2 - xdn(zup(Bz))**2
    fzzy = ddy*ydn(zup(Uz+Urotz))**2 - ydn(zup(Bz))**2
    pp  = (fxx+fyy+fzz)
    ! Diagonal (negative) Maxwell stress (borrowing E for scratch)
    Ex = xup(Bx)**2
    Ey = yup(By)**2
    Ez = zup(Bz)**2
    ! Magnetic pressure
    pb = 0.5*(Ex+Ey+Ez)
    ! Fast mode speed + velocity
    pa = self%nu(2+6*(k-1))*d*du**2
    cs = sqrt((self%gamma*pg+2.0*pb + pa)/d)
    if (io%verbose>2) then
      u_max = self%fmaxval(cs)
      self%u_max = max(self%u_max,u_max)
      print 1, 'u_max(ca)  ', u_max, self%u_max
    1 format(1x,a,1p,2g12.3)
    end if
    u_max = self%fmaxval(cs + sqrt(pp))
    self%u_max = max(self%u_max,u_max)
    if (io%verbose>2) print 1, 'u_max(ca+u)', u_max, self%u_max
    ! Isotropic pressure part
    pp = pg + pa + pb
    ! Physical viscosity
    if (associated(self%viscosity%selected) .and. self%do_viscosity) then
      dnu = self%viscosity%selected(self%time, self%mem(:,:,:,:,self%it,1), U, cs, self%mesh)
    else
      dnu = 0.0
    end if
    ! Total diagonal stress
    du = dsmax*self%nu(1+6*(k-1))*cs
    fxx = d*(fxx - (du+2.0*dnu)*Txx + (2./3.)*dnu*divU) - Ex
    fyy = d*(fyy - (du+2.0*dnu)*Tyy + (2./3.)*dnu*divU) - Ey
    fzz = d*(fzz - (du+2.0*dnu)*Tzz + (2./3.)*dnu*divU) - Ez
    ! Addition of stress tensor terms to the centrifugal terms
    fyyx = fyyx - ddx*xdn(du+2.0*dnu)*xdn(Tyy) + (2./3.)*ddx*xdn(dnu)*xdn(divU) - xdn(Ey)
    fzzx = fzzx - ddx*xdn(du+2.0*dnu)*xdn(Tzz) + (2./3.)*ddx*xdn(dnu)*xdn(divU) - xdn(Ez)
    fzzy = fzzy - ddy*ydn(du+2.0*dnu)*ydn(Tzz) + (2./3.)*ddy*ydn(dnu)*ydn(divU) - ydn(Ez)
  else
    ! Diagonal Reynolds stress
    fxx = xup(Ux)**2
    fyy = yup(Uy)**2
    fzz = zup(Uz)**2
    ! Face centered centrifugal terms
    fyyx = ddx* xdn(yup(Uy))**2
    fzzx = ddx*(xdn(zup(Uz+Urotz)))**2
    fzzy = ddy*(ydn(zup(Uz+Urotz)))**2
    pp = (fxx+fyy+fzz)
    ! Isotropic pressure part
    pa = self%nu(2+6*(k-1))*d*du**2                               ! dimension: pressure
    ! Sound speed + velocity
    cs = sqrt((self%gamma*pg+pa)/d)
    if (io%verbose>2) then
      u_max = self%fmaxval(cs)
      self%u_max = max(self%u_max,u_max)
      print 1, 'u_max(cs)  ', u_max, self%u_max
    end if
    u_max = self%fmaxval(cs + sqrt(pp))
    self%u_max = max(self%u_max,u_max)
    if (io%verbose>2) print 1, 'u_max(cs+u)', u_max, self%u_max
    ! Total diagonal stress, with density factor
    du = dsmax*self%nu(1+6*(k-1))*cs
    pp = pg + pa
    ! Physical viscosity
    if (associated(self%viscosity%selected) .and. self%do_viscosity) then
      dnu = self%viscosity%selected(self%time, self%mem(:,:,:,:,self%it,1), U, cs, self%mesh)
    else
      dnu = 0.0
    end if
    !---------------------------------------------------------------------------
    ! Incoming fxx is valid in (2:gn-1), as is du and ddxup(Ux)
    !---------------------------------------------------------------------------
    fxx = d*(fxx - (du+2.0*dnu)*Txx + (2./3.)*dnu*divU)
    fyy = d*(fyy - (du+2.0*dnu)*Tyy + (2./3.)*dnu*divU)
    fzz = d*(fzz - (du+2.0*dnu)*Tzz + (2./3.)*dnu*divU)
    ! Addition of stress tensor terms to the centrifugal terms
    fyyx = fyyx - ddx*xdn(du+2.0*dnu)*xdn(Tyy) + (2./3.)*ddx*xdn(dnu)*xdn(divU)
    fzzx = fzzx - ddx*xdn(du+2.0*dnu)*xdn(Tzz) + (2./3.)*ddx*xdn(dnu)*xdn(divU)
    fzzy = fzzy - ddy*ydn(du+2.0*dnu)*ydn(Tzz) + (2./3.)*ddy*ydn(dnu)*ydn(divU)
  end if
  !-----------------------------------------------------------------------------
  ! Viscous time scale evaluations.  Viscous terms are ~ 1/ds^2, and while the
  ! value of the viscosity scales with dsmax, the diffusion speed scales with
  ! 1/dsmin.
  !-----------------------------------------------------------------------------
  u_max = self%fmaxval(du)*dsmax/dsmin**2
  self%u_max = max(self%u_max,u_max)
  if (io%verbose>1) print 1, 'u_max(du)', u_max, self%u_max
  u_max = self%fmaxval(dnu)*dsmax/dsmin**2
  self%u_max = max(self%u_max,u_max)
  if (io%verbose>1) print 1, 'u_max(dnu)', u_max, self%u_max
  !-----------------------------------------------------------------------------
  ! du has units area per unit time, which gives velocity squared times mass
  ! density per unit time = energy per unit volume and time
  !-----------------------------------------------------------------------------
  Q = d*du*(Txx**2 + Tyy**2 + Tzz**2) + &
      2.0*d*dnu*( (Txx-divU/3.0)**2 + (Tyy-divU/3.0)**2 + (Tzz-divU/3.0)**2 )
  call deallocate_scalars (Txx, Tyy, Tzz)
  !-----------------------------------------------------------------------------
  ! Conservation of mass
  !-----------------------------------------------------------------------------
  fd = -(self%nu(4+6*(k-1))*dd*down(du))*grad(self%mesh,lnd)   ! diffusive mass flux
  fm = plin + fd
  fm(:,:,:,3) = fm(:,:,:,3) - ddz*Urotz
  dddt = -div(self%mesh,fm)                                    ! density derivative
  ! viscosity related diffusivity
  if (associated(self%viscosity%selected2) .and. self%do_viscosity) then
    dddt = dddt + self%viscosity%selected2(self%time, self%mem(:,:,:,:,self%it,1), U, self%mesh, self%force)
  end if
  u_max = self%cdtd*dsmax*self%fmaxval(abs(dddt/d))             ! ops: 0 0 1
  self%u_max = max(self%u_max,u_max)
  if (io%verbose>2) print 1, 'u_max(dddt)', u_max, self%u_max
  if (self%tracer .and. k==1) then                                     ! ADDED COLOR EQUATION FOR TRACER DENSITY
    lndtr  = log(dtr)
    ldtr   = down(lndtr)
    ddtr   = exp(ldtr)
    fdtr   = -(self%nu(4+6*(k-1))*ddtr*down(du))*grad(self%mesh,lndtr) ! diffusive mass flux of tracer
    fmtr   = (U*ddtr)+fdtr                                             ! mass flux of tracer
    ddtrdt = -div(self%mesh,fmtr) + dtr*divU                           ! density derivative of tracer
  end if
  call deallocate_scalars (divU)
  !---------------------------------------------------------------------------
  ! Off-diagonal part of Reynolds stress, for now without the density factor.
  ! Note that, because of the symmetry, we define only three components
  !---------------------------------------------------------------------------
  fxy = ydn(Ux)*xdn(Uy)
  fyz = zdn(Uy)*ydn(Uz)
  fzx = xdn(Uz)*zdn(Ux)
  !---------------------------------------------------------------------------
  ! Add viscous part -- edge-centered
  !---------------------------------------------------------------------------
  call allocate_scalars (self%gn, Txy, Tyz, Tzx)
  Txy = T(:,:,:,1,2) + T(:,:,:,2,1)
  Tyz = T(:,:,:,2,3) + T(:,:,:,3,2)
  Tzx = T(:,:,:,3,1) + T(:,:,:,1,3)
  deallocate(T)
  fxy = fxy - xdn1(ydn1(self%nu(3+6*(k-1))*0.5*du+dnu))*Txy
  fyz = fyz - ydn1(zdn1(self%nu(3+6*(k-1))*0.5*du+dnu))*Tyz
  fzx = fzx - zdn1(xdn1(self%nu(3+6*(k-1))*0.5*du+dnu))*Tzx
  Q   = Q   + d*(self%nu(3+6*(k-1))*0.5*du+dnu)*(xup(yup(Txy))**2 + &
                                                 yup(zup(Tyz))**2 + &
                                                 zup(xup(Tzx))**2)
  call deallocate_scalars (Txy, Tyz, Tzx)
  !---------------------------------------------------------------------------
  ! Multiply with edge-centered density
  !---------------------------------------------------------------------------
  fxy = exp(ydn(ldx))*fxy
  fyz = exp(zdn(ldy))*fyz
  fzx = exp(xdn(ldz))*fzx
  !---------------------------------------------------------------------------
  ! Off-diagonal Maxwell stress
  !---------------------------------------------------------------------------
  if (self%mhd) then
    fxy = fxy - ydn(Bx)*xdn(By)
    fyz = fyz - zdn(By)*ydn(Bz)
    fzx = fzx - xdn(Bz)*zdn(Bx)
  end if
  !-----------------------------------------------------------------------------
  ! Consistent addition of mass flux related momentum flux, symmetrized. The
  ! diagonal part is analogous, a centered version of \nu \rho U_i \grad\ln\rho,
  ! while the off-diagonal part needs to be properly symmetrized.
  !-----------------------------------------------------------------------------
  if (self%nu(4+6*(k-1)) > 0.0) then
    pa = d*du*self%nu(4+6*(k-1))
    !---------------------------------------------------------------------------
    ! fxx, pa, ddxup(ldx), and xup(Ux) are all valid in (2:gn-1)
    !---------------------------------------------------------------------------
    fxx = fxx - pa*ddxup1(self%mesh,ldx)*xup(Ux)
    fyy = fyy - pa*ddyup1(self%mesh,ldy)*yup(Uy)
    fzz = fzz - pa*ddzup1(self%mesh,ldz)*zup(Uz+Urotz)
    fxy = fxy + (xdn1(fdy)*ydn(Ux)+ydn1(fdx)*xdn(Uy))
    fyz = fyz + (ydn1(fdz)*zdn(Uy)+zdn1(fdy)*ydn(Uz+Urotz))
    fzx = fzx + (zdn1(fdx)*xdn(Uz)+xdn1(fdz)*zdn(Ux))
    ! Addition of mass flux related momentum flux to the centrifugal terms
    pa   = ddx*xdn(du)*self%nu(4+6*(k-1))
    fyyx = fyyx - pa*ddyup1(self%mesh,xdn(ldy))*xdn(yup(Uy))
    fzzx = fzzx - pa*ddzup1(self%mesh,xdn(ldz))*xdn(zup(Uz+Urotz))
    pa   = ddy*ydn(du)*self%nu(4+6*(k-1))
    fzzy = fzzy - pa*ddzup1(self%mesh,ydn(ldz))*ydn(zup(Uz+Urotz))
  end if
  !-----------------------------------------------------------------------------
  ! Momentum equations.  The ddxdn(fxx) part is valid in (3:gn-1), and hence px
  ! only needs guard zones in (1:2) and (gn:gn)
  !-----------------------------------------------------------------------------
  dpdtx = - divx(self%mesh,fxx,fxy,fzx,fyyx,fzzx)
  dpdty = - divy(self%mesh,fxy,fyy,fyz,fzzy)
  dpdtz = - divz(self%mesh,fzx+zdn(plin(:,:,:,1)          )*Urotx, &
                           fyz+zdn(plin(:,:,:,2)          )*Uroty, &
                           fzz+zup(plin(:,:,:,3)-ddz*Urotz)*Urotz)
  dpdt  =   dpdt - vector_scalar_grad(self%mesh,pp)
  !-----------------------------------------------------------------------------
  ! External force; e.g. point source gravity.  Note that the density factor is
  ! applied here, and should not be included in the function definition
  !-----------------------------------------------------------------------------
  if (associated(self%force%selected) .and. self%do_force) then
    if (self%force_coupled) then
      dpdt = dpdt + dd*self%force%selected (self%time, d, self%mem(:,:,:,:,self%it,1), Ux, Uy, Uz, self%mesh)
      u_max = dsmax/self%force%tmin
      self%u_max = max(self%u_max,u_max)
      if (io%verbose>2) print 1, 'u_max(drag)', u_max, self%u_max
    else
      dpdt = dpdt + dd*self%force%selected (self%time, d, p, Ux, Uy, Uz, self%mesh)
    end if
  end if
  if (self%selfgravity) then
    phi => self%mem(:,:,:,self%idx%phi,self%it,1)
    dpdt = dpdt - dd*vector_scalar_grad(self%mesh,phi)
  end if
  ! viscosity related dispersion
  if (associated(self%viscosity%selected3) .and. self%do_viscosity) then
    dpdt = dpdt + self%viscosity%selected3(self%time, self%mem(:,:,:,:,self%it,1), U, self%mesh, self%force)
  end if
  !-----------------------------------------------------------------------------
  ! Conservation of entropy.  The formulation that works best so far is to make
  ! the diffusive entropy flux directly proportional to the diffusive mass flux,
  ! plus a contribution proportional to the gradient of specific entropy ss.
  !-----------------------------------------------------------------------------
  if (any(self%no_pressure.EQ.k)) then                           ! ADDED no_pressure CONDITION
    dsdt(:,:,:) = 0d0
  else
    if (self%gamma==1.0) then
        dsdt = 0.0
    else
        fm = fm*down(ss) - (self%nu(5+6*(k-1))*dd*down(du))*grad(self%mesh,ss)
        !---------------------------------------------------------------------------
        ! Q is heating per unit volume (energy per unit volume and time), which when
        ! divided by presssure (energy per unit volume) becomes entropy (unitless)
        ! per unit mass.  Multiplied by mass density it becomes mass-entropy per
        ! unit volume, which is our units for s.
        !---------------------------------------------------------------------------
        dsdt = -div(self%mesh,fm) + Q*d/pg
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Induction equation
  !-----------------------------------------------------------------------------
  if (self%mhd) then
    !---------------------------------------------------------------------------
    ! Electric current, resistive electric field and magnetic dissipation
    !---------------------------------------------------------------------------
    J = curl_dn(self%mesh,B)
    E = self%nu(6+6*(k-1))*edge(du)*J                                           ! mhd: 0 6 0
    Q = yup(zup(Jx*Ex)) + zup(xup(Jy*Ey)) + xup(yup(Jz*Ez))                     ! mhd: 2 3 0
    !---------------------------------------------------------------------------
    ! Advective electric field on cell edges
    !---------------------------------------------------------------------------
    Ex = Ex + ydn(Uz)*zdn(By) - zdn(Uy)*ydn(Bz)                                 ! mhd: 2 2 0
    Ey = Ey + zdn(Ux)*xdn(Bz) - xdn(Uz)*zdn(Bx)                                 ! mhd: 2 2 0
    Ez = Ez + xdn(Uy)*ydn(Bx) - ydn(Ux)*xdn(By)                                 ! mhd: 2 2 0
    call non_ideal%update (self%gn, minval(ds_min), d, Jx, Jy, Jz, Bx, By, Bz, Q, E, self%u_max)
    if (self%gamma/=1.0) then
      dsdt = dsdt + Q*d/pg                                                      ! mhd: 1 1 1
    end if
    dBdt = dBdt - curl_up(self%mesh,E)
  end if
  !-----------------------------------------------------------------------------
  ! Time step limitation for rate of change dS/dt = (ds/dt-(s/d)*dddt)/d
  !-----------------------------------------------------------------------------
  u_max = self%cdtd*dsmax*self%fmaxval(abs((dsdt-s*dddt/d)/d))  ! ops: 1 1 2
  self%u_max = max(self%u_max,u_max)
  if (io%verbose>2) print 1, 'u_max(dsdt)', u_max, self%u_max
  !-----------------------------------------------------------------------------
  call deallocate_vectors (ld, dd, E, J, fd, fx, fy, fz, fm)
  call deallocate_vectors (Urot, plin)
  call deallocate_scalars (lnd, pa, du, ss, cs, pp, pb, Q)
  call deallocate_scalars (fyyx, fzzx, fzzy, dnu)
  if (self%tracer .and. k==1) then                             ! ADDED DEALLOCATION OF PARAMETERS IF tracer IS TRUE
    call deallocate_vectors (ldtr, ddtr, fdtr, fmtr)
    call deallocate_scalars (lndtr)
  end if
  if (.not.self%do_pebbles) then
    call deallocate_vectors (U)
    call deallocate_scalars (pg)
  end if
  end do
  call trace%end (itimer)
END SUBROUTINE pde

!===============================================================================
!> Take a complete update step.  This can be overloaded in higher levels, e.g.
!> in experiment_mod.
!===============================================================================
SUBROUTINE update (self)
  class(mhd_t):: self
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%update',itimer=itimer)
  !call self%gravity_solve
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
  !call self%gravity_update
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Add extra, MHD-related information to the patch .txt output file
!===============================================================================
SUBROUTINE output (self)
  class(mhd_t):: self
  character(len=64):: filename
  logical:: exists
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (self%time < self%out_next .or. .not.io%do_output) return
  call trace%begin('mhd_t%output',itimer=itimer)
  call self%gpatch_t%output
  if (io%do_legacy) then
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
  call trace%end (itimer)
END SUBROUTINE output

!===============================================================================
FUNCTION gas_pressure (self, d_opt, s_opt) RESULT (pg)
  class(mhd_t):: self
  real, dimension(:,:,:), pointer:: d, s
  real, dimension(:,:,:), pointer, optional:: d_opt, s_opt
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: pg
  !-----------------------------------------------------------------------------
  if (trim(self%eos)=='ideal') then
    if (present(d_opt)) then
      d => d_opt
    else
      d => self%mem(:,:,:,self%idx%d,self%it,1)
    end if
    if (present(s_opt)) then
      s => s_opt
    else
      s => self%mem(:,:,:,self%idx%s,self%it,1)
    end if
    if (self%gamma==1.0) then
      pg = d*self%csound**2
    else
      pg = d**self%gamma*exp(s/d*(self%gamma-1))
    end if
  else
    call eos%lookup (shape(d), x=s, y=d, pg=pg)
  end if
END FUNCTION gas_pressure
END MODULE mhd_mod
