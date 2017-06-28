!> =============================================================================
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
  USE stagger_mod
  USE gpatch_mod
  USE force_mod
  USE eos_mod
  USE index_mod
  USE timestep_mod
  USE planet_mod
  implicit none
  private
  type, extends(gpatch_t), public:: mhd_t
    real, dimension(:,:,:),   pointer:: d, dddt, e, dedt
    real, dimension(:,:,:,:), pointer:: p, dpdt, B, dBdt
    real:: nu(6), csound
    real:: cdtd
    logical:: mhd=.true.
    type(force_t):: force
  contains
    procedure:: init
    procedure:: pde
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
    procedure:: cell_gas_pressure
    procedure:: gas_temperature
    procedure:: cell_gas_temperature
    procedure:: velocity_magnitude => void
    procedure:: magnetic_field_magnitude => void
    procedure:: grav_potential => void
  end type
  type(index_t):: idx
CONTAINS

!> =============================================================================
!> The construction with a double test on first_time avoids entering a critical
!> region once the variable has been set.
!> =============================================================================
SUBROUTINE init (self)
  class(mhd_t):: self
  real, save:: nu(6)=[0.1,1.0,0.0,0.5,0.5,0.5], gamma=1.0, csound=1.0, cdtd=0.5
  integer:: i
  logical, save:: first_time=.true.
  namelist /mhd_params/ nu, gamma, csound, cdtd, hardwire
  !----------------------------------------------------------------------------
  call trace%begin('mhd_t%init')
  if (self%mhd) then
    self%kind = 'stagger2e_mhd_patch'
  else
    self%kind = 'stagger2e_hd_patch'
  end if
  call self%gpatch_t%init
  self%i_d  = idx%d
  self%i_px = idx%px
  self%i_py = idx%py
  self%i_pz = idx%pz
  self%i_e  = idx%e
  self%i_bx = idx%bx
  self%i_by = idx%by
  self%i_bz = idx%bz
  !$omp critical (mhd_input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, mhd_params)
    if (io%master) write (*, mhd_params)
    call timestep%init
  end if
  !$omp end critical (mhd_input_cr)
  self%nu = nu
  self%gamma = gamma
  self%csound = csound
  self%cdtd = cdtd
  do i=1,3
    self%mesh(i)%h(2+i) = -0.5
    if (self%nv > 5) self%mesh(i)%h(5+i) = -0.5
  end do
  !self%unsigned(idx%d) = .true.
  !self%unsigned(idx%e) = .true.
  call self%force%init(self%kind)
  call trace%end
END SUBROUTINE init

!> =============================================================================
SUBROUTINE pde (self)
  class(mhd_t):: self
  !.............................................................................
  real, dimension(:,:,:,:), pointer:: U, emf, J, p, B, dpdt, dBdt
  real, dimension(:,:,:,:), pointer:: fx, fy, fz, fd, ld, dd, fm
  real, dimension(:,:,:),   pointer:: d, e, dddt, dedt
  real, dimension(:,:,:),   pointer:: du, lnd, ee, cs, pg, pa, pp, pb, Q
  real, dimension(:,:,:),   pointer:: dpdtx, dpdty, dpdtz, dBdtx, dBdty, dBdtz
  real, dimension(:,:,:),   pointer:: Bx, By, Bz, px, py, pz
  real, dimension(:,:,:),   pointer:: Ux, Uy, Uz, Ex, Ey, Ez, Jx, Jy, Jz
  real, dimension(:,:,:),   pointer:: fdx, fdy, fdz, ldx, ldy, ldz, ddx, ddy, ddz
  real, dimension(:,:,:),   pointer:: fxx, fxy
  real, dimension(:,:,:),   pointer::      fyy, fyz
  real, dimension(:,:,:),   pointer:: fzx,      fzz
  !.............................................................................
  real:: ds(3), dx
  integer, save:: itimer=0
  type(stagger_t):: st
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%pde',itimer=itimer)
  st%m => self%mesh
  self%d    => self%mem(:,:,:,self%i_d,self%it,1)
  self%dddt => self%mem(:,:,:,self%i_d,self%it,2)
  self%e    => self%mem(:,:,:,self%i_e,self%it,1)
  self%dedt => self%mem(:,:,:,self%i_e,self%it,2)
  self%p    => self%mem(:,:,:,self%i_px:self%i_pz,self%it,1)
  self%dpdt => self%mem(:,:,:,self%i_px:self%i_pz,self%it,2)
  d=>self%d; dddt=>self%dddt
  e=>self%e; dedt=>self%dedt
  p=>self%p; dpdt=>self%dpdt
  px=>p(:,:,:,1); py=>p(:,:,:,2); pz=>p(:,:,:,3);
  dpdtx=>dpdt(:,:,:,1); dpdty=>dpdt(:,:,:,2); dpdtz=>dpdt(:,:,:,3);
  if (self%mhd) then
    self%B    => self%mem(:,:,:,self%i_bx:self%i_bz,self%it,1)
    self%dBdt => self%mem(:,:,:,self%i_bx:self%i_bz,self%it,2)
    B=>self%B; dBdt=>self%dBdt
    Bx   =>   B(:,:,:,1); By   =>   B(:,:,:,2); Bz   =>   B(:,:,:,3);
    dBdtx=>dBdt(:,:,:,1); dBdty=>dBdt(:,:,:,2); dBdtz=>dBdt(:,:,:,3);
  end if
  call allocate_vectors (self%gn, ld, dd, emf, J, fd, fx, fy, fz, fm)
  call allocate_scalars (self%gn, lnd, pa, du, ee, cs, pp, pb, Q)
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
  ! Compute velocity U=momentum/density, leaving down shifted density and log
  ! density for reuse below
  !-----------------------------------------------------------------------------
  lnd = log(d)
  ld = sUP(lnd)
  dd = exp(ld)
  U = p/dd
  !-----------------------------------------------------------------------------
  ! Gas pressure
  !-----------------------------------------------------------------------------
  ee = e/d
  if (self%gamma==1.0) then
    pg = d*self%csound**2
  else
    call eos%lookup (shape(ee), x=ee, lny=lnd, pg=pg)
  end if
  !-----------------------------------------------------------------------------
  ! Divergence of velocity, and artificial pressure
  !-----------------------------------------------------------------------------
  ds = self%mesh%d
  dx = minval(ds)                                       ! smallest {dx,dy,dz}
  du = st%div_DN(U)                                     ! save div(U) for now
  if (self%gamma /= 1.0) then                           ! unless isothermal
    dedt = -pg*du                                       ! -Pg*div(U) term
  else
    dedt = 0.0
  end if
  du = -dx*du                                           ! convergence rate
  du = max(du,0.0)                                      ! positive part
  !-----------------------------------------------------------------------------
  ! Conservation of momentum -- first the diagonal part of the stress tensor
  !-----------------------------------------------------------------------------
  if (self%mhd) then
    ! Diagonal Reynolds stress
    fxx = xDN(Ux)**2
    fyy = yDN(Uy)**2
    fzz = zDN(Uz)**2
    pp = (fxx+fyy+fzz)
    ! Diagonal (negative) Maxwell stress (borrowing emf for scratch)
    Ex = yDN(zDN(Bx))**2
    Ey = zDN(xDN(By))**2
    Ez = xDN(yDN(Bz))**2
    ! Magnetic pressure
    pb = 0.5*(Ex+Ey+Ez)
    ! Fast mode speed + velocity
    pa = self%nu(2)*d*du**2
    cs = sqrt((self%gamma*pg+2.0*pb + pa)/d)
    if (io%verbose>2) then
      self%u_max = self%fmaxval(cs)
      print*,'u_max',self%u_max
      self%u_max = self%fmaxval(cs + sqrt(pp))
      print*,'u_max',self%u_max
    end if
    self%u_max = self%fmaxval(cs + sqrt(pp) + du)
    if (io%verbose>2) print*,'u_max',self%u_max
    ! Isotropic pressure part
    pp = pg + pa + pb
    ! Total diagonal stress
    du = dx*self%nu(1)*cs
    fxx = pp + d*(fxx - du*st%ddxDN(Ux)) - Ex
    fyy = pp + d*(fyy - du*st%ddyDN(Uy)) - Ey
    fzz = pp + d*(fzz - du*st%ddzDN(Uz)) - Ex
    Q = d*du*(st%ddxDN(Ux)**2+st%ddyDN(Uy)**2+st%ddzDN(Uz)**2)
  else
    ! Diagonal Reynolds stress
    fxx = xDN(Ux)**2
    fyy = yDN(Uy)**2
    fzz = zDN(Uz)**2
    pp = (fxx+fyy+fzz)
    ! Isotropic pressure part
    pa = self%nu(2)*d*du**2                               ! dimension: pressure
    ! Sound speed + velocity
    cs = sqrt((self%gamma*pg+pa)/d)
    if (io%verbose>2) then
      self%u_max = self%fmaxval(cs)
      print*,'u_max',self%u_max
      self%u_max = self%fmaxval(cs + sqrt(pp))
      print*,'u_max',self%u_max
    end if
    self%u_max = self%fmaxval(cs + sqrt(pp) + du)
    if (io%verbose>2) print*,'u_max',self%u_max
    pp = pg + pa
    ! Total diagondal stress, with density factor
    du = dx*self%nu(1)*cs
    fxx = pp + d*(fxx - du*st%ddxDN(Ux))
    fyy = pp + d*(fyy - du*st%ddyDN(Uy))
    fzz = pp + d*(fzz - du*st%ddzDN(Uz))
    Q = d*du*(st%ddxDN(Ux)**2+st%ddyDN(Uy)**2+st%ddzDN(Uz))
  end if
  !-----------------------------------------------------------------------------
  ! Conservation of mass
  !-----------------------------------------------------------------------------
  fd = -(self%nu(4)*dd*sUP(du))*st%dUP(lnd)             ! diffusive mass flux
  fm = p+fd                                             ! mass flux
  dddt = -st%div_DN(fm)                                 ! density derivative
  self%u_max = max(self%u_max,self%cdtd*dx*self%fmaxval(abs(dddt/d)))
  if (io%verbose>2) print*,'u_max',self%u_max
  !---------------------------------------------------------------------------
  ! Off-diagonal part of Reynolds stress, for now without the density factor.
  ! Note that, because of the symmmetry, we define only three components
  !---------------------------------------------------------------------------
  fxy = yUP(Ux)*xUP(Uy)
  fyz = zUP(Uy)*yUP(Uz)
  fzx = xUP(Uz)*zUP(Ux)
  if (self%id==io%id_debug) print *, 'DPDTZ1', fzx(19:20,19,5)
  !---------------------------------------------------------------------------
  ! Add viscous part -- edge-centered
  !---------------------------------------------------------------------------
  fxy = fxy - self%nu(3)*0.5*xUP1(yUP1(du))*(st%ddxUP(Uy)+st%ddyUP(Ux))
  fyz = fyz - self%nu(3)*0.5*yUP1(zUP1(du))*(st%ddyUP(Uz)+st%ddzUP(Uy))
  fzx = fzx - self%nu(3)*0.5*zUP1(xUP1(du))*(st%ddzUP(Ux)+st%ddxUP(Uz))
  if (self%id==io%id_debug) then
    fzx = st%ddzUP(Uy)
    print *, 'DPDTZ2', fzx(19:20,19,5)
    fzx = self%nu(3)*0.5*st%ddzUP(Uy)
    print *, 'DPDTZ2', fzx(19:20,19,5)
    fzx = self%nu(3)*0.5*yUP1(zUP1(du))*st%ddzUP(Uy)
    print *, 'DPDTZ2', fzx(19:20,19,5)
    fzx = xUP(Uz)*zUP(Ux)
    fzx = fzx - self%nu(3)*0.5*zUP1(xUP1(du))*(st%ddzUP(Ux)+st%ddxUP(Uz))
  if (self%id==io%id_debug) print *, 'DPDTZ2', fzx(19:20,19,5), Uy(19:20,19,4:5), Uz(19:20,18:19,5)
  end if
  !---------------------------------------------------------------------------
  ! Multiply with edge-centered density
  !---------------------------------------------------------------------------
  fxy = exp(yUP(ldx))*fxy
  fyz = exp(zUP(ldy))*fyz
  fzx = exp(xUP(ldz))*fzx
  !---------------------------------------------------------------------------
  ! Off-diagonal Maxwell stress
  !---------------------------------------------------------------------------
  if (self%mhd) then
    fxy = fxy - xUP(zDN(Bx))*yUP(zDN(By))
    fyz = fyz - yUP(xDN(By))*zUP(xDN(Bz))
    fzx = fzx - zUP(yDN(Bz))*xUP(yDN(Bx))
  end if
  !-----------------------------------------------------------------------------
  ! Consistent addition of mass flux related momentum flux, symmetrized. The
  ! diagonal part is analogous, a centered version of \nu \rho U_i \grad\ln\rho,
  ! while the off-diagonal part needs to be properly symmetrized.
  !-----------------------------------------------------------------------------
  if (self%nu(4) > 0.0) then
    pa = d*du*self%nu(4)
    fxx = fxx - pa*st%ddxDN(ldx)*xDN(Ux)
    fyy = fyy - pa*st%ddyDN(ldy)*yDN(Uy)
    fzz = fzz - pa*st%ddzDN(ldz)*zDN(Uz)
    fxy = fxy + (xUP1(fdy)*yUP(Ux)+yUP1(fdx)*xUP(Uy))
    fyz = fyz + (yUP1(fdz)*zUP(Uy)+zUP1(fdy)*yUP(Uz))
    fzx = fzx + (zUP1(fdx)*xUP(Uz)+xUP1(fdz)*zUP(Ux))
    !fxx = fxx + xDN1(fdx)*xDN(Ux)
    !fyy = fyy + yDN1(fdy)*yDN(Uy)
    !fzz = fzz + zDN1(fdz)*zDN(Uz)
    !fxy = fxy + 0.5*(xUP1(fdy)*yUP(Ux)+yUP1(fdx)*xUP(Uy))
    !fyz = fyz + 0.5*(yUP1(fdz)*zUP(Uy)+zUP1(fdy)*yUP(Uz))
    !fzx = fzx + 0.5*(zUP1(fdx)*xUP(Uz)+xUP1(fdz)*zUP(Ux))
    if (self%id==io%id_debug) print *, 'DPDTZ4', fzx(19:20,19,5)
  end if
  !-----------------------------------------------------------------------------
  ! Momentum equations
  !-----------------------------------------------------------------------------
  dpdtx = - st%ddxUP1(fxx) - st%ddyDN1(fxy) - st%ddzDN1(fzx)
  dpdty = - st%ddxDN1(fxy) - st%ddyUP1(fyy) - st%ddzDN1(fyz)
  dpdtz = - st%ddxDN1(fzx) - st%ddyDN1(fyz) - st%ddzUP1(fzz)
  !-----------------------------------------------------------------------------
  ! External force; e.g. point source gravity.  Note that the density factor is
  ! applied here, and should not be included in the function definition
  !-----------------------------------------------------------------------------
  if (associated(self%force%selected)) &
    dpdt = dpdt + dd*self%force%selected (self%time, d, p, Ux, Uy, Uz, self%mesh)
  !-----------------------------------------------------------------------------
  ! Conservation of energy.  The formulation that works best so far is to make
  ! the diffusive energy flux directly proportional to the diffusive mass flux,
  ! plus a contribution proportional to the gradient of specific energy ee.
  !-----------------------------------------------------------------------------
  if (self%gamma == 1.0) then
    dedt = dddt*e/d
  else
    ee = e/d
    fm = fm*sUP(ee) - (self%nu(5)*dd*sUP(du))*st%dUP(ee)
    dedt = dedt - st%div_DN(fm) + Q
  end if
  !-----------------------------------------------------------------------------
  ! Induction equation
  !-----------------------------------------------------------------------------
  if (self%mhd) then
    !---------------------------------------------------------------------------
    ! Electric current, resistive electric field and magnetic dissipation
    ! emf = {Ex,Ey,Ez} is the vector with components Ex, Ey, Ez
    !---------------------------------------------------------------------------
    J = st%curl_DN(B)
    emf = self%nu(6)*edge(du)*J
    if (self%gamma /= 1.0) then
      Q = xDN(Jx*Ex) + yDN(Jy*Ey) + zDN(Jz*Ez)
      dedt = dedt + Q
    end if
    !---------------------------------------------------------------------------
    ! Advective electric field on cell FACES
    !---------------------------------------------------------------------------
    Ex = Ex + xUP(yDN(Uz))*zDN(By) - xUP(zDN(Uy))*yDN(Bz)
    Ey = Ey + yUP(zDN(Ux))*xDN(Bz) - yUP(xDN(Uz))*zDN(Bx)
    Ez = Ez + zUP(xDN(Uy))*yDN(Bx) - zUP(yDN(Ux))*xDN(By)
    dBdt = -st%curl_UP(emf)
  end if
  !-----------------------------------------------------------------------------
  call deallocate_vectors (ld, dd, emf, J, fd, fx, fy, fz, fm)
  call deallocate_scalars (lnd, pa, du, ee, cs, pp, pb, Q)
  if (.not.self%do_pebbles) then
    call deallocate_vectors (U)
    call deallocate_scalars (pg)
  end if
  call trace%end (itimer)
END SUBROUTINE pde

!> =============================================================================
!> Take a complete update step.  This can be overloaded in higher level ud
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

!> =====================================================================
!> Get velocities from momentum
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
!> Put velocities to momentum
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
!> Get thermal energy from total/specific energy
!> =====================================================================
SUBROUTINE e2E_th(self, E_th, it)
  class(mhd_t) :: self
  real, dimension(:,:,:), pointer:: E_th
  integer:: it
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%i_d,it,1), &
            e => self%mem(:,:,:,self%i_e,it,1))
  E_th = e/d
  end associate
END SUBROUTINE e2E_th

!> =====================================================================
!> Put thermal energy to total/specific energy
!> =====================================================================
SUBROUTINE E_th2e(self, E_th, it)
  class(mhd_t) :: self
  real, dimension(:,:,:), pointer:: E_th
  integer:: it
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%i_d,it,1), &
            e => self%mem(:,:,:,self%i_e,it,1))
  e = E_th*d
  end associate
END SUBROUTINE E_th2e

!> =====================================================================
!> Get thermal energy per unit mass from entropy
!> =====================================================================
SUBROUTINE e2S(self, s, it)
  class(mhd_t) :: self
  real, dimension(:,:,:), pointer:: s
  integer:: it
  !.....................................................................
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
  real, dimension(:,:,:), pointer:: s
  integer:: it
  !.....................................................................
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%i_d,it,1), &
            e => self%mem(:,:,:,self%i_s,it,1))
  g1 = self%gamma-1.0
  e = d**self%gamma*exp(s*g1)/g1
  end associate
END SUBROUTINE S2e

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
    self%id, self%iout-1
  inquire(file=filename, exist=exists)
  if (exists) then
    open (io%data_unit, file=trim(filename), form='formatted', status='old', &
      position='append')
    write (io%data_unit,'(8f10.4)') self%nu
    close (io%data_unit)
  end if
  !$omp end critical (out_cr)
END SUBROUTINE output

!> =============================================================================
FUNCTION gas_pressure (self) RESULT (pg)
  class(mhd_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: pg
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_mod::gas_pressure',2)
  if (self%gamma==1.0) then
    pg = self%d*self%csound**2
  else
    pg = self%e*(self%gamma-1.0)
  end if
  call trace%end
END FUNCTION gas_pressure

FUNCTION cell_gas_pressure(self, c) RESULT (pg)
  class(mhd_t) :: self
  integer :: c(3)
  real :: pg
  !---------------------------------------------------------------------
  associate (d => self%mem(c(1),c(2),c(3),self%i_d,self%it,1), &
             e => self%mem(c(1),c(2),c(3),self%i_e,self%it,1))
  if (self%gamma==1.0) then
    pg = self%d(c(1),c(2),c(3))*self%csound**2
  else
    pg = self%e(c(1),c(2),c(3))*(self%gamma-1.0)
  end if
  end associate
END FUNCTION cell_gas_pressure

!> =============================================================================
!> temperature from ideal gas law, taking into account that the temperature is
!> normalised by k_B and mu (so that specific gas constant = 1)
!> =============================================================================
FUNCTION gas_temperature (self) RESULT (tmp)
  class(mhd_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: tmp
  !-----------------------------------------------------------------------------
  associate (d => self%mem(:,:,:,self%i_d,self%it,1))
  tmp = gas_pressure(self)/d
  end associate
END FUNCTION gas_temperature

!> =============================================================================
!> =============================================================================
FUNCTION cell_gas_temperature (self, c) RESULT (tmp)
  class(mhd_t):: self
  integer :: c(3)
  real :: tmp(c(1),c(2),c(3))
  !-----------------------------------------------------------------------------
  associate (d => self%mem(c(1),c(2),c(3),self%i_d,self%it,1))
  tmp = cell_gas_pressure(self,c)/self%d(c(1),c(2),c(3))
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
SUBROUTINE log_pressure (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = log(gas_pressure(self))
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

END MODULE mhd_mod
