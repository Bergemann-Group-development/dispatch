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
    procedure:: update
    procedure:: output
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
  call self%patch_t%init
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
  call self%force%init
  call trace%end
END SUBROUTINE init

!> =============================================================================
SUBROUTINE update (self)
  USE scalar_mod
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
  real, dimension(:,:,:),   pointer:: fxx, fxy, fxz
  real, dimension(:,:,:),   pointer:: fyx, fyy, fyz
  real, dimension(:,:,:),   pointer:: fzx, fzy, fzz
  !.............................................................................
  real:: ds(3), dx
  integer, save:: itimer=0
  type(stagger_t):: st
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%update',itimer=itimer)
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
  fxx=>fx(:,:,:,1); fxy=>fx(:,:,:,2); fxz=>fx(:,:,:,3)
  fyx=>fy(:,:,:,1); fyy=>fy(:,:,:,2); fyz=>fy(:,:,:,3)
  fzx=>fz(:,:,:,1); fzy=>fz(:,:,:,2); fzz=>fz(:,:,:,3)
  fdx=>fd(:,:,:,1); fdy=>fd(:,:,:,2); fdz=>fd(:,:,:,3)
  ldx=>ld(:,:,:,1); ldy=>ld(:,:,:,2); ldz=>ld(:,:,:,3)
  ddx=>dd(:,:,:,1); ddy=>dd(:,:,:,2); ddz=>dd(:,:,:,3)
  !-----------------------------------------------------------------------------
  ! Compute velocity U=momentum/density, leaving down shifted density and log
  ! density for reuse below
  !-----------------------------------------------------------------------------
  lnd = log(d)
  ld = down(lnd)
  dd = exp(ld)
  U = p/dd
  !-----------------------------------------------------------------------------
  ! Gas pressure
  !-----------------------------------------------------------------------------
  ee = e/d
!  call eos%lookup (shape(ee), ee=ee, lnd=lnd, pg=pg)
  if (self%gamma==1.0) then
    pg = d*self%csound**2
  else
    pg = e*(self%gamma-1.0)
  end if
  !-----------------------------------------------------------------------------
  ! Divergence of velocity, and artificial pressure
  !-----------------------------------------------------------------------------
  ds = self%mesh%d
  dx = minval(ds)                                       ! smallest {dx,dy,dz}
  du = st%div(U)                                        ! save div(U) for now
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
    fxx = xup(Ux)**2
    fyy = yup(Uy)**2
    fzz = zup(Uz)**2
    pp = (fxx+fyy+fzz)
    ! Diagonal (negative) Maxwell stress (borrowing emf for scratch)
    Ex = xup(Bx)**2
    Ey = yup(By)**2
    Ez = zup(Bz)**2
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
    fxx = pp + d*(fxx - du*st%ddxup(Ux)) - Ex
    fyy = pp + d*(fyy - du*st%ddyup(Uy)) - Ey
    fzz = pp + d*(fzz - du*st%ddzup(Uz)) - Ex
    Q = d*du*(st%ddxup(Ux)**2+st%ddyup(Uy)**2+st%ddzup(Uz)**2)
  else
    ! Diagonal Reynolds stress
    fxx = xup(Ux)**2
    fyy = yup(Uy)**2
    fzz = zup(Uz)**2
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
    fxx = pp + d*(fxx - du*st%ddxup(Ux))
    fyy = pp + d*(fyy - du*st%ddyup(Uy))
    fzz = pp + d*(fzz - du*st%ddzup(Uz))
    Q = d*du*(st%ddxup(Ux)**2+st%ddyup(Uy)**2+st%ddzup(Uz))
  end if
  !-----------------------------------------------------------------------------
  ! Conservation of mass
  !-----------------------------------------------------------------------------
  fd = -(self%nu(4)*dd*down(du))*st%ddown(lnd)          ! diffusive mass flux
  fm = p+fd                                             ! mass flux
  dddt = -st%div(fm)                                    ! density derivative
  self%u_max = max(self%u_max,self%cdtd*dx*self%fmaxval(abs(dddt/d)))
  if (io%verbose>2) print*,'u_max',self%u_max
  !---------------------------------------------------------------------------
  ! Off-diagonal part of Reynolds stress, for now without the density factor.
  ! Note that, because of the symmmetry, we define only three components
  !---------------------------------------------------------------------------
  fxy = ydn(Ux)*xdn(Uy)
  fyz = zdn(Uy)*ydn(Uz)
  fzx = xdn(Uz)*zdn(Ux)
  !---------------------------------------------------------------------------
  ! Add viscous part -- edge-centered
  !---------------------------------------------------------------------------
  fxy = fxy - self%nu(3)*0.5*xdn1(ydn1(du))*(st%ddxdn(Uy)+st%ddydn(Ux))
  fyz = fyz - self%nu(3)*0.5*ydn1(zdn1(du))*(st%ddydn(Uz)+st%ddzdn(Uy))
  fzx = fzx - self%nu(3)*0.5*zdn1(xdn1(du))*(st%ddzdn(Ux)+st%ddxdn(Uz))
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
  if (self%nu(4) > 0.0) then
    pa = d*du*self%nu(4)
    fxx = fxx - pa*st%ddxup(ldx)*xup(Ux)
    fyy = fyy - pa*st%ddyup(ldy)*yup(Uy)
    fzz = fzz - pa*st%ddzup(ldz)*zup(Uz)
    fxy = fxy + (xdn1(fdy)*ydn(Ux)+ydn1(fdx)*xdn(Uy))
    fyz = fyz + (ydn1(fdz)*zdn(Uy)+zdn1(fdy)*ydn(Uz))
    fzx = fzx + (zdn1(fdx)*xdn(Uz)+xdn1(fdz)*zdn(Ux))
    !fxx = fxx + xup1(fdx)*xup(Ux)
    !fyy = fyy + yup1(fdy)*yup(Uy)
    !fzz = fzz + zup1(fdz)*zup(Uz)
    !fxy = fxy + 0.5*(xdn1(fdy)*ydn(Ux)+ydn1(fdx)*xdn(Uy))
    !fyz = fyz + 0.5*(ydn1(fdz)*zdn(Uy)+zdn1(fdy)*ydn(Uz))
    !fzx = fzx + 0.5*(zdn1(fdx)*xdn(Uz)+xdn1(fdz)*zdn(Ux))
  end if
  !-----------------------------------------------------------------------------
  ! Momentum equations
  !-----------------------------------------------------------------------------
  dpdtx = - st%ddxdn1(fxx) - st%ddyup1(fxy) - st%ddzup1(fzx)
  dpdty = - st%ddxup1(fxy) - st%ddydn1(fyy) - st%ddzup1(fyz)
  dpdtz = - st%ddxup1(fzx) - st%ddyup1(fyz) - st%ddzdn1(fzz)
  !-----------------------------------------------------------------------------
  ! External force; e.g. point source gravity.  Note that the density factor is
  ! applied here, and should not be included in the function definition
  !-----------------------------------------------------------------------------
  if (associated(self%force%selected)) &
    dpdt = dpdt + dd*self%force%selected (self%time, Ux, Uy, Uz, self%mesh)
  !-----------------------------------------------------------------------------
  ! Conservation of energy.  The formulation that works best so far is to make
  ! the diffusive energy flux directly proportional to the diffusive mass flux,
  ! plus a contribution proportional to the gradient of specific energy ee.
  !-----------------------------------------------------------------------------
  if (self%gamma /= 1.0) then
    ee = e/d
    fm = fm*down(ee) - (self%nu(5)*dd*down(du))*st%ddown(ee)
    dedt = dedt - st%div(fm) + Q
  end if
  !-----------------------------------------------------------------------------
  ! Induction equation
  !-----------------------------------------------------------------------------
  if (self%mhd) then
    !---------------------------------------------------------------------------
    ! Electric current, resistive electric field and magnetic dissipation
    ! emf = {Ex,Ey,Ez} is the vector with components Ex, Ey, Ez
    !---------------------------------------------------------------------------
    J = st%curl_dn(B)
    emf = self%nu(6)*edge(du)*J
    if (self%gamma /= 1.0) then
      Q = yup(zup(Jx*Ex)) + zup(xup(Jy*Ey)) + xup(yup(Jz*Ez))
      dedt = dedt + Q
    end if
    !---------------------------------------------------------------------------
    ! Advective electric field on cell edges
    !---------------------------------------------------------------------------
    Ex = Ex + ydn(Uz)*zdn(By) - zdn(Uy)*ydn(Bz)
    Ey = Ey + zdn(Ux)*xdn(Bz) - xdn(Uz)*zdn(Bx)
    Ez = Ez + xdn(Uy)*ydn(Bx) - ydn(Ux)*xdn(By)
    dBdt = -st%curl_up(emf)
  end if
  !-----------------------------------------------------------------------------
  call deallocate_vectors (ld, dd, emf, J, fd, fx, fy, fz, fm)
  call deallocate_scalars (lnd, pa, du, ee, cs, pp, pb, Q)
  if (.not.self%do_pebbles) then
    call deallocate_vectors (U)
    call deallocate_scalars (pg)
  end if
  call self%courant_condition
  if (io%verbose>2) &
    call self%check_density('before timestep')
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call self%check_density(' after timestep')
  call trace%end (itimer)
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
  !$omp critical (output_cr)
  write (filename,'(a,i4.4,"_",i4.4,".txt")') trim(io%outputname), &
    self%id, self%iout-1
  inquire(file=filename, exist=exists)
  if (exists) then
    open (io%data_unit, file=trim(filename), form='formatted', status='old', &
      position='append')
    write (io%data_unit,'(8f10.4)') self%nu
    close (io%data_unit)
  end if
  !$omp end critical (output_cr)
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
