!> =============================================================================
!> =============================================================================
MODULE mhd_mod
  USE io_mod
  USE bits_mod
  USE omp_mod
  USE mpi_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  USE scalar_mod
  USE vector_mod
  USE stagger_mod
  USE patch_mod
  USE force_mod
  USE download_mod
  USE timestep_mod
  implicit none
  private
  type, extends(patch_t), public:: mhd_t
    real, dimension(:,:,:),   pointer:: d, dddt, e, dedt, lnd
    real, dimension(:,:,:,:), pointer:: p, dpdt, B, dBdt
    real, dimension(:,:,:,:), pointer:: U, ld, dd
    real:: nu(6), csound, cdtd=0.3
    logical:: mhd=.true.
    type(force_t):: force
  contains
    procedure:: init
    procedure:: p2U
    procedure:: U2p
    procedure:: update
    procedure:: dnload
    procedure:: pde
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
CONTAINS

!> =============================================================================
!> The construction with a double test on first_time avoids entering a critical
!> region once the variable has been set.
!> =============================================================================
SUBROUTINE init (self)
  class(mhd_t):: self
  real, save:: nu(6)=[0.1,1.0,0.0,0.5,0.5,0.5], gamma=1.0, csound=1.0
  integer:: i
  logical, save:: first_time=.true.
  namelist /mhd_params/ nu, gamma, csound
  !----------------------------------------------------------------------------
  call trace%begin('mhd_t%init')
  call self%patch_t%init
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, mhd_params)
    if (io%master) write (*, mhd_params)
  end if
  !$omp end critical (input_cr)
  self%nu = nu
  self%gamma = gamma
  self%csound = csound
  self%mesh(1)%h(1:2) = 0.0
  self%mesh(2)%h(1:2) = 0.0
  self%mesh(3)%h(1:2) = 0.0
  !----------------------------------------------------------------------------
  ! Stagger offsets for (up-) face centered mass flux components
  !----------------------------------------------------------------------------
  self%mesh(1)%h(3) = +0.5
  self%mesh(2)%h(4) = +0.5
  self%mesh(3)%h(5) = +0.5
  !----------------------------------------------------------------------------
  ! Stagger offsets for (up-up) edge centered magnetic field components
  !----------------------------------------------------------------------------
  self%mesh(2)%h(6) = +0.5
  self%mesh(3)%h(6) = +0.5
  self%mesh(1)%h(7) = +0.5
  self%mesh(3)%h(7) = +0.5
  self%mesh(1)%h(8) = +0.5
  self%mesh(2)%h(8) = +0.5
  self%kind='stagger2e_pic'
  call self%force%init
  call timestep%init
  ! fool-proofing
  if (timestep%time_order <= 0) call mpi%abort("The Stagger solvers require `time_order` > 0. Abort!")
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> =============================================================================
SUBROUTINE dnload (self)
  class(mhd_t):: self
  call trace%begin ('experiment_mod::dnload')
  call download%download_link (self%link, all_cells=.false.)
  call trace%end
END SUBROUTINE

!> =============================================================================
!> Compute velocity
!> =============================================================================
SUBROUTINE p2U (self)
  class(mhd_t):: self
  real, dimension(:,:,:), pointer:: lnd
  real, dimension(:,:,:,:), pointer:: ld, dd
  !.............................................................................
  call trace%begin ('experiment_mod::p2U')
  call allocate_scalars (self%gn, lnd)
  call allocate_vectors (self%gn, ld, self%dd, self%U)
  self%d => self%mem(:,:,:,self%i_d,self%it,1)
  self%p => self%mem(:,:,:,self%i_px:self%i_pz,self%it,1)
  lnd = log(self%d)
  ld  = vup(lnd)
  self%dd  = exp(ld)
  self%U   = self%p/self%dd
  call deallocate_scalars (lnd)
  call deallocate_vectors (ld)
  call trace%end
END SUBROUTINE p2U

!> =============================================================================
!> Back to momentum, assuming previous call to p2U
!> =============================================================================
SUBROUTINE U2p (self)
  class(mhd_t):: self
  real, dimension(:,:,:), pointer:: lnd
  real, dimension(:,:,:,:), pointer:: ld, dd
  !.............................................................................
  call trace%begin ('experiment_mod::U2p')
  self%p = self%U*self%dd
  call deallocate_vectors (self%dd, self%U)
  call trace%end
END SUBROUTINE U2p

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
!> =============================================================================
SUBROUTINE pde (self)
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
  real, dimension(:,:,:),   pointer:: fxx, fxy
  real, dimension(:,:,:),   pointer::      fyy, fyz
  real, dimension(:,:,:),   pointer:: fzx,      fzz
  !.............................................................................
  real:: ds(3), dx
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('mhd_t%pde',itimer=itimer)
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
  ! Compute velocity U=momentum/density, leaving up shifted density and log
  ! density for reuse below
  !-----------------------------------------------------------------------------
  lnd = log(d)
  ld  = vup(lnd)
  dd  = exp(ld)
  U   = p/dd
  !-----------------------------------------------------------------------------
  ! Gas pressure
  !-----------------------------------------------------------------------------
  ee = e/d
  if (self%gamma==1.0) then
    pg = d*self%csound**2
  else
    pg = e*(self%gamma-1.)
  end if
  !-----------------------------------------------------------------------------
  ! Divergence of velocity, and artificial pressure
  !-----------------------------------------------------------------------------
  ds = self%mesh%d
  dx = minval(ds)                                       ! smallest {dx,dy,dz}
  du = div_dn(ds,U)                                     ! save div(U) for now
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
    fxx = xdn(Ux)**2
    fyy = ydn(Uy)**2
    fzz = zdn(Uz)**2
    pp = (fxx+fyy+fzz)
    ! Diagonal (negative) Maxwell stress (borrowing emf for scratch)
    Ex = ydn(zdn(Bx))**2
    Ey = zdn(xdn(By))**2
    Ez = xdn(ydn(Bz))**2
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
    fxx = pp + d*(fxx - du*ddxdn(ds,Ux)) - Ex
    fyy = pp + d*(fyy - du*ddydn(ds,Uy)) - Ey
    fzz = pp + d*(fzz - du*ddzdn(ds,Uz)) - Ez
    Q = d*du*(ddxdn(ds,Ux)**2+ddydn(ds,Uy)**2+ddzdn(ds,Uz)**2)
  else
    ! Diagonal Reynolds stress
    fxx = xdn(Ux)**2
    fyy = ydn(Uy)**2
    fzz = zdn(Uz)**2
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
    fxx = pp + d*(fxx - du*ddxdn(ds,Ux))
    fyy = pp + d*(fyy - du*ddydn(ds,Uy))
    fzz = pp + d*(fzz - du*ddzdn(ds,Uz))
    Q = d*du*(ddxdn(ds,Ux)**2+ddydn(ds,Uy)**2+ddzdn(ds,Uz))
  end if
  !-----------------------------------------------------------------------------
  ! Conservation of mass
  !-----------------------------------------------------------------------------
  fd = -(self%nu(4)*dd*vup(du))*dupv(ds,lnd)            ! diffusive mass flux
  fm = p+fd                                             ! mass flux
  dddt = -div_dn(ds,fm)                                 ! density derivative
  self%u_max = max(self%u_max,self%cdtd*dx*self%fmaxval(abs(dddt/d)))
  if (io%verbose>2) print*,'u_max',self%u_max
  !---------------------------------------------------------------------------
  ! Off-diagonal part of Reynolds stress, for now without the density factor.
  ! Note that, because of the symmmetry, we define only three components
  !---------------------------------------------------------------------------
  fxy = yup(Ux)*xup(Uy)
  fyz = zup(Uy)*yup(Uz)
  fzx = xup(Uz)*zup(Ux)
  !---------------------------------------------------------------------------
  ! Add viscous part -- edge-centered
  !---------------------------------------------------------------------------
  fxy = fxy - self%nu(3)*0.5*xup1(yup1(du))*(ddxup(ds,Uy)+ddyup(ds,Ux))
  fyz = fyz - self%nu(3)*0.5*yup1(zup1(du))*(ddyup(ds,Uz)+ddzup(ds,Uy))
  fzx = fzx - self%nu(3)*0.5*zup1(xup1(du))*(ddzup(ds,Ux)+ddxup(ds,Uz))
  !---------------------------------------------------------------------------
  ! Multiply with edge-centered density
  !---------------------------------------------------------------------------
  fxy = exp(yup(ldx))*fxy
  fyz = exp(zup(ldy))*fyz
  fzx = exp(xup(ldz))*fzx
  !---------------------------------------------------------------------------
  ! Off-diagonal Maxwell stress
  !---------------------------------------------------------------------------
  if (self%mhd) then
    fxy = fxy - xup(zdn(Bx))*yup(zdn(By))
    fyz = fyz - yup(xdn(By))*zup(xdn(Bz))
    fzx = fzx - zup(ydn(Bz))*xup(ydn(Bx))
  end if
  !-----------------------------------------------------------------------------
  ! Consistent addition of mass flux related momentum flux, symmetrized. The
  ! diagonal part is analogous, a centered version of \nu \rho U_i \grad\ln\rho,
  ! while the off-diagonal part needs to be properly symmetrized.
  !-----------------------------------------------------------------------------
  if (self%nu(4) > 0.0) then
    pa = d*du*self%nu(4)
    fxx = fxx - pa*ddxdn(ds,ldx)*xdn(Ux)
    fyy = fyy - pa*ddydn(ds,ldy)*ydn(Uy)
    fzz = fzz - pa*ddzdn(ds,ldz)*zdn(Uz)
    fxy = fxy + (xup1(fdy)*yup(Ux)+yup1(fdx)*xup(Uy))
    fyz = fyz + (yup1(fdz)*zup(Uy)+zup1(fdy)*yup(Uz))
    fzx = fzx + (zup1(fdx)*xup(Uz)+xup1(fdz)*zup(Ux))
  end if
  !-----------------------------------------------------------------------------
  ! Momentum equations
  !-----------------------------------------------------------------------------
  dpdtx = - ddxup(ds,fxx) - ddydn(ds,fxy) - ddzdn(ds,fzx)
  dpdty = - ddxdn(ds,fxy) - ddyup(ds,fyy) - ddzdn(ds,fyz)
  dpdtz = - ddxdn(ds,fzx) - ddydn(ds,fyz) - ddzup(ds,fzz)
  !-----------------------------------------------------------------------------
  ! External force; e.g. point source gravity.  Note that the density factor is
  ! applied here, and should not be included in the function definition
  !-----------------------------------------------------------------------------
  if (associated(self%force%selected)) then
    call self%force%selected (self%time, self%mesh, dd, Ux, Uy, Uz, dpdt)
  end if
  !-----------------------------------------------------------------------------
  ! Conservation of energy.  The formulation that works best so far is to make
  ! the diffusive energy flux directly proportional to the diffusive mass flux,
  ! plus a contribution proportional to the gradient of specific energy ee.
  !-----------------------------------------------------------------------------
  if (self%gamma == 1.0) then
    dedt = dddt*e/d
  else
    ee = e/d
    fm = fm*vup(ee) - (self%nu(5)*dd*vup(du))*dupv(ds,ee)
    dedt = dedt - div_dn(ds,fm) + Q
  end if
  !-----------------------------------------------------------------------------
  ! Induction equation
  !-----------------------------------------------------------------------------
  if (self%mhd) then
    !---------------------------------------------------------------------------
    ! Electric current, resistive electric field and magnetic dissipation
    ! emf = {Ex,Ey,Ez} is the vector with components Ex, Ey, Ez
    !---------------------------------------------------------------------------
    J = curl_dn(ds,B)
    emf = self%nu(6)*edge(du)*J
    if (self%gamma /= 1.0) then
      Q = xdn(Jx*Ex) + ydn(Jy*Ey) + zdn(Jz*Ez)
      dedt = dedt + Q
    end if
    !---------------------------------------------------------------------------
    ! Advective electric field on cell FACES
    !---------------------------------------------------------------------------
    Ex = Ex + xup(ydn(Uz))*zdn(By) - xup(zdn(Uy))*ydn(Bz)
    Ey = Ey + yup(zdn(Ux))*xdn(Bz) - yup(xdn(Uz))*zdn(Bx)
    Ez = Ez + zup(xdn(Uy))*ydn(Bx) - zup(ydn(Ux))*xdn(By)
    dBdt = -curl_up(ds,emf)
  end if
  !-----------------------------------------------------------------------------
  call deallocate_vectors (ld, dd, emf, J, fd, fx, fy, fz, fm)
  call deallocate_scalars (lnd, pa, du, ee, cs, pp, pb, Q)
  if (.not.self%do_pebbles) then
    call deallocate_vectors (U)
    call deallocate_scalars (pg)
  end if
  call trace_end (itimer)
END SUBROUTINE pde

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
  call self%patch_t%output (experiment_name)
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
  if (self%gamma==1.0) then
    pg = self%d*self%csound**2
  else
    pg = self%e*(self%gamma-1.0)
  end if
END FUNCTION gas_pressure

FUNCTION cell_gas_pressure(self, c) RESULT (pg)
  class(mhd_t) :: self
  integer :: c(3)
  real :: pg(c(1),c(2),c(3))
  !---------------------------------------------------------------------
  associate (d => self%mem(c(1),c(2),c(3), 1,self%it,1), &
             e => self%mem(c(1),c(2),c(3), 2,self%it,1))
  if (self%gamma==1.0) then
    pg = d*self%csound**2
  else
    pg = self%e*(self%gamma-1.0)
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
  associate (d => self%mem(:,:,:,  1,self%it,1))
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
  v = norm(self%p/exp(vup(log(self%d))))
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
