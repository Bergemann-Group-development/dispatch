!> =============================================================================
!> =============================================================================
MODULE mhd_mod
  USE io_mod
  USE trace_mod
  USE scalar_mod
  USE vector_mod
  USE stagger_mod
  USE patch_mod
  USE forcing_mod
  implicit none
  private
  type, extends(patch_t), public:: mhd_t
    real, dimension(:,:,:),   pointer:: d, dddt, s, dsdt
    real, dimension(:,:,:,:), pointer:: p, dpdt, B, dBdt
    real:: nu(6)
    logical:: mhd=.true.
    type(force_t):: force
  contains
    procedure:: init
    procedure:: update
    procedure:: log_density => void
    procedure:: log_pressure => void
    procedure:: gas_pressure
    procedure:: velocity_magnitude => void
    procedure:: magnetic_field_magnitude => void
    procedure:: grav_potential => void
  end type
  logical, save:: first_time=.true.
CONTAINS

!> =============================================================================
SUBROUTINE init (self)
  class(mhd_t):: self
  real, save:: nu(6)=1.0, courant=0.2, gamma=1.0
  namelist /mhd_params/ nu, courant, gamma
  !----------------------------------------------------------------------------
  call trace_begin('mhd_t%init')
  !$omp critical (mhd_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, mhd_params)
    if (io%master) write (*, mhd_params)
  end if
  !$omp end critical (mhd_cr)
  self%courant = courant
  self%nu = nu
  self%gamma = gamma
  call self%force%init
  call trace_end
END SUBROUTINE init

!> =============================================================================
SUBROUTINE update (self)
  USE scalar_mod
  class(mhd_t):: self
  !.............................................................................
  real, dimension(:,:,:,:), pointer:: U, E, p, B, dpdt, dBdt
  real, dimension(:,:,:,:), pointer:: fx, fy, fz, fd, ld, dd, fm
  real, dimension(:,:,:),   pointer:: d, s, dddt, dsdt
  real, dimension(:,:,:),   pointer:: du, lnd, ss, cs, pg, pa, pp, pb
  real, dimension(:,:,:),   pointer:: dpdtx, dpdty, dpdtz, dBdtx, dBdty, dBdtz
  real, dimension(:,:,:),   pointer:: Bx, By, Bz, px, py, pz
  real, dimension(:,:,:),   pointer:: Ux, Uy, Uz, Ex, Ey, Ez
  real, dimension(:,:,:),   pointer:: fdx, fdy, fdz, ldx, ldy, ldz, ddx, ddy, ddz
  real, dimension(:,:,:),   pointer:: fxx, fxy, fxz
  real, dimension(:,:,:),   pointer:: fyx, fyy, fyz
  real, dimension(:,:,:),   pointer:: fzx, fzy, fzz
  !.............................................................................
  real:: ds(3), dx
  !-----------------------------------------------------------------------------
  call trace_begin('mhd_t%update')
  self%d => self%mem(:,:,:,  1,self%it,1); self%dddt => self%mem(:,:,:,  1,self%it,2)
  self%s => self%mem(:,:,:,  2,self%it,1); self%dsdt => self%mem(:,:,:,  2,self%it,2)
  self%p => self%mem(:,:,:,3:5,self%it,1); self%dpdt => self%mem(:,:,:,3:5,self%it,2)
  self%B => self%mem(:,:,:,6:8,self%it,1); self%dBdt => self%mem(:,:,:,6:8,self%it,2)
  d=>self%d; dddt=>self%dddt; s=>self%s; dsdt=>self%dsdt;
  p=>self%p; dpdt=>self%dpdt; B=>self%B; dBdt=>self%dBdt
  px=>p(:,:,:,1); py=>p(:,:,:,2); pz=>p(:,:,:,3);
  Bx=>B(:,:,:,1); By=>B(:,:,:,2); Bz=>B(:,:,:,3);
  dpdtx=>dpdt(:,:,:,1); dpdty=>dpdt(:,:,:,2); dpdtz=>dpdt(:,:,:,3);
  dBdtx=>dBdt(:,:,:,1); dBdty=>dBdt(:,:,:,2); dBdtz=>dBdt(:,:,:,3);
  call allocate_vectors (self%gn, ld, dd, U, E, fd, fx, fy, fz, fm)
  call allocate_scalars (self%gn, lnd, pg, pa, du, ss, cs, pp, pb)
  Ux=>U(:,:,:,1); Uy=>U(:,:,:,2); Uz=>U(:,:,:,3);
  Ex=>E(:,:,:,1); Ey=>E(:,:,:,2); Ez=>E(:,:,:,3);
  fxx=>fx(:,:,:,1); fxy=>fx(:,:,:,2); fxz=>fx(:,:,:,3)
  fyx=>fy(:,:,:,1); fyy=>fy(:,:,:,2); fyz=>fy(:,:,:,3)
  fzx=>fz(:,:,:,1); fzy=>fz(:,:,:,2); fzz=>fz(:,:,:,3)
  fdx=>fd(:,:,:,1); fdy=>fd(:,:,:,2); fdz=>fd(:,:,:,3)
  ldx=>ld(:,:,:,1); ldy=>ld(:,:,:,2); ldz=>ld(:,:,:,3)
  ddx=>dd(:,:,:,1); ddy=>dd(:,:,:,2); ddz=>dd(:,:,:,3)
  call self%init
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
  pg = self%gas_pressure()
  !-----------------------------------------------------------------------------
  ! Divergence of velocity, and artificial pressure
  !-----------------------------------------------------------------------------
  ds = self%mesh%d
  dx = minval(ds)                                       ! smallest {dx,dy,dz}
  du = -div(ds,U)                                       ! convergence rate
  du = dx*max(du,1e-20)                                 ! dimension: U
  pa = self%nu(3)*d*du**2                               ! dimension: pressure
  !-----------------------------------------------------------------------------
  ! Conservation of momentum -- first the diagonal part of the stress tensor
  !-----------------------------------------------------------------------------
  if (self%mhd) then
    ! Diagonal Reynolds stress
    fxx = xup(Ux)**2
    fyy = yup(Uy)**2
    fzz = zup(Uz)**2
    pp = (fxx+fyy+fzz)
    ! Diagonal (negative) Maxwell stress
    fxx = xup(Bx)**2
    fyy = yup(By)**2
    fzz = zup(Bz)**2
    ! Magnetic pressure
    pb = 0.5*(fxx+fyy+fzz)
    ! Fast mode speed + velocity
    cs = sqrt((self%gamma*pg+2.0*pb)/d)
    self%u_max = maxval(cs + sqrt(pp))
    ! Isotropic pressure part
    pp = d*pp + pg + pa + pb
    ! Total diagondal stress
    fxx = pp - fxx - d*du*ddxup(ds,Ux)
    fyy = pp - fyy - d*du*ddyup(ds,Uy)
    fzz = pp - fzz - d*du*ddzup(ds,Uz)
  else
    ! Diagonal Reynolds stress
    fxx = xup(Ux)**2
    fyy = yup(Uy)**2
    fzz = zup(Uz)**2
    pp = (fxx+fyy+fzz)
    ! Sound speed + velocity
    cs = sqrt(self%gamma*pg/d)
    self%u_max = maxval(cs + sqrt(pp))
    ! Isotropic pressure part
    pp = d*pp + pg + pa
    ! Total diagondal stress, with density factor
    fxx = pp - d*du*ddxup(ds,Ux)
    fyy = pp - d*du*ddyup(ds,Uy)
    fzz = pp - d*du*ddzup(ds,Uz)
  end if
  du = (self%nu(1)*cs + self%nu(2)*du)*dx               ! dimension: U*ds
  !-----------------------------------------------------------------------------
  ! Conservation of mass
  !-----------------------------------------------------------------------------
  fd = -(self%nu(5)*dd*down(du))*ddown(ds,lnd)          ! diffusive mass flux
  fm = p+fd                                             ! mass flux
  dddt = -div(ds,fm)                                    ! density derivative
  !-----------------------------------------------------------------------------
  ! Conservation of entropy.  The formulation that works best so far is to make
  ! the diffusive entropy flux directly proportional to the diffusive mass flux,
  ! plus a contribution proportional to the gradient of specific entropy ss.
  !-----------------------------------------------------------------------------
  if (self%gamma==1.0) then
    dsdt = 0.0
  else
    ss = s/d
    fm = fm*down(ss) - (self%nu(6)*dd*down(du))*ddown(ds,ss)
    dsdt = -div(ds,fm)
  end if
  !---------------------------------------------------------------------------
  ! Off-diagonal part of Reynolds stress, for now without the density factor.
  ! Note that, because of the symmmetry, we define only three components
  !---------------------------------------------------------------------------
  fxy = ydn(Ux)*xdn(Uy)
  fxz = zdn(Ux)*xdn(Uz)
  fyz = zdn(Uy)*ydn(Uz)
  !---------------------------------------------------------------------------
  ! Add viscous part -- edge-centered
  !---------------------------------------------------------------------------
  fxy = fxy - self%nu(4)*xdn1(ydn1(du))*(ddxdn(ds,Uy)+ddydn(ds,Ux))
  fxz = fxz - self%nu(4)*ydn1(zdn1(du))*(ddxdn(ds,Uz)+ddzdn(ds,Ux))
  fyz = fyz - self%nu(4)*zdn1(xdn1(du))*(ddydn(ds,Uz)+ddzdn(ds,Uy))
  !---------------------------------------------------------------------------
  ! Multiply with edge-centered density
  !---------------------------------------------------------------------------
  fxy = exp(ydn(ldx))*fxy
  fxz = exp(zdn(ldx))*fxz
  fyz = exp(zdn(ldy))*fyz
  !---------------------------------------------------------------------------
  ! Off-diagonal Maxwell stress
  !---------------------------------------------------------------------------
  if (self%mhd) then
    fxy = fxy - ydn(Bx)*xdn(By)
    fxz = fxz - zdn(Bx)*xdn(Bz)
    fyz = fyz - zdn(By)*ydn(Bz)
  end if
  !-----------------------------------------------------------------------------
  ! Consistent addition of mass flux related momentum flux, symmetrized. The
  ! diagonal part is analogous, a centered version of \nu \rho U_i \grad\ln\rho,
  ! while the off-diagonal part needs to be properly symmetrized.
  !-----------------------------------------------------------------------------
  pa = d*du*self%nu(5)
  fxx = fxx - pa*ddxup(ds,ldx)*xup(Ux)
  fyy = fyy - pa*ddyup(ds,ldy)*yup(Uy)
  fzz = fzz - pa*ddzup(ds,ldz)*zup(Uz)
  fxy = fxy + (xdn1(fdy)*ydn(Ux)+ydn1(fdx)*xdn(Uy))
  fxz = fxz + (xdn1(fdz)*zdn(Ux)+zdn1(fdx)*xdn(Uz))
  fyz = fyz + (zdn1(fdy)*ydn(Uz)+ydn1(fdz)*zdn(Uy))
  !-----------------------------------------------------------------------------
  ! Momentum equations
  !-----------------------------------------------------------------------------
  dpdtx = - ddxdn1(ds,fxx) - ddyup1(ds,fxy) - ddzup1(ds,fxz)
  dpdty = - ddxup1(ds,fxy) - ddydn1(ds,fyy) - ddzup1(ds,fyz)
  dpdtz = - ddxup1(ds,fxz) - ddyup1(ds,fyz) - ddzdn1(ds,fzz)
  if (associated(self%force%selected)) &
    dpdt = dpdt + dd*self%force%selected (self%time, self%mesh)
  if (self%mhd) then
    !---------------------------------------------------------------------------
    ! Advective electric field on cell edges
    !---------------------------------------------------------------------------
    Ex = ydn(Uz)*zdn(By) - zdn(Uy)*ydn(Bz)
    Ey = zdn(Ux)*xdn(Bz) - xdn(Uz)*zdn(Bx)
    Ez = xdn(Uy)*ydn(Bx) - ydn(Ux)*xdn(By)
    !---------------------------------------------------------------------------
    ! Induction equation, with resistive term
    !---------------------------------------------------------------------------
    E = E + self%nu(6)*edge(du)*curl_dn(ds,B)
    dBdt = -curl_up(ds,E)
  end if
  !-----------------------------------------------------------------------------
  call deallocate_vectors (ld, dd, U, E, fd, fx, fy, fz, fm)
  call deallocate_scalars (lnd, pg, pa, du, ss, cs, pp, pb)
  call trace_end
END SUBROUTINE update

!> =============================================================================
FUNCTION gas_pressure (self) RESULT (pg)
  class(mhd_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: pg
  !-----------------------------------------------------------------------------
  if (self%gamma==1.0) then
    pg = self%d
  else
    pg = self%d**self%gamma*exp(self%s/self%d*(self%gamma-1))
  end if
END FUNCTION gas_pressure

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
