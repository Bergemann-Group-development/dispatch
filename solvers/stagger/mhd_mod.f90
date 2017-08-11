!> =============================================================================
!> =============================================================================
MODULE mhd_mod
  USE io_mod
  USE mpi_mod
  USE trace_mod
  USE timestep_mod
  USE scalar_mod
  USE vector_mod
  USE stagger_mod
  USE patch_mod
  USE force_mod
  USE download_mod
  implicit none
  private
  type, extends(patch_t), public:: mhd_t
    type(scalar_t):: d, dddt, s, dsdt
    type(vector_t):: p, dpdt, B, dBdt
    type(force_t):: force
    real:: nu(6)
    logical:: mhd=.true.
  contains
    procedure:: init
    procedure:: update
    procedure:: log_density
    procedure:: log_pressure
    procedure:: gas_pressure
    procedure:: velocity_magnitude
    procedure:: magnetic_field_magnitude
    procedure:: grav_potential => void
  end type
  logical, save:: first_time=.true.
CONTAINS

!> =============================================================================
SUBROUTINE init (self)
  class(mhd_t):: self
  real:: nu(6)=1.0, courant=0.2, gamma=1.0
  namelist /mhd_parameters/ nu, courant
  namelist /mhd_params/ nu, courant, gamma
  !----------------------------------------------------------------------------
  call trace_begin('mhd_t%init')
  !$omp critical (mhd_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, mhd_params)
    if (io%master) write (*, mhd_params)
    self%courant = courant
    self%nu = nu
    self%gamma = gamma
    call timestep%init
    ! fool-proofing
    if (timestep%time_order <= 0) call mpi%abort("The Stagger solvers require `time_order` > 0. Abort!")
  end if
  !$omp end critical (mhd_cr)
  call trace_end
END SUBROUTINE init

!> =============================================================================
SUBROUTINE update (self)
  class(mhd_t):: self
  !.............................................................................
  type(vector_t):: U, E, J                              ! velocity, electric field
  type(vector_t):: fx, fy, fz, fd, ld, dd, fm           ! extra vector temporaries
  type(scalar_t):: du, lnd, ss, cs, pg, pa, pp, pb      ! extra scalar temporaries
  real, dimension(:,:,:),   pointer:: s
  real, dimension(:,:,:,:), pointer:: v
  !.............................................................................
  integer:: n(3)
  real:: dx, ds(3)
  !-----------------------------------------------------------------------------
  call trace_begin('mhd_t%update')
  if (first_time) call self%init
  !----------------------------------------------------------------------------
  ! Fill guard zones and checks whether to do output
  !----------------------------------------------------------------------------
  call download%download_link (self%link, all_cells=.true.)
  call self%output
  !----------------------------------------------------------------------------
  ! Map variables to memory
  !----------------------------------------------------------------------------
  n = self%mesh%gn
  s=>self%mem(:,:,:,  1,self%it,1); call self%d%scalar_map    (s, n, self%ds)
  s=>self%mem(:,:,:,  2,self%it,1); call self%s%scalar_map    (s, n, self%ds)
  v=>self%mem(:,:,:,3:5,self%it,1); call self%p%vector_map    (v, n, self%ds)
  v=>self%mem(:,:,:,6:8,self%it,1); call self%B%vector_map    (v, n, self%ds)
  s=>self%mem(:,:,:,  1,self%it,2); call self%dddt%scalar_map (s, n, self%ds)
  s=>self%mem(:,:,:,  2,self%it,2); call self%dsdt%scalar_map (s, n, self%ds)
  v=>self%mem(:,:,:,3:5,self%it,2); call self%dpdt%vector_map (v, n, self%ds)
  v=>self%mem(:,:,:,6:8,self%it,2); call self%dBdt%vector_map (v, n, self%ds)
  call init_vectors (self%mesh%gn, self%ds, ld, dd, U, E, fd, fm, fx, fy, fz)
  call init_scalars (self%mesh%gn, self%ds, lnd, pg, pa, du, ss, cs, pp, pb)
  !-----------------------------------------------------------------------------
  ! Compute velocity U=momentum/density, leaving down shifted density and log
  ! density for reuse below
  !-----------------------------------------------------------------------------
  associate(d=>self%d,dddt=>self%dddt,s=>self%s,dsdt=>self%dsdt,p=>self%p,dpdt=>self%dpdt,B=>self%B,dBdt=>self%dBdt)
  lnd = log(d)
  ld = down(lnd)
  dd = exp(ld)
  U = p/dd
  !-----------------------------------------------------------------------------
  ! Gas pressure
  !-----------------------------------------------------------------------------
  pg = self%gas_pressure()
  !-----------------------------------------------------------------------------
  ! Diffusion coefficient; the result 'nu' is cell centered, and may be staggered
  ! down without increasing the number of ghost zones needed
  !-----------------------------------------------------------------------------
  ds = self%mesh%d
  dx = minval(ds)                                       ! smallest {dx,dy,dz}
  du = -dx*div(U)                                    ! convergence rate
  !-----------------------------------------------------------------------------
  ! Conservation of momentum -- first the diagonal part of the stress tensor
  !-----------------------------------------------------------------------------
  if (self%mhd) then
    ! Diagonal Reynolds stress
    fx%x = xup(U%x)**2
    fy%y = yup(U%y)**2
    fz%z = zup(U%z)**2
    pp = (fx%x+fy%y+fz%z)
    ! Diagonal (negative) Maxwell stress
    fx%x = xup(B%x)**2
    fy%y = yup(B%y)**2
    fz%z = zup(B%z)**2
    ! Magnetic pressure
    pb = 0.5*(fx%x+fy%y+fz%z)
    ! Fast mode speed + velocity
    pa = self%nu(3)*d*max(du,0.0)**2
    cs = sqrt((self%gamma*pg+2.0*pb + pa)/d)
    self%u_max = maxval(cs + sqrt(pp) + abs(du))
    du = max(du,0.0)
    ! Isotropic pressure part
    pp = d*pp + pg + pa + pb
    ! Total diagondal stress
    fx%x = pp - fx%x - d*du*ddxup(U%x)
    fy%y = pp - fy%y - d*du*ddyup(U%y)
    fz%z = pp - fz%z - d*du*ddzup(U%z)
  else
    ! Diagonal Reynolds stress
    fx%x = xup(U%x)**2
    fy%y = yup(U%y)**2
    fz%z = zup(U%z)**2
    pp = (fx%x+fy%y+fz%z)
    ! Sound speed + velocity
    cs = sqrt(self%gamma*pg/d)
    self%u_max = maxval(cs + sqrt(pp) + du)
    ! Isotropic pressure part
    pa = self%nu(3)*d*du**2                               ! dimension: pressure
    pp = d*pp + pg + pa
    ! Total diagondal stress, with density factor
    fx%x = pp - d*du*ddxup(U%x)
    fy%y = pp - d*du*ddyup(U%y)
    fz%z = pp - d*du*ddzup(U%z)
  end if
  du = (self%nu(1)*cs + self%nu(2)*du)*dx               ! dimension: U*ds
  !-----------------------------------------------------------------------------
  ! Conservation of mass
  !-----------------------------------------------------------------------------
  fd = -(self%nu(5)*dd*down(du))*ddown(lnd)          ! diffusive mass flux
  fm = p+fd                                             ! mass flux
  dddt = -div(fm)                                    ! density derivative
  !-----------------------------------------------------------------------------
  ! Conservation of entropy.  The formulation that works best so far is to make
  ! the diffusive entropy flux directly proportional to the diffusive mass flux,
  ! plus a contribution proportional to the gradient of specific entropy ss.
  !-----------------------------------------------------------------------------
  if (self%gamma==1.0) then
    dsdt = 0.0
  else
    ss = s/d
    fm = fm*down(ss) - (self%nu(3)*dd*down(du))*ddown(ss)
    dsdt = -div(fm)
  end if
  !---------------------------------------------------------------------------
  ! Off-diagonal part of Reynolds stress, for now without the density factor.
  ! Note that, because of the symmmetry, we define only three components
  !---------------------------------------------------------------------------
  fx%y = ydn(U%x)*xdn(U%y)
  fx%z = zdn(U%x)*xdn(U%z)
  fy%z = zdn(U%y)*ydn(U%z)
  !---------------------------------------------------------------------------
  ! Add viscous part -- edge-centered
  !---------------------------------------------------------------------------
  fx%y = fx%y - self%nu(4)*xdn1(ydn1(du))*(ddxdn(U%y)+ddydn(U%x))
  fx%z = fx%z - self%nu(4)*ydn1(zdn1(du))*(ddxdn(U%z)+ddzdn(U%x))
  fy%z = fy%z - self%nu(4)*zdn1(xdn1(du))*(ddydn(U%z)+ddzdn(U%y))
  !---------------------------------------------------------------------------
  ! Multiply with edge-centered density
  !---------------------------------------------------------------------------
  fx%y = exp(ydn(ld%x))*fx%y
  fx%z = exp(zdn(ld%x))*fx%z
  fy%z = exp(zdn(ld%y))*fy%z
  !---------------------------------------------------------------------------
  ! Off-diagonal Maxwell stress
  !---------------------------------------------------------------------------
  if (self%mhd) then
    fx%y = fx%y - ydn(B%x)*xdn(B%y)
    fx%z = fx%z - zdn(B%x)*xdn(B%z)
    fy%z = fy%z - zdn(B%y)*ydn(B%z)
  end if
  !-----------------------------------------------------------------------------
  ! Consistent addition of mass flux related momentum flux, symmetrized. The
  ! diagonal part is analogous, a centered version of \nu \rho U_i \grad\ln\rho,
  ! while the off-diagonal part needs to be properly symmetrized.
  !-----------------------------------------------------------------------------
  pa = d*du*self%nu(5)
  fx%x = fx%x - du*ddxup(ld%x)*xup(U%x)
  fy%y = fy%y - du*ddyup(ld%y)*yup(U%y)
  fz%z = fz%z - du*ddzup(ld%z)*zup(U%z)
  fx%y = fx%y + (xdn1(fd%y)*ydn(U%x)+ydn1(fd%x)*xdn(U%y))
  fx%z = fx%z + (xdn1(fd%z)*zdn(U%x)+zdn1(fd%x)*xdn(U%z))
  fy%z = fy%z + (zdn1(fd%y)*ydn(U%z)+ydn1(fd%z)*zdn(U%y))
  !-----------------------------------------------------------------------------
  ! Momentum equations
  !-----------------------------------------------------------------------------
  dpdt%x = - ddxdn1(fx%x) - ddyup1(fx%y) - ddzup1(fx%z)
  dpdt%y = - ddxup1(fx%y) - ddydn1(fy%y) - ddzup1(fy%z)
  dpdt%z = - ddxup1(fx%z) - ddyup1(fy%z) - ddzdn1(fz%z)
  if (associated(self%force%selected)) &
    dpdt = dpdt + dd*self%force%selected (self%time, self%mesh)
  if (self%mhd) then
    !---------------------------------------------------------------------------
    ! Advective electric field on cell edges
    !---------------------------------------------------------------------------
    E%x = ydn(U%z)*zdn(B%y) - zdn(U%y)*ydn(B%z)
    E%y = zdn(U%x)*xdn(B%z) - xdn(U%z)*zdn(B%x)
    E%z = xdn(U%y)*ydn(B%x) - ydn(U%x)*xdn(B%y)
    !---------------------------------------------------------------------------
    ! Induction equation, with resistive term
    !---------------------------------------------------------------------------
    E = E + self%nu(6)*edge(du)*curl_dn(B)
    dBdt = -curl_up(E)
  end if
  !-----------------------------------------------------------------------------
  call deallocate_vectors (ld, dd, U, E, fd, fx, fy, fz, fm)
  call deallocate_scalars (lnd, pg, pa, du, ss, cs, pp, pb)
  call trace_end
  end associate
END SUBROUTINE update

!> =============================================================================
FUNCTION gas_pressure (self) RESULT (pg)
  class(mhd_t):: self
  type(scalar_t):: pg
  !-----------------------------------------------------------------------------
  if (self%gamma==1.0) then
    pg = self%d
  else
    pg = self%d**self%gamma*exp(self%s/self%d*(self%gamma-1))
  end if
END FUNCTION gas_pressure

!> =============================================================================
SUBROUTINE log_density (self, v)
  class (mhd_t):: self
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = log(self%d%f)
END SUBROUTINE

!> =============================================================================
SUBROUTINE log_pressure (self, v)
  class (mhd_t):: self
  type (scalar_t):: s
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  call s%init (self%d%n, self%d%ds)
  !s = self%gas_pressure()
  v = log(s%f)
  call s%deallocate()
END SUBROUTINE

!> =============================================================================
SUBROUTINE velocity_magnitude (self, v)
  class (mhd_t):: self
  real, dimension(:,:,:), pointer:: v
  type (scalar_t):: s
  type(vector_t):: u
  !-----------------------------------------------------------------------------
  call s%init(self%d%n, self%d%ds)
  call u%init(self%d%n, self%d%ds)
  s = log(self%d)
  U = down(s)
  U = self%p/U
  s = norm(U)
  v = s%f
  call s%deallocate()
  call u%deallocate()
END SUBROUTINE

!> =============================================================================
SUBROUTINE magnetic_field_magnitude (self, v)
  class (mhd_t):: self
  real, dimension(:,:,:), pointer:: v
  type (scalar_t):: s
  type(vector_t):: u
  !-----------------------------------------------------------------------------
  call s%init(self%d%n, self%d%ds)
  s = norm(self%B)
  v = s%f
  call s%deallocate()
END SUBROUTINE

!> =============================================================================
SUBROUTINE void (self, v)
  class (mhd_t):: self
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = 0.0
END SUBROUTINE

END MODULE mhd_mod
