
!*******************************************************************************
SUBROUTINE pde (d, s, p, B, dddt, dsdt, dpdt, dBdt)     ! conserved variables
  USE scalar_mod
  USE vector_mod
  USE stagger_mod
  implicit none
  type(scalar_t):: d, s, dddt, dsdt                     ! density, entropy
  type(vector_t):: p, B, dpdt, dBdt                     ! momentum, magnetic field
  !.............................................................................
  type(vector_t):: U, E, J                              ! velocity, electric field
  type(vector_t):: Us, Bs, fx, fy, fz, fd, ld, dd       ! extra vector temporaries
  type(scalar_t):: du, lnd, ss, cs, pg, pm, pa          ! extra scalar temporaries
  !.............................................................................
  real:: nu(6)=1.0, ds, gamma=1.0
  logical:: mhd=.true.
  !-----------------------------------------------------------------------------
  call init_vectors (d%n, d%ds, ld, dd, U, E, fd, fx, fy, fz)
  call init_scalars (d%n, d%ds, lnd, pg, pa, du, ss, cs, pm)
  !-----------------------------------------------------------------------------
  ! Compute velocity U=momentum/density, leaving down shifted density and log
  ! density for reuse below
  !-----------------------------------------------------------------------------
  lnd = log(d)
  ld = down(lnd)
  dd = exp(lnd)
  U = p/dd
  !-----------------------------------------------------------------------------
  ! Advective electric field on cell edges
  !-----------------------------------------------------------------------------
  E%x = ydn(U%z)*zdn(B%y) - zdn(U%y)*ydn(B%z)
  E%y = zdn(U%x)*xdn(B%z) - xdn(U%z)*zdn(B%x)
  E%z = xdn(U%y)*ydn(B%x) - ydn(U%x)*xdn(B%y)
  !-----------------------------------------------------------------------------
  ! Gas pressure
  !-----------------------------------------------------------------------------
  if (gamma==1.0) then                                  ! isothermal
    pg = d                                              ! temperature = 1
    cs = 1.0                                            ! sound speed
  else
    pg = d**gamma*exp(s/d*(gamma-1))                    ! ideal gas pressure
    cs = sqrt(gamma*pg/d)                               ! sound speed
  end if
  !-----------------------------------------------------------------------------
  ! Diffusion coefficient; the result 'nu' is cell centered, and may be staggered
  ! down without increasing the number of ghost zones needed
  !-----------------------------------------------------------------------------
  ds = minval(d%ds)                                     ! smallest {dx,dy,dz}
  du = -div(U)
  du = ds*max(du,1e-20)                                 ! dimension: U
  pa = nu(3)*d*du**2                                    ! dimension: pressure
  du = (nu(1)*cs + nu(2)*du)*ds                         ! dimension: U*ds
  !-----------------------------------------------------------------------------
  ! Conservation of mass
  !-----------------------------------------------------------------------------
  fd = nu(5)*dd*ddown(lnd)*down(du)                     ! diffusive mass flux
  fx = p-fd
  dddt = -div (fx)                                 ! density derivative
  !-----------------------------------------------------------------------------
  ! Conservation of entropy.  The formulation that works best so far is to make
  ! the diffusive entropy flux directly proportional to the diffusive mass flux.
  ! For consistency, the same terms are added in the equations of motion below
  !-----------------------------------------------------------------------------
  if (gamma==1.0) then
    dsdt = 0.0
  else
    ss = s/d
    fx = fx*down(ss) - nu(6)*dd*ddown(ss)*down(du)
    dsdt = -div (fx)
  end if
  !-----------------------------------------------------------------------------
  ! Conservation of momentum
  !-----------------------------------------------------------------------------
  if (mhd) then
    pm = 0.5*norm2(B)
    fx = d*norm2(U) + pg - pm
  else
    fx = d*norm2(U) + pg
  end if
  !---------------------------------------------------------------------------
  ! Add diagonal part of viscous stress
  !---------------------------------------------------------------------------
  fx = fx + pa - d*du*dup(U)
  fy = fy + pa - d*du*dup(U)
  fz = fz + pa - d*du*dup(U)
  !---------------------------------------------------------------------------
  ! Off-diagonal part of Reynolds stress, without the density factor
  !---------------------------------------------------------------------------
  fx%y = ydn(U%x)*xdn(U%y)
  fx%z = zdn(U%x)*xdn(U%z)
  fy%z = zdn(U%y)*ydn(U%z)
  !---------------------------------------------------------------------------
  ! Add viscous part. nu(4) is a multiplier ! (relative strength) in front of
  ! the off-diagonal stresses.
  !---------------------------------------------------------------------------
  if (nu(4)>0.0) then
    fx%y = fx%y - nu(4)*xdn1(ydn1(du))*(ddxdn(U%y)+ddydn(U%x))
    fx%z = fx%z - nu(4)*ydn1(zdn1(du))*(ddxdn(U%z)+ddzdn(U%x))
    fy%z = fy%z - nu(4)*zdn1(xdn1(du))*(ddydn(U%z)+ddzdn(U%y))
  end if
  !---------------------------------------------------------------------------
  ! Multiply with edge-centered density
  !---------------------------------------------------------------------------
  fx%y = exp(ydn(ld%x))*fx%y
  fx%z = exp(zdn(ld%x))*fx%z
  fy%z = exp(zdn(ld%y))*fy%z
  !---------------------------------------------------------------------------
  ! Off-diagonal part of Maxwell stress
  !---------------------------------------------------------------------------
  if (mhd) then
    fx%y = fx%y - ydn(B%x)*xdn(B%y)
    fx%z = fx%z - zdn(B%x)*xdn(B%z)
    fy%z = fy%z - zdn(B%y)*ydn(B%z)
  end if
  !-----------------------------------------------------------------------------
  ! Consistent addition of mass flux related momentum flux, symmetrized
  !-----------------------------------------------------------------------------
  pa = d*du*nu(5)
  fx = fx - pa*ddxup(ld%x)*xup(U%x)
  fy = fy - pa*ddyup(ld%y)*yup(U%y)
  fz = fz - pa*ddzup(ld%z)*zup(U%z)
  fx%y = fx%y - (xdn1(fd%y)*ydn(U%x)+ydn1(fd%x)*xdn(U%y))
  fx%z = fx%z - (xdn1(fd%z)*zdn(U%x)+zdn1(fd%x)*xdn(U%z))
  fy%z = fy%z - (zdn1(fd%y)*ydn(U%z)+ydn1(fd%z)*zdn(U%y))
  !-----------------------------------------------------------------------------
  ! Momentum equations
  !-----------------------------------------------------------------------------
  dpdt%x = - ddxdn1(fx%x) - ddyup1(fx%y) - ddzup1(fx%z)
  dpdt%y = - ddydn1(fy%y) - ddxup1(fx%y) - ddzup1(fy%z)
  dpdt%z = - ddzdn1(fz%z) - ddxup1(fx%z) - ddyup1(fy%z)
  !-----------------------------------------------------------------------------
  ! Induction equation
  !-----------------------------------------------------------------------------
  E = E + nu(6)*edge(du)*curl(B)
  dBdt = -curl(E)
  !-----------------------------------------------------------------------------
  call deallocate_vectors (ld, dd, U, E, fd, fx, fy, fz)
  call deallocate_scalars (lnd, pg, pa, du, ss, cs, pm)
END SUBROUTINE pde
