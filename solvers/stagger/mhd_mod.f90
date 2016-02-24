!> =============================================================================
!> =============================================================================
MODULE mhd_mod
  USE io_mod
  USE trace_mod
  USE scalar_mod
  USE vector_mod
  USE stagger_mod
  USE patch_mod
  implicit none
  private
  type, extends(patch_t), public:: mhd_t
    type(scalar_t):: d, dddt, s, dsdt
    type(vector_t):: p, dpdt, B, dBdt
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
  end if
  !$omp end critical (mhd_cr)
  call trace_end
END SUBROUTINE init

!> =============================================================================
SUBROUTINE update (self)
  class(mhd_t):: self
  !.............................................................................
  type(vector_t):: U, E, J                              ! velocity, electric field
  type(vector_t):: fx, fy, fz, fd, ld, dd               ! extra vector temporaries
  type(scalar_t):: du, lnd, ss, cs, pg, pa, pp, pb      ! extra scalar temporaries
  real, dimension(:,:,:),   pointer:: s
  real, dimension(:,:,:,:), pointer:: v
  !.............................................................................
  integer:: n(3)
  real:: ds
  !-----------------------------------------------------------------------------
  call trace_begin('mhd_t%update')
  if (first_time) call self%init
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
  call init_vectors (self%mesh%gn, self%ds, ld, dd, U, E, fd, fx, fy, fz)
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
  cs = sqrt(self%gamma*pg/d)                            ! sound speed
  !-----------------------------------------------------------------------------
  ! Diffusion coefficient; the result 'nu' is cell centered, and may be staggered
  ! down without increasing the number of ghost zones needed
  !-----------------------------------------------------------------------------
  ds = minval(d%ds)                                     ! smallest {dx,dy,dz}
  du = -div(U)                                          ! convergence rate
  du = ds*max(du,1e-20)                                 ! dimension: U
  pa = self%nu(3)*d*du**2                               ! dimension: pressure
  du = (self%nu(1)*cs + self%nu(2)*du)*ds               ! dimension: U*ds
  !-----------------------------------------------------------------------------
  ! Conservation of mass
  !-----------------------------------------------------------------------------
  fd = -(self%nu(5)*dd*down(du))*ddown(lnd)             ! diffusive mass flux
  fx = p+fd                                             ! mass flux
  dddt = -div(fx)                                       ! density derivative
  block
    integer:: l(5)
    l=0
    call self%check_error('mhd_t%update dddt',l)
    if (any(l>0)) then
      print*,'dddt'
      print 1,dddt%f(l(1),l(2),:)
      print*,'fd'
      print 1,fd%z%f(l(1),l(2),:)
      print*,'fx'
      print 1,fx%z%f(l(1),l(2),:)
      print*,'dd'
      print 1,dd%z%f(l(1),l(2),:)
      print*,'p'
      print 1,p%z%f(l(1),l(2),:)
      print*,'du'
      print 1,du%f(l(1),l(2),:)
      print*,'lnd'
      print 1,lnd%f(l(1),l(2),:)
      print*,'U'
      print 1,u%x%f(l(1),l(2),:)
      print 1,u%y%f(l(1),l(2),:)
      print 1,u%z%f(l(1),l(2),:)
    1 format(1p,12e10.2)
      stop
    end if
  end block
  !-----------------------------------------------------------------------------
  ! Conservation of entropy.  The formulation that works best so far is to make
  ! the diffusive entropy flux directly proportional to the diffusive mass flux,
  ! plus a contribution proportional to the gradient of specific entropy ss.
  !-----------------------------------------------------------------------------
  if (self%gamma==1.0) then
    dsdt = 0.0
  else
    ss = s/d
    fx = fx*down(ss) - (self%nu(6)*dd*down(du))*ddown(ss)
    dsdt = -div(fx)
  end if
  call self%check_error('mhd_t%update 1')
  !-----------------------------------------------------------------------------
  ! Conservation of momentum -- first the diagonal part of the stress tensor
  !-----------------------------------------------------------------------------
  if (self%mhd) then
    pb = 0.5*norm2(B)
    pp = norm2(U)
    self%u_max = maxval(sqrt((self%gamma*pg+2.*pb)/d) + sqrt(pp))
    pp = d*pp + pg + pa - pb
  else
    pp = d*norm2(U) + pg + pa
  end if
  fx%x = pp - d*du*ddxup(U%x)
  fy%y = pp - d*du*ddyup(U%y)
  fz%z = pp - d*du*ddzup(U%z)
  call self%check_error('mhd_t%update 2')
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
  call self%check_error('mhd_t%update 3')
  !---------------------------------------------------------------------------
  ! Off-diagonal part of Maxwell stress
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
  fx%x = fx%x - pa*ddxup(ld%x)*xup(U%x)
  fy%y = fy%y - pa*ddyup(ld%y)*yup(U%y)
  fz%z = fz%z - pa*ddzup(ld%z)*zup(U%z)
  fx%y = fx%y + (xdn1(fd%y)*ydn(U%x)+ydn1(fd%x)*xdn(U%y))
  fx%z = fx%z + (xdn1(fd%z)*zdn(U%x)+zdn1(fd%x)*xdn(U%z))
  fy%z = fy%z + (zdn1(fd%y)*ydn(U%z)+ydn1(fd%z)*zdn(U%y))
  !-----------------------------------------------------------------------------
  ! Momentum equations
  !-----------------------------------------------------------------------------
  call self%check_error('mhd_t%update 7')
  dpdt%x = - ddxdn1(fx%x) - ddyup1(fx%y) - ddzup1(fx%z)
  call self%check_error('mhd_t%update 7x')
  dpdt%y = - ddxup1(fx%y) - ddydn1(fy%y) - ddzup1(fy%z)
  call self%check_error('mhd_t%update 7y')
  print*,dpdt%z%id,dpdt%z%temporary,dpdt%z%allocated,dpdt%z%n,dpdt%z%l,dpdt%z%u,dpdt%z%ds
  dpdt%x = - ddxup1(fx%z)
  call self%check_error('mhd_t%update 7a')
  dpdt%x =- ddyup1(fy%z)
  call self%check_error('mhd_t%update 7b')
  dpdt%x = - ddzdn1(fz%z)
  call self%check_error('mhd_t%update 7c')
  dpdt%z = - ddxup1(fx%z) - ddyup1(fy%z) - ddzdn1(fz%z)
  call self%check_error('mhd_t%update 7z')
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
    E = E + self%nu(6)*edge(du)*curl(B)
    dBdt = -curl(E)
    call self%check_error('mhd_t%update 8')
  end if
  !-----------------------------------------------------------------------------
  call deallocate_vectors (ld, dd, U, E, fd, fx, fy, fz)
  call deallocate_scalars (lnd, pg, pa, du, ss, cs, pp, pb)
  call self%check_error('mhd_t%update 9')
  !call cleanup_temp
  end associate
  call trace_end
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
