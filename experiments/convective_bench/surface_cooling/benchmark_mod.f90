!===============================================================================
!> Generic module implementing initial conditions, boundare conditions, cooling,
!> and initial perturbations.  To include in your code and call procedures:
!>
!>   USE benchmark_mod
!>   ...
!>   call benchmark%initial (z, d, e, u, v, w)
!>   ...
!>   call benchmark%cooling (x, y, z, d)
!>   ...
!>   call benchmark%top (z(1:lb), d(:,:,1:lb), e(:,:,1:lb), &
!>                   u(:,:,1:lb), v(:,:,1:lb), w(:,:,1:lb))
!>   ...
!>   call benchmark%bottom (z(ub:), d(:,:,ub:), e(:,:,ub:), &
!>                      u(:,:,ub:), v(:,:,ub:), w(:,:,ub:))
!>   ...
!> Note that this assumes depth index increasing downwards.  We use the notation
!> "top" for the lower index boundary region, with lb the first internal index,
!> and "bottom" for the upper index boundaries, with ub the last internal index.
!===============================================================================
MODULE benchmark_mod
  implicit none
  real, private:: p_bb=1., d_bb=1., z_bb=1.
  type benchmark_t
    real::  n_poly = 10       ! polytropic index
    real::      z0 = 0.2      ! break to constant T
    real::      z1 = 1.0      ! extended range of polytrope
    real::     uz0 = 0.0      ! bulk vertical velocity (for tests)
    real::  q_rate = 1.       ! cooling rate
    real:: perturb = 1e-3     ! perturbation amplitude.
    real:: gamma, g           ! derived quantities
  contains
    procedure:: initial
    procedure:: top
    procedure:: bottom
    procedure:: cooling
    procedure:: perturbations
  end type
  type(benchmark_t), public:: benchmark
CONTAINS

!===============================================================================
!> Initial condition -- assumes there is no patch entirely above z0
!===============================================================================
SUBROUTINE initial (self, z, d, e, u, v, w)
  class(benchmark_t):: self
  real(8), dimension(:):: z
  real, dimension(:,:,:):: d, e, u, v, w
  real:: t, p, p1, t1, h, z1, s, zp
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  self%gamma = 1.0 + 1.0/self%n_poly
  self%g = (self%n_poly + 1.)/(1. + self%z1)
  n = shape(e)
  do k=n(3),1,-1
    !---------------------------------------------------------------------------
    ! Below z0, use polytrope
    !---------------------------------------------------------------------------
    if (z(k) > self%z0) then
      zp = (z(k)+self%z1)/(1+self%z1)
      t = zp
      p = zp**(self%n_poly+1)
      do j=1,n(2)
      do i=1,n(1)
        d(i,j,k) = p/t
        e(i,j,k) = p/(self%gamma-1.)
      end do
      end do
      p1 = p
      t1 = t
      z1 = z(k)
      h = t/self%g
    !---------------------------------------------------------------------------
    ! Above z0, use constant temperature
    !---------------------------------------------------------------------------
    else
      p = p1*exp((z(k)-z1)/h)
      do j=1,n(2)
      do i=1,n(1)
        d(i,j,k) = p/t1
        e(i,j,k) = p/(self%gamma-1.)
      end do
      end do
    end if
  end do
  u = 0.0
  v = 0.0
  w = self%uz0
END SUBROUTINE initial

!===============================================================================
!> Top BC: vanishing vertical derivatives of velocity components, constant
!> local temperature and scale heights in the guard zones
!===============================================================================
SUBROUTINE top (self, z, d, e, u, v, w)
  class(benchmark_t):: self
  real(8), dimension(:):: z
  real, dimension(:,:,:):: d, e, u, v, w
  real:: s, t1, h, p, p1, z1
  integer:: i, j, k, n(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  n = shape(e)
  do j=1,n(2)
  do i=1,n(1)
    p1 = e(i,j,n(3))*(self%gamma-1.)
    z1 = z(n(3))
    t1 = p1/d(i,j,n(3))
    h = t1/self%g
    do k=1,n(3)-1
      p = p1*exp((z(k)-z1)/h)
      e(i,j,k) = p/(self%gamma-1.)
      d(i,j,k) = p/t1
      u(i,j,k) = u(i,j,n(3))
      v(i,j,k) = v(i,j,n(3))
      w(i,j,k) = w(i,j,n(3))
    end do
  end do
  end do
END SUBROUTINE top

!===============================================================================
!> Bottom BC:  Constant gas pressure p_bb=1, constant entropy s_bb=0 for
!> ascending gas, while descending gas keeps its entropy
!===============================================================================
SUBROUTINE bottom (self, z, d, e, u, v, w)
  class(benchmark_t):: self
  real(8), dimension(:):: z
  real, dimension(:,:,:):: d, e, u, v, w
  real:: t, h, p, p1, d1, t1, s1
  integer:: i, j, k, n(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  n = shape(e)
  do j=1,n(2)
  do i=1,n(1)
  do k=2,n(3)
    !---------------------------------------------------------------------------
    ! Compute the local pressure scale height and extrapolate pressure with that
    !---------------------------------------------------------------------------
    p1 = e(i,j,1)*(self%gamma-1.)
    d1 = d(i,j,1)
    t1 = p1/d1
    h = t1/self%g
    p = p_bb*exp((z(k)-z_bb)/h)
    e(i,j,k) = p/(self%gamma-1.)
    !---------------------------------------------------------------------------
    ! Set the entropy per unit mass of ascending gas to the initial s_bb=0,
    ! while leaving entropy of descending gas untouched
    !---------------------------------------------------------------------------
    if (w(i,j,1) < 0.0) then
      d(i,j,k) = d_bb*(p/p_bb)**(1.0/self%gamma)
    else
      d(i,j,k) = d1*(p/p1)**(1.0/self%gamma)
    end if
    !---------------------------------------------------------------------------
    ! Vanishing vertical derivative of velocity components
    !---------------------------------------------------------------------------
    u(i,j,k) = u(i,j,1)
    v(i,j,k) = v(i,j,1)
    w(i,j,k) = w(i,j,1)
  end do
  end do
  end do
END SUBROUTINE bottom

!===============================================================================
!> Surface cooling -- Newtonian cooling towards T0 in surface layers
!===============================================================================
SUBROUTINE cooling (self, z, d, e, q)
  class(benchmark_t):: self
  real(8), dimension(:):: z
  real, dimension(:,:,:):: d, e, q
  integer:: n(3), i, j, k
  real:: dz, rate, T, T0
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  n = shape(d)
  dz = 0.1*self%z0
  T0 = 0.5*(self%z0+self%z1)/(1+self%z1)
  do k=1,n(3)
    rate = self%q_rate/(self%gamma-1.0)*(0.5 + 0.5*tanh((self%z0-z(k))/dz))
    do j=1,n(2)
    do i=1,n(1)
      T = (self%gamma-1.0)*e(i,j,k)/d(i,j,k)
      q(i,j,k) = (T0-T)*rate
    end do
    end do
  end do
END SUBROUTINE cooling

!===============================================================================
!> Horizontal perturbations for predictive triggering of convection; using all
!> wavenumbers within a radius of 4, with phase and amplitudes chosen to avoid
!> symmetries
!===============================================================================
SUBROUTINE perturbations (self, x, y, z, d)
  class(benchmark_t):: self
  real(8), pointer:: x(:), y(:), z(:)
  real, dimension(:,:,:):: d
  real, parameter:: pi = 4.*atan(1.)
  integer:: i, j, k, k1, k2, n(3)
  real:: a, kx, ky, d1, f
  !-----------------------------------------------------------------------------
  n = shape(d)
  do k=1,n(3)
    d1 = d(1,1,k)
    do j=1,n(2)
    do i=1,n(1)
      f = 1.0
      do k1=1,4
      do k2=1,4
        if (sqrt(real(k1)**2+real(k2)**2) <= 4.) then
          kx = 2.*pi/4.*k1*x(i) - 11.3*k1**0.5
          ky = 2.*pi/4.*k2*y(j) - 50.1*k2**0.2
          a = self%perturb/(1+(k1-1)+3*(k2-1))
          f = f + a*sin(kx)*sin(ky)*exp(-(z(k)-self%z0)**2/0.01)
        end if
      end do
      end do
      d(i,j,k) = d1*f
    end do
    end do
  end do
END SUBROUTINE perturbations

END MODULE benchmark_mod
