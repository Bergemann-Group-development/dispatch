!===============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE vector_mod
  USE mesh_mod
  USE atmosphere_mod
  USE planet_mod
  USE lagrange_mod
  USE units_mod
  USE scaling_mod
  USE index_mod
  implicit none
  private
  real(8), parameter:: pi2=8d0*atan(1d0)
  type, public:: initial_t
    real:: b0(3)
    real:: gamma=1.4
    real:: surface_temperature, surface_pressure
    real(8):: size(3)
    logical:: mhd
    procedure(planet_atmosphere), pointer:: condition => null()
  contains
    procedure init
  end type
  public planet_atmosphere
  type(lagrange_t):: interp
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self, solver)
   class(initial_t):: self
   character(len=*):: solver
   !............................................................................
   logical, save:: first_time=.true.
   real, save:: b0(3) = 0.0
   namelist /IC_params/ b0
   !----------------------------------------------------------------------------
   call trace_begin('initial_t%init')
   call scaling%init
   if (first_time) then
     first_time = .false.
     rewind (io%input); read (io%input, initial_params)
     if (io%master) write (*, initial_params)
   end if
   self%b0 = b0
   self%condition => planet_atmosphere
   call trace_end
END SUBROUTINE init

!===============================================================================
!> Interpolate in a hydrostatic planet atmosphere, computed in atmosphere_mod
!===============================================================================
SUBROUTINE planet_atmosphere (self, m, mem, idx)
  class(initial_t)             :: self
  type(index_t):: idx
  real, dimension(:,:,:,:)     :: mem
  class(mesh_t), pointer, dimension(:) :: m
  !.............................................................................
  integer :: ix, iy, iz, nr
  real(8) :: x, y, z, r, d, p
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t%planet_atmosphere')
  nr = atmosphere%nr
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + m(3)%r(iz) - planet%position(3)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + m(2)%r(iy) - planet%position(2)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + m(1)%r(ix) - planet%position(1)
        r = sqrt(x**2+y**2+z**2)
        d = interpolate (r, nr ,atmosphere%rr, atmosphere%rho)
        p = interpolate (r, nr, atmosphere%rr, atmosphere%pg)
        mem(ix,iy,iz,1) = d
        mem(ix,iy,iz,2) = d*(log(p)-self%gamma*log(d))/(self%gamma-1.0)
        mem(ix,iy,iz,3) = 0.0
        mem(ix,iy,iz,4) = 0.0
        mem(ix,iy,iz,5) = 0.0
        mem(ix,iy,iz,6) = self%b0(1)
        mem(ix,iy,iz,7) = self%b0(2)
        mem(ix,iy,iz,8) = self%b0(3)
      end do
    end do
  end do
  call trace_end
END SUBROUTINE planet_atmosphere

!===============================================================================
!> Use default Lagrange interpolation
!===============================================================================
FUNCTION interpolate (x0, n, x, y) RESULT(f)
  implicit none
  integer:: n
  real(kind=8):: x0, x(n), y(n), f
  !.............................................................................
  interp%verbose = max(io%verbose-4,0)
  if (x0 < x(1)) then
    f = y(1)
  else if (x0 > x(n)) then
    f = y(n)
  else
    f = interp%sequence (x0, x, y, n)
  end if
END FUNCTION

!===============================================================================
!> Solve the equations for hydrostatic balance, initially neglecting the mass
!> of the pseudo-atmosphere -- later on we will include it.
!===============================================================================
SUBROUTINE hydrostatic_balance (nr, rr, rho, pg, temp, mass)
  implicit none
  integer:: nr, i, iter
  real(kind=8), dimension(nr):: lnr, rho, pg, temp, rr
  real(kind=8):: rho1, dr, dlnpg, mass
  integer, parameter:: niter=10
  !-----------------------------------------------------------------------------
  ! Start with the external pressure and integrate inwards
  !-----------------------------------------------------------------------------
  rho(nr) = pg(nr)/temp(nr)
  mass = 0.0
  do i=nr-1,1,-1
    pg(i) = pg(i+1)
    do iter=1,niter
      dr = rr(i+1)-rr(i)
      rho(i) = pg(i)/temp(i)
      dlnpg = dr*code%grav*mass* &
        0.5*(rho(i+1)/pg(i+1)/rr(i+1)**2 + rho(i)/pg(i)/rr(i)**2)
      pg(i) = pg(i+1)*exp(dlnpg)
    end do
    mass = mass + 4.*cgs%pi*(rr(i+1)**3-rr(i)**3)/3.*(rho(i+1)+rho(i))*0.5
!    print'("hydrostatic",i5,1p,8e11.3)',i,rr(i),pg(i),temp(i)*code%temp,rho(i),dr,dlnpg,planet%mass*code%grav/rr(i)**2
  end do
END SUBROUTINE hydrostatic_balance

END MODULE
