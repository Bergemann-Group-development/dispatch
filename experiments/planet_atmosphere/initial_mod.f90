!> =============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!> =============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE vector_mod
  USE mesh_mod
  USE lagrange_mod, only: lagrange_sequence
  USE units_mod
  USE scaling_mod
  implicit none
  private
  real(8), parameter:: pi2=8d0*atan(1d0)
  type, public:: initial_t
    real:: b0(3)
    real:: gamma=1.4
    real:: surface_temperature, surface_pressure
    real(8):: size(3)
    logical:: mhd
    procedure(planet_atmosphere), pointer:: selected => null()
  contains
    procedure init
  end type
  logical, save:: first_time=.true.
  public planet_atmosphere
CONTAINS

!> =============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!> =============================================================================
SUBROUTINE init (self)
   class(initial_t):: self
   logical, save:: mhd=.true.
   real, save:: b0(3) = 0.0
   real, save:: size_re=1.0, surface_temperature_k=3000., surface_pressure_cgs=1e10
   namelist /initial_params/ size_re, surface_temperature_k, surface_pressure_cgs, b0
   !----------------------------------------------------------------------------
   call trace_begin('initial_t%init')
   if (first_time) then
    !$omp critical (experiment_cr)
     if (first_time) then
       first_time = .false.
       rewind (io%input)
       read (io%input, initial_params)
       mhd = any(self%b0 /= 0.0)
       if (io%master) write (*, initial_params)
     end if
    !$omp end critical (experiment_cr)
   end if
   self%size = size_re*cgs%r_earth/scaling%l
   self%b0 = b0
   self%surface_temperature = surface_temperature_k/scaling%temp
   self%surface_pressure = surface_pressure_cgs/scaling%p
   self%selected => planet_atmosphere
   call trace_end
END SUBROUTINE init

!> =============================================================================
!> Integrate a hydrostatic and constant entropy planet atmosphere, with initially
!> just a central gravitational field.
!> =============================================================================
FUNCTION planet_atmosphere (self, m) RESULT (ff)
  class(initial_t)      :: self
  type(mesh_t)          :: m(3)
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,8) :: ff
  !.............................................................................
  integer               :: ix, iy, iz
  real(8)               :: x, y, z, r, d, p
  integer, parameter    :: nrad=1000
  real(8), dimension(nrad) :: dd, pp, rr
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t%planet_atmosphere')
  if (first_time) call self%init
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + m(3)%r(iz)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + m(2)%r(iy)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + m(1)%r(ix)
        r = sqrt(x**2+y**2+z**2)
        d = interpolate(r,nrad,dd,rr)
        p = interpolate(p,nrad,pp,rr)
        ff(ix,iy,iz,1) = d
        ff(ix,iy,iz,3) = (log(p)-self%gamma*log(d))/(self%gamma-1.0)
        ff(ix,iy,iz,3) = 0.0
        ff(ix,iy,iz,4) = 0.0
        ff(ix,iy,iz,5) = 0.0
        ff(ix,iy,iz,6) = self%b0(1)
        ff(ix,iy,iz,7) = self%b0(2)
        ff(ix,iy,iz,8) = self%b0(3)
      end do
    end do
  end do
  do ix=m(1)%lb,m(1)%ub
    print'(a,i4,f22.16)', 'initial:', ix, ff(ix,1,1,3)
  end do
  call trace_end
END FUNCTION planet_atmosphere

!> =============================================================================
!> Use default Lagrange interpolation
!> =============================================================================
FUNCTION interpolate (x0, n, x, y) RESULT(f)
  implicit none
  integer:: n
  real(kind=8):: x0, x(n), y(n), f
  !.............................................................................
  f = lagrange_sequence (x0, x, y, n)
END FUNCTION

!> =============================================================================
!> Solve the equations for hydrostatic balance, initially neglecting the mass
!> of the pseudo-atmosphere -- later on we will include it.
!> =============================================================================
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
