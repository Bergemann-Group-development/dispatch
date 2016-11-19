!> =============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!> =============================================================================
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
  implicit none
  private
  real(8), parameter:: pi2=8d0*atan(1d0)
  type, public:: initial_t
    real:: b0(3)
    real:: gamma
    real:: v_head, v_accrete
    real:: surface_temperature, surface_pressure
    real(8):: size(3)
    logical:: mhd
    procedure(planet_atmosphere), pointer:: selected => null()
  contains
    procedure init
  end type
  public planet_atmosphere
  type(lagrange_t):: interp
CONTAINS

!> =============================================================================
!> Setup an initial state, with head wind speed and accretion speed
!> =============================================================================
SUBROUTINE init (self)
   class(initial_t):: self
   logical, save:: first_time=.true.
   real, save:: b0(3) = 0.0, gamma=1.4, v_head_cgs=0.0, v_accrete_cgs=0.0
   namelist /initial_params/ b0, gamma, v_head_cgs, v_accrete_cgs
   !----------------------------------------------------------------------------
   call trace_begin('initial_t%init')
   call scaling%init
   if (first_time) then
     first_time = .false.
     rewind (io%input); read (io%input, initial_params)
     if (io%master) write (*, initial_params)
   end if
   self%b0 = b0
   self%gamma = gamma
   self%v_head = v_head_cgs/scaling%u
   self%v_accrete = v_accrete_cgs/scaling%u
   self%selected => planet_atmosphere
   call trace_end
END SUBROUTINE init

!> =============================================================================
!> Interpolate in a hydrostatic planet atmosphere, computed in atmosphere_mod
!> =============================================================================
SUBROUTINE planet_atmosphere (self, mem, m, omega)
  class(initial_t)             :: self
  real, dimension(:,:,:,:,:,:) :: mem
  class(mesh_t), allocatable, dimension(:):: m
  real:: omega(3), gamma
  !.............................................................................
  integer :: ix, iy, iz, nr
  real(8) :: x, y, z, r, d, d0, p, csound2
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t%planet_atmosphere')
  call self%init
  csound2 = self%gamma*atmosphere%pg(1)/atmosphere%rho(1)
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
        csound2 = self%gamma*p/d
        !-----------------------------------------------------------------------
        ! Vertical stratification is added by adding Gaussian profile and
        ! subtracting off the otherwise constant external density
        !-----------------------------------------------------------------------
        d0 = d
        d = d + &
          atmosphere%external%density*(exp(-(omega(3)**2/csound2*0.5*z**2))-1d0)
        p = p*(d/d0)
        mem(ix,iy,iz,1,1,1) = d
        mem(ix,iy,iz,2,1,1) = d*(log(p)-self%gamma*log(d))/(self%gamma-1.0)
        !-----------------------------------------------------------------------
        ! A rotation profile, where the rotation center is at y=+1 AU, and
        ! the direction of rotation is +x, so the head wind is in the -x
        ! direction.  The slower rotation at larger radii means that Ux is
        ! proportional to +y.  This gives a speed difference of only 130 cm/s
        ! across the Earth, and plus/minus 0.15 km/s across two Hill radii.
        ! The effect is not enough to significantly flatten the outer part of
        ! the Hill sphere, so we will probably only be seeing effects of the
        ! Coriolis force beyond the Hill radius.
        !-----------------------------------------------------------------------
        mem(ix,iy,iz,3,1,1) = d*(0.5*omega(3)*y + self%v_head)  ! azimuthal speed
        mem(ix,iy,iz,4,1,1) = d*self%v_accrete                  ! accretion speed
        mem(ix,iy,iz,5,1,1) = 0.0                               ! vertical speed
        mem(ix,iy,iz,6,1,1) = self%b0(1)
        mem(ix,iy,iz,7,1,1) = self%b0(2)
        mem(ix,iy,iz,8,1,1) = self%b0(3)
      end do
    end do
  end do
  call trace_end
END SUBROUTINE planet_atmosphere

!> =============================================================================
!> Use default Lagrange interpolation
!> =============================================================================
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
    f = interp%sequence (x0, x, y)
  end if
END FUNCTION

END MODULE
