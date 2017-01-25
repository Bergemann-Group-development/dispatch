!> =============================================================================
MODULE parameters_mod
  USE io_mod
  USE trace_mod
  USE scaling_mod
  implicit none
  private

  type, public :: experiment_parameters_t
    real(8) :: mtotal  = 1.0     ! mass of star + planet
    real(8) :: aplanet = 1.0     ! fixed orbital radius of the planet
    real(8) :: mplanet = 1.0e-3  ! mass of the planet
    real(8) :: gpp_soften = 0.6  ! softening coefficient for the planets grav. pot.
    real(8) :: mstar = 1.0 - 1.0e-3  ! mass of the star
    real(4) :: kill_rmin = 0.5   ! inner and outer radii where one starts to
    real(4) :: kill_rmax = 2.1   ! apply the "wave killing" boundary conditions
    real(8) :: aspect_rat = 0.05 ! disk aspect ratio (= H/R)
    logical :: ramp_mass = .true. ! increase the mass of the planet from zero to mplanet slowly
    real(8) :: phiplanet = 0.0   ! initial phi position of the planet
    real(8) :: roffset = 0.4     ! offset of the inner edge of the domain from the (cylindrical) origin
  end type
  type(experiment_parameters_t),  public :: experiment_params

  ! the embedded planet
  ! (not that it has no physical size, thus the need for softening)
  type, public :: planet_t
    real(8) :: massratio, mass, period, mplanet, softening
    real(8), dimension(3) :: position, velocity
  contains
    procedure :: init => init_planet
    procedure :: grow => grow_planet
    procedure :: move => move_planet
  end type
  type(planet_t), public :: planet

  real(8), parameter :: pi = 3.14159265358979323846d0

CONTAINS
!> =============================================================================
SUBROUTINE init_planet (self)
  class(planet_t) :: self
  !-----------------------------------------------------------------------------
  call trace_begin('planet_t%init')

  self%position = [0.0d0, experiment_params%aplanet, experiment_params%phiplanet]
  self%softening = experiment_params%gpp_soften
  self%massratio = experiment_params%mplanet / experiment_params%mtotal
  self%period = 2.0 * pi * sqrt(experiment_params%aplanet**3 &
              / (scaling%grav * experiment_params%mstar))
  self%velocity = [0.0d0, 0.0d0, 2.0d0 * pi / self%period]

  ! set the initial mass of the planet (at t = 0), to zero.
  ! the planet then grows using self%grow (eq. 9 of de Val-Borro et al.)
  self%mass = 0.0
  self%mplanet = experiment_params%mplanet

  call trace_end
END SUBROUTINE init_planet

!> =============================================================================
!> If `ramp_mass` is true, grow the planet mass from zero to `mplanet` over
!> roughly 100 orbits; this is a standard prescription and taken from
!> de Val-Borro et al. (2004).
!> =============================================================================
SUBROUTINE grow_planet (self, time)
  class(planet_t) :: self
  real(8), intent(in) :: time
  !-----------------------------------------------------------------------------
  if (experiment_params%ramp_mass) then
    self%mass = self%mplanet * (sin(0.1d0 * pi * time / self%period))**2
  else
    self%mass = self%mplanet
  end if
END SUBROUTINE grow_planet

!> =============================================================================
!> Move the planet in its orbit (assuming it is perfectly circular and it does
!> not migrate).
!> This subroutine also (unfortunately) assumes the phi domain is [-pi, pi].
!> =============================================================================
SUBROUTINE move_planet (self, dt, omega)
  class(planet_t) :: self
  real(8), intent(in) :: dt, omega
  !-----------------------------------------------------------------------------
  associate (phipos => self%position(3))
  phipos = phipos + (self%velocity(3) - omega) * dt
  if (phipos >  pi) phipos = phipos - 2.0d0 * pi
  if (phipos < -pi) phipos = phipos + 2.0d0 * pi
  ! Account for machine-round off error (relative to omega).
  if (abs(phipos) <= omega * epsilon(omega)) phipos = 0.0d0
  end associate

END SUBROUTINE move_planet

END MODULE parameters_mod
