!> =============================================================================
!> =============================================================================
MODULE planetesimal_mod
  USE io_mod
  USE trace_mod
  USE bits_mod
  USE patch_mod
  USE mhd_mod
  USE mesh_mod
  USE units_mod
  USE scaling_mod
  implicit none
  private

  ! an immersed boundary that acts like a solid sphere (our planetesimal)
  type, public :: planetesimal_t
    real(8) :: radius, mass, vol_density, period, massratio, mass_final
    real(8) :: mstar ! alias for the same variable in `experiment_params`
    real(8) :: cyl_pos(3) ! position of the planet in cylindrical coords.
    real(8) :: position(3) ! position in Cartesian coords.
    real(8) :: velocity(3)
    logical :: ramp_mass
  contains
    procedure :: init => init_planetesimal
    procedure :: grow => grow_planetesimal
    procedure :: move => move_planetesimal
    procedure :: update => update_planetesimal
  end type

  real(8), parameter :: pi = asin(1.0) * 2.0

CONTAINS

!> =============================================================================
!> Initialize a solid sphere with a given radius and position
!> =============================================================================
SUBROUTINE init_planetesimal (self)
  class(planetesimal_t) :: self
  !-----------------------------------------------------------------------------
  ! Namelist parameters for the component
  !-----------------------------------------------------------------------------
  real(8) :: radius, aplanet, phiplanet, zplanet, mass, mass_final
  logical :: ramp_mass = .true.
  namelist / planetesimal_params / radius, mass, mass_final, ramp_mass, aplanet, &
                                   phiplanet, zplanet
  !-----------------------------------------------------------------------------
  call trace_begin ('solid_sphere_t%init')

  ! default values
  radius = cgs%r_earth / cgs%AU ! radius in AU
  aplanet = 1.0 ! semi-major axis in units of Jupiter's semi-major axis
  phiplanet = 0.0
  zplanet = 0.0
  mass = 0.0 ! initial mass in Solar mass units
  mass_final = 3.0e-6 ! (final) Earth mass in Solar mass units

  rewind (io%input); read (io%input, planetesimal_params)
  if (io%master) write (*, planetesimal_params)

  self%cyl_pos(1) = zplanet
  self%cyl_pos(2) = aplanet
  self%cyl_pos(3) = phiplanet
  self%position = CylindricalToCartesian(zplanet, aplanet, phiplanet)
  self%radius = radius
  self%mass_final = mass_final
  if (.not. ramp_mass) then
    self%mass = mass_final
  else
    self%mass = mass
  end if
  self%vol_density = self%mass * 0.75 / pi / self%radius**3
  self%period = 2.0 * pi * sqrt(aplanet**3 &
              / (scaling%grav * self%mstar))
  self%velocity = [0.0d0, 0.0d0, sqrt(scaling%grav * self%mstar / aplanet**3)]
  self%massratio = mass_final + (self%mstar + mass_final)

  call trace_end
END SUBROUTINE init_planetesimal

!> =============================================================================
!> If `ramp_mass` is true, grow the planet mass from zero to `mplanet` over
!> roughly 100 orbits; this is a standard prescription and taken from
!> de Val-Borro et al. (2004).
!> =============================================================================
SUBROUTINE grow_planetesimal (self, time)
  class(planetesimal_t) :: self
  real(8), intent(in) :: time
  !-----------------------------------------------------------------------------

  if (self%ramp_mass) then
    self%mass = self%mass_final * (sin(0.1d0 * pi * time / self%period))**2
    self%radius = self%mass * 0.75 / pi / self%vol_density
  else
    self%mass = self%mass_final
  end if

END SUBROUTINE grow_planetesimal

!> =============================================================================
!> Move the planet in its orbit (assuming it is perfectly circular and it does
!> not migrate).
!> This subroutine also (unfortunately) assumes the phi domain is [-pi, pi].
!> =============================================================================
SUBROUTINE move_planetesimal (self, dt, omega)
  class(planetesimal_t) :: self
  real(8), intent(in) :: dt, omega
  !-----------------------------------------------------------------------------

  associate (phipos => self%cyl_pos(3))
  phipos = phipos + (self%velocity(3) - omega) * dt
  if (phipos >  pi) phipos = phipos - 2.0d0 * pi
  if (phipos < -pi) phipos = phipos + 2.0d0 * pi
  ! Account for machine-round off error (relative to omega).
  if (abs(phipos) <= omega * epsilon(omega)) phipos = 0.0d0
  self%position = CylindricalToCartesian(self%cyl_pos(1), self%cyl_pos(2), self%cyl_pos(3))
  end associate

END SUBROUTINE move_planetesimal

!> =============================================================================
!> =============================================================================
SUBROUTINE update_planetesimal (self, time, dt)
  class(planetesimal_t):: self
  real(8), intent(in):: time, dt
  !-----------------------------------------------------------------------------
  call trace_begin('planetesimal_t%update')

  ! update the planet mass
  call self%grow(time)
  ! update the planet position
  call self%move(dt, self%velocity(3))

  call trace_end
END SUBROUTINE update_planetesimal

END MODULE planetesimal_mod
