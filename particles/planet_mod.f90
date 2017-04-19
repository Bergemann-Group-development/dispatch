!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for planet task
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE planet_mod
  USE io_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  USE task_mod
  USE particle_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  type, public, extends(particles_t):: planet_t
    real(8):: mass, radius, temperature
    real:: orbit_radius, hill_radius
    real:: omega(3)=0.0
    real:: friction=0.0
    real:: trap=0.0
    real:: v_head=0.0, v_accrete=0.0
    contains
    procedure:: init
    procedure:: inside
  end type
  type(planet_t), public:: planet
CONTAINS

!> =============================================================================
SUBROUTINE init (self)
  class(planet_t):: self
  !.............................................................................
  real:: position_orbit_radii(3)=[0.0, 1.0, 0.0]
  real:: radius_r_earth=1.0
  real:: rad_mult = 1.0
  real:: mass_m_earth=1.0
  real:: temperature_k=5500.
  real:: orbit_radius_au=1.0
  real:: artificial_friction=0.0
  real:: trap=0.0
  real(8):: t_inj = 0.5
  namelist /planet_params/ orbit_radius_au, rad_mult, position_orbit_radii, radius_r_earth, &
    mass_m_earth, temperature_k, artificial_friction, trap
  !-----------------------------------------------------------------------------
  rewind (io%input); read (io%input, planet_params); write (*,planet_params)
  self%radius = radius_r_earth*rad_mult*cgs%r_earth/scaling%l
  self%mass = mass_m_earth*cgs%m_earth/scaling%m
  self%orbit_radius = orbit_radius_au*cgs%au/scaling%l
  self%position = position_orbit_radii*self%orbit_radius
  self%temperature = temperature_k/scaling%temp
  self%friction = artificial_friction
  self%trap = trap
  print *, "============= PLANET INIT POSITION ===========", self%position
END SUBROUTINE init


!> =============================================================================
FUNCTION inside (self, p) RESULT (out)
  class (planet_t):: self
  real:: p(3)
  logical:: out
  !.............................................................................
  out = sqrt(sum(p**2)) < self%radius
END FUNCTION inside

!> =============================================================================
END MODULE
