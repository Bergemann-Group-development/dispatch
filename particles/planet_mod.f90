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
  real(8):: position_au(3)=[0.0, 1.0, 0.0], radius_r_earth=1.0, mass_m_earth=1.0
  real(8):: temperature_k=5500.
  namelist /planet_params/ position_au, radius_r_earth, mass_m_earth, temperature_k
  !-----------------------------------------------------------------------------
  rewind (io%input); read (io%input, planet_params); write (*,planet_params)
  self%radius = radius_r_earth*cgs%r_earth/scaling%l
  self%mass = mass_m_earth*cgs%m_earth/scaling%m
  self%position = position_au*cgs%au/scaling%l
  self%temperature = temperature_k/scaling%temp
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
