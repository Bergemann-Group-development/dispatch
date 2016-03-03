!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for star task
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE star_mod
  USE io_mod
  USE task_mod
  USE particle_mod
  USE units_mod
  USE scaling_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  type, public, extends(particles_t):: star_t
    real:: mass, radius, surface_temperature
  contains
    procedure:: init
  end type
  type(star_t):: star
CONTAINS

!> =============================================================================
SUBROUTINE init (self)
  class(star_t):: self
  real:: mass_solar=1.0, radius_solar=1.0, surface_temperature=5900.0
  namelist /star_params/ mass_solar, radius_solar, surface_temperature
  !.............................................................................
  rewind (io%input); read (io%input, star_params); write (*, star_params)
  self%mass   =   mass_solar*cgs%m_sun/scaling%m
  self%radius = radius_solar*cgs%r_sun/scaling%l
  call self%task_t%init
END SUBROUTINE init

!> =============================================================================
END MODULE
