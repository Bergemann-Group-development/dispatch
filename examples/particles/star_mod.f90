!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Simple star module
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE star_mod
  USE io_mod
  USE particle_mod
  implicit none
  private
  type, public, extends(particle_t):: star_t
    real:: g
    real:: omega(3)
  contains
    procedure:: init
    procedure:: kepler_speed
    procedure:: central_force
  end type
  integer:: nprint=10
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self)
  class(star_t):: self
  real:: grav=1.0, mass=1.0
  real:: omega(3) = [.0, 0.0, 1.0]
  logical, save:: first_time=.true.
  namelist /star_params/ grav, mass, omega
  !.............................................................................
  call io%init
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, star_params)
    if (io%master) write (*, star_params)
  end if
  self%m = mass
  self%g = grav
  self%omega = omega
END SUBROUTINE init

!> =============================================================================
!> =============================================================================
FUNCTION cross (a, b) RESULT (out)
  real, intent(in):: a(3), b(3)
  real:: out(3)
  !.............................................................................
  out(1) = a(2)*b(3)-a(3)*b(2)
  out(2) = a(3)*b(1)-a(1)*b(3)
  out(3) = a(1)*b(2)-a(2)*b(1)
END FUNCTION cross

!> =============================================================================
!> =============================================================================
FUNCTION kepler_speed (self, r) RESULT (out)
  class(star_t):: self
  real, intent(in):: r(3)
  real:: out(3)
  real:: rm, v
  !.............................................................................
  rm = sqrt(sum(r**2))
  v = self%g*self%m/rm
  out = cross (r, self%omega)
  out = v*out/sqrt(sum(out**2))
  if (nprint>0) then
    print*,'kepler_speed',out
    nprint = nprint-1
  end if
END FUNCTION kepler_speed

!> =============================================================================
!> =============================================================================
FUNCTION central_force (self, r) RESULT (out)
  class(star_t):: self
  real, intent(in):: r(3)
  real:: out(3)
  !.............................................................................
  out = -r*(self%g*self%m/sqrt(sum(r**2))**3)
END FUNCTION central_force

END MODULE star_mod
