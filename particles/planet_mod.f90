!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for planet task
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE planet_mod
  USE io_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  USE task_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  type, public:: planet_t
    real(8):: mass, radius, temperature
    real(8), dimension(3):: position
    real(8):: t_growth= 0d0
    real:: orbit_radius, hill_radius
    real:: bondi_radius = 0.0
    real:: omega(3)=0.0
    real:: friction=0.0
    real:: trap=0.0
    real:: v_head=0.0, v_accrete=0.0
    real:: r_smooth=0.0
    contains
    procedure:: init
    procedure:: inside
    procedure:: newton
  end type
  type(planet_t), public:: planet
CONTAINS

!===============================================================================
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
  real:: r_smooth=0.0
  real(8):: t_growth=0d0                ! in orbit periods
  real(8):: t_inj = 0.5
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /planet_params/ orbit_radius_au, rad_mult, position_orbit_radii, radius_r_earth, &
    mass_m_earth, temperature_k, artificial_friction, trap, r_smooth
  !-----------------------------------------------------------------------------
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, planet_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, planet_params)
    end if
    write (io%namelist,planet_params)
  end if
  call io%lock%unset
  !
  self%radius = radius_r_earth*rad_mult*cgs%r_earth/scaling%l
  self%mass = mass_m_earth*cgs%m_earth/scaling%m
  self%orbit_radius = orbit_radius_au*cgs%au/scaling%l
  self%position = position_orbit_radii*self%orbit_radius
  self%temperature = temperature_k/scaling%temp
  self%friction = artificial_friction
  self%trap = trap
  self%r_smooth = r_smooth
  self%t_growth = t_growth * self%orbit_radius * 2 * cgs%pi
  !print '(a,3g14.5)', ' planet postion:', self%position
END SUBROUTINE init


!===============================================================================
FUNCTION inside (self, p) RESULT (out)
  class (planet_t):: self
  real:: p(3)
  logical:: out
  !.............................................................................
  out = sqrt(sum(p**2)) < self%radius
END FUNCTION inside

FUNCTION newton (self, rad, sp) RESULT (out)
  class (planet_t):: self
  real(8) :: rad, sp
  real    :: out
  !.....................................................................
  ! --------------------------------------------------------------------
  if (rad > self%r_smooth) then
    out = self%mass*scaling%grav/rad**2
  else
    out = self%mass*scaling%grav/self%r_smooth**2*(4.*(rad/self%r_smooth)-3.*(rad/self%r_smooth)**2)
  end if
  if (sp >= 0) out = -out

END FUNCTION newton

!===============================================================================
END MODULE
