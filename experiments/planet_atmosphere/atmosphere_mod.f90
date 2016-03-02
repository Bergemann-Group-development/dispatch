!> =============================================================================
!> =============================================================================
MODULE atmosphere_mod
  USE io_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  USE particle_mod
  implicit none
  private
  type, public:: conditions_t
    real:: density, pressure, temperature, entropy, gamma
  end type
  type, public, extends(particle_t):: atmosphere_t
    integer:: nr
    real(kind=8), dimension(:), pointer:: rho, pg, temp, rr
    real:: planet_mass, planet_radius, dlnr, hill_radius, orbit_radius
    type(conditions_t):: external, surface
  contains
    procedure hydrostatic
  end type
  type(atmosphere_t), public:: atmosphere
CONTAINS

!> =============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!> =============================================================================
SUBROUTINE init (self)
  class(atmosphere_t):: self
  real:: density_cgs=1e-10, temperature_K=200, delta_r=0.01, orbit_au=1.0
  real:: mass_m_earth=0.9
  integer:: nr=500
  logical, save:: first_time=.true.
  namelist /atmosphere_params/ orbit_au, density_cgs, temperature_K, nr
  !----------------------------------------------------------------------------
  call trace_begin('atmosphere_t%init')
  call scaling%init
  if (first_time) then
    first_time = .false.
    rewind (io%input); read (io%input, atmosphere_params)
    if (io%master) write (*, atmosphere_params)
  end if
  self%nr = nr
  allocate (self%rr(self%nr))
  allocate (self%rho(self%nr))
  allocate (self%pg(self%nr))
  allocate (self%temp(self%nr))
  self%orbit_radius = orbit_au*cgs%au/scaling%l
  self%planet_mass = mass_m_earth*cgs%m_earth/scaling%m
  self%external%density = density_cgs/scaling%d
  self%external%temperature = temperature_K/scaling%temp
  self%external%pressure = density_cgs*cgs%k_B*temperature_K/(cgs%mu*cgs%m_u) &
    /scaling%p
  self%hill_radius = self%orbit_radius*(self%planet_mass*scaling%m/(3.*cgs%m_sun))**0.33333
  self%dlnr = log(self%hill_radius/self%planet_radius)/(self%nr-1.0)
  call self%hydrostatic
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Solve the equations for hydrostatic balance, initially neglecting the mass
!> of the pseudo-atmosphere -- later on we will include it.
!> =============================================================================
SUBROUTINE hydrostatic (self)
  class(atmosphere_t):: self
  real:: dlnr, dlnpg, temp_K, dlnT, dr, gamma
  integer:: nr, i, iter
  integer, parameter:: niter=10
  !-----------------------------------------------------------------------------
  ! Start with the external pressure and integrate inwards
  !-----------------------------------------------------------------------------
  call trace_begin('experiment_t%hydrostatic_balance')
  nr = self%nr
  self%temp(nr) = self%external%temperature
  self%rho(nr)  = self%external%density
  self%pg(nr)   = self%external%pressure
  self%rr(nr)   = self%hill_radius
  do i=nr-1,1,-1
    self%pg(i) = self%pg(i+1)
    do iter=1,niter
      self%rr(i) = self%rr(i+1)*exp(-self%dlnr)
      dr = self%rr(i+1)-self%rr(i)
      self%rho(i) = self%pg(i)/self%temp(i)
      dlnpg = dr*code%grav*self%planet_mass* &
        0.5*(self%rho(i+1)/self%pg(i+1)/self%rr(i+1)**2 + &
             self%rho(i  )/self%pg(i  )/self%rr(i  )**2)
      temp_K = self%temp(i+1)*scaling%temp
      if (temp_K > 2000.) then
        gamma = 5./3.
      else
        gamma = 7./5.
      end if
      dlnT = dlnpg*(gamma-1.0)/gamma
      self%pg  (i) = self%pg  (i+1)*exp(dlnpg)
      self%temp(i) = self%temp(i+1)*exp(dlnT)
    end do
    print'("hydrostatic",i5,1p,8e11.3)', &
      i,self%rr(i),self%pg(i),self%temp(i)*self%temp,self%rho(i), &
      dr,dlnpg,self%planet_mass*code%grav/self%rr(i)**2
  end do
  self%surface%density     = self%rho(1)
  self%surface%pressure    = self%pg(1)
  self%surface%temperature = self%temp(1)
  call trace_end
END SUBROUTINE hydrostatic

END MODULE
