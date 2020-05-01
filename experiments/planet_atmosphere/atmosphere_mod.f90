!===============================================================================
!> Solve the hydrostatic equilibrium for an extended planet atmosphere, 
!> stretching from the planet surface to beyond the Hill radius.  The mass
!> of the atmosphere is neglected, because the Python examples show that most
!> of the mass is anyway at large radii, where the solution will not be exact.
!===============================================================================
MODULE atmosphere_mod
  USE io_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  USE particle_mod
  USE planet_mod
  implicit none
  private
  type, public:: conditions_t
    real:: density, pressure, temperature, entropy, gamma
  end type
  type, public, extends(particle_t):: atmosphere_t
    integer:: nr, niter
    real:: gamma(2), temp_gamma
    real(kind=8), dimension(:), pointer:: rho, pg, temp, rr
    real:: dlnr, hill_radius, orbit_radius
    type(conditions_t):: external, surface
  contains
    procedure init
    procedure hydrostatic
  end type
  type(atmosphere_t), public:: atmosphere
CONTAINS

!===============================================================================
!> Read parameters
!===============================================================================
SUBROUTINE init (self)
  class(atmosphere_t):: self
  logical, save:: first_time=.true.
  real:: density_cgs=1e-10, temperature_K=200, delta_r=0.01, orbit_au=1.0
  real:: mass_m_earth=0.9
  real:: gamma(2)=1.4, temp_gamma=2000.
  integer:: nr=500, niter=10
  namelist /atmosphere_params/ orbit_au, density_cgs, temperature_K, &
    mass_m_earth, nr, niter, gamma, temp_gamma
  !----------------------------------------------------------------------------
  call trace_begin('atmosphere_t%init')
  call scaling%init
  !$omp critical (atm_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input); read (io%input, atmosphere_params)
    if (io%master) write (*, atmosphere_params)
  end if
  !$omp end critical (atm_cr)
  self%nr = nr
  self%niter = niter
  self%gamma = gamma
  self%temp_gamma = temp_gamma
  allocate (self%rr(self%nr))
  allocate (self%rho(self%nr))
  allocate (self%pg(self%nr))
  allocate (self%temp(self%nr))
  self%orbit_radius = orbit_au*cgs%au/scaling%l
  self%external%density = density_cgs/scaling%d
  self%external%temperature = temperature_K/scaling%temp
  self%external%pressure = density_cgs*cgs%k_B*temperature_K/(scaling%mu*cgs%m_u) &
    /scaling%p
  self%hill_radius = self%orbit_radius*(planet%mass*scaling%m/(3d0*cgs%m_sun))**(1d0/3d0)
  planet%radius = cgs%r_earth/scaling%l*(planet%mass*scaling%m/cgs%m_earth)**(1d0/3d0)
  self%dlnr = log(self%hill_radius/planet%radius)/(self%nr-1.0)
  if (io%verbose>0) print*, 'dlnr:', planet%radius, self%hill_radius/planet%radius, self%dlnr
  call self%hydrostatic
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Solve the equations for hydrostatic balance, initially neglecting the mass
!> of the pseudo-atmosphere.
!===============================================================================
SUBROUTINE hydrostatic (self)
  class(atmosphere_t):: self
  logical, save:: first_time=.true.
  real(8):: dr, dlnr, dlnpg, temp_K, dlnT, gamma, s
  integer:: nr, i, iter
  !-----------------------------------------------------------------------------
  ! Start with the external pressure and integrate inwards
  !-----------------------------------------------------------------------------
  call io%lock%set
  if (first_time) then
    first_time = .false.
    call trace_begin('experiment_t%hydrostatic_balance')
    nr = self%nr
    self%temp(nr) = self%external%temperature
    self%rho(nr)  = self%external%density
    self%pg(nr)   = self%external%pressure
    self%rr(nr)   = self%hill_radius
    if (io%verbose>1) print*, &
      '              i         rr         pg       temp        rho         dr      dlnpg       grav'
    do i=nr-1,1,-1
      self%rr(i) = exp(log(planet%radius) + (i-1)*log(self%hill_radius/planet%radius)/(nr-1))
      self%temp(i) = self%temp(i+1)
      self%rho(i) = self%rho(i+1)
      self%pg(i) = self%pg(i+1)
      self%rr(i) = self%rr(i+1)*exp(-self%dlnr)
      do iter=1,self%niter
        self%rho(i) = self%pg(i)/self%temp(i)
        dr = self%rr(i+1) - self%rr(i)
        dlnpg = dr*code%grav*planet%mass* &
          0.5*(self%rho(i+1)/self%pg(i+1)/self%rr(i+1)**2 + &
               self%rho(i  )/self%pg(i  )/self%rr(i  )**2)
        temp_K = self%temp(i)*scaling%temp
        if (temp_K > self%temp_gamma) then
          gamma = self%gamma(2)
        else
          gamma = self%gamma(1)
        end if
        dlnT = dlnpg*(gamma-1.0)/gamma
        self%pg  (i) = self%pg  (i+1)*exp(dlnpg)
        self%temp(i) = self%temp(i+1)*exp(dlnT)
        if (io%verbose>2) print*, i, iter, self%rr(i), self%pg(i), self%rho(i), self%temp(i)
      end do
      if (self%gamma(2)==1.0) then
        s = 1.0
      else
        s = (log(self%pg(i))-gamma*log(self%rho(i)))/(self%gamma(2)-1.0)
      end if
      if (io%verbose>1) print'("hydrostatic",i5,1p,9e11.3)', &
        i,self%rr(i)*scaling%l/cgs%r_earth, self%pg(i)*scaling%p, self%temp(i)*scaling%temp, &
        self%rho(i)*scaling%d, s, dlnpg, planet%mass*code%grav/self%rr(i)**2
    end do
    if (io%verbose>1) print*, &
      '              i         rr         pg       temp        rho      dlnpg       dlnT       grav'
    self%surface%density     = self%rho(1)
    self%surface%pressure    = self%pg(1)
    self%surface%temperature = self%temp(1)
    print '(a,f7.0," K")', 'atmosphere_mod::hydrostatic: surface%temperature =', &
      self%surface%temperature*scaling%temp
    call trace_end
  end if
  call io%lock%unset
END SUBROUTINE hydrostatic

END MODULE atmosphere_mod
