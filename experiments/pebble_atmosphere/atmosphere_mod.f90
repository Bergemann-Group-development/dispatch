!> =============================================================================
!> Solve the hydrostatic equilibrium for an extended planet atmosphere, 
!> stretching from the planet surface to beyond the Hill radius.  The mass
!> of the atmosphere is neglected, because the Python examples show that most
!> of the mass is anyway at large radii, where the solution will not be exact.
!> =============================================================================
MODULE atmosphere_mod
  USE io_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  USE particle_mod
  USE planet_mod
  USE lagrange_mod
  implicit none
  private
  type, public:: conditions_t
    real:: density, pressure, temperature, energy, entropy, gamma, tau
  end type
  type, public, extends(particle_t):: atmosphere_t
    integer:: nr, niter
    real:: gamma(2), temp_gamma
    real:: mass
    real:: t_heat
    real(kind=8), dimension(:), pointer:: rho, pg, temp, rr, tau
    real(8):: dlnr, hill_radius, orbit_radius, kappa, L
    logical :: do_rad
    type(conditions_t):: external, surface
  contains
    procedure init
    procedure hydrostatic
  end type
  integer, save:: verbose = 0
  type(atmosphere_t), public:: atmosphere
  real(8), parameter:: pi = asin(1.0d0) * 2.0d0
CONTAINS

!> =============================================================================
!> Read parameters
!> =============================================================================
SUBROUTINE init (self)
  class(atmosphere_t):: self
  logical, save:: first_time=.true., do_rad=.false.
  real,save:: density_cgs=1e-10, temperature_K=200, delta_r=0.01, orbit_au=1.0
  real,save:: t_heat = 1.0
  real,save:: gamma(2)=1.4, temp_gamma=2000., kappa_cgs=0.1, m_acc_me=1e-5
  real(8):: inner_radius
  integer:: iostat
  integer, save:: nr=500, niter=10
  namelist /atmosphere_params/ orbit_au, density_cgs, temperature_K, &
    nr, niter, gamma, temp_gamma, kappa_cgs, m_acc_me, do_rad, &
    t_heat, verbose
  !-----------------------------------------------------------------------------
  call trace_begin('atmosphere_t%init')
  call scaling%init
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, atmosphere_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, atmosphere_params)
    end if
    write (io%namelist, atmosphere_params)
  end if
  call io%lock%unset
  self%nr = nr
  self%niter = niter
  self%gamma = gamma
  self%temp_gamma = temp_gamma
  self%external%gamma = gamma(1)
  self%surface%gamma  = gamma(2)
  self%t_heat = t_heat
  self%mass = planet%mass
  self%L = scaling%grav*self%mass/(planet%radius*0.5)*m_acc_me/cgs%yr*scaling%t
  self%do_rad = do_rad
  allocate (self%rr(self%nr))
  allocate (self%rho(self%nr))
  allocate (self%pg(self%nr))
  allocate (self%temp(self%nr))
  allocate (self%tau(self%nr))
  self%kappa = kappa_cgs*scaling%m/scaling%l**2
  self%orbit_radius = orbit_au*cgs%au/scaling%l
  self%external%density = density_cgs/scaling%d
  self%external%temperature = temperature_K/scaling%temp
  self%external%pressure = &
    density_cgs*cgs%k_B*temperature_K/(scaling%mu*cgs%m_u)/scaling%p
  self%external%energy = self%external%pressure/self%external%density &
    /(self%external%gamma-1.0)
  self%hill_radius = self%orbit_radius&
                   * (planet%mass*scaling%m/(3d0*cgs%m_sun))**(1d0/3d0)
  planet%hill_radius = self%hill_radius
  !-----------------------------------------------------------------------------
  ! Make sure the pressure force is contiuous at the surface, by continuing the
  ! integration to half the planet radius
  !-----------------------------------------------------------------------------
  !planet%radius = cgs%r_earth/scaling%l*(planet%mass*scaling%m/cgs%m_earth)**(1d0/3d0)
  inner_radius = planet%radius*0.5
  self%dlnr = log(self%hill_radius/inner_radius)/(self%nr-1.0)
  if (verbose > 0) print*, 'dlnr:', planet%radius, self%hill_radius/planet%radius, self%dlnr
  call self%hydrostatic
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Solve the equations for hydrostatic balance, initially neglecting the mass
!> of the pseudo-atmosphere.
!> =============================================================================
SUBROUTINE hydrostatic (self)
  class(atmosphere_t):: self
  logical, save:: first_time=.true.
  real(8):: dr, dlnr, dlnpg, temp_K, dlnT, gamma, s, dlnT_rad
  real:: rb
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
    self%tau(nr)  = 0.0
    if (verbose > 1) print*, &
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
        dlnpg = dr*scaling%grav*self%mass* &
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
        if (self%do_rad) then
          dlnT_rad=3.*self%kappa*self%L*self%pg(i)/(64.*scaling%stefan*scaling%grav*self%mass*self%temp(i)**4)
          dlnT = min(dlnT,dlnT_rad)
        end if
        !print *, 'dlnt, dlnt_rad', dlnt, dlnt_rad,'exp rad', exp(dlnt_rad)
        self%temp(i) = self%temp(i+1)*exp(dlnT)
        if (verbose > 2) print*, i, iter, self%rr(i), self%pg(i), self%rho(i), self%temp(i)
      end do
      self%tau(i) = self%tau(i+1) + 0.5*(self%rho(i)+self%rho(i+1))*self%kappa*dr
      if (self%gamma(2)==1.0) then
        s = 1.0
      else
        s = (log(self%pg(i))-gamma*log(self%rho(i)))/(self%gamma(2)-1.0)
      end if
      rb = self%mass*scaling%grav/self%temp(i)
      if (rb > planet%bondi_radius) planet%bondi_radius = rb
      if (verbose > 1) print'("hydrostatic",i5,1p,10e11.3)', &
        i, self%rr(i)*scaling%l/cgs%r_earth, self%pg(i)*scaling%p, &
        self%temp(i)*scaling%temp, self%rho(i)*scaling%d, s, dlnpg, &
        self%mass*scaling%grav/self%rr(i)**2, &
        (log(self%pg(i))-self%gamma*log(self%rho(i)))/(self%gamma-1.0)
    end do
    !---------------------------------------------------------------------------
    ! Add a another layer of optical depth at the outer boundary, as an estimate
    !---------------------------------------------------------------------------
    self%tau = self%tau + self%tau(nr-1)
    if (verbose > 0) then
      print*, &
        '              i     rr         pg       temp        rho          s       grav        tau'
      do i=1,nr
        if (self%gamma(2)==1.0) then
          s = 1.0
        else
          s = (log(self%pg(i))-gamma*log(self%rho(i)))/(self%gamma(2)-1.0)
        end if
        print'("hydrostatic",i5,1p,9g11.3)',i,self%rr(i)*scaling%l/cgs%r_earth, &
          self%pg(i)*scaling%p, self%temp(i)*scaling%temp,self%rho(i)*scaling%d, &
          s, self%mass*scaling%grav/self%rr(i)**2*scaling%u**2/scaling%l, &
          self%tau(i)
      end do
      print*, &
        '              i     rr         pg       temp        rho          s       grav        tau'
    end if
    self%surface%density     = lagrange%interpolate1d (planet%radius, self%rr, self%rho)
    self%surface%pressure    = lagrange%interpolate1d (planet%radius, self%rr, self%pg)
    self%surface%temperature = lagrange%interpolate1d (planet%radius, self%rr, self%temp)
    self%surface%tau         = lagrange%interpolate1d (planet%radius, self%rr, self%tau)
    print '(a,1p,g12.4," K")', ' atmosphere_t%hydrostatic: surface%temperature =', &
      self%surface%temperature*scaling%temp
    print '(a,1p,g12.4," g cm-3")', ' atmosphere_t%hydrostatic: surface%density =', &
      self%surface%density*scaling%d
    call trace_end
  end if
  call io%lock%unset
END SUBROUTINE hydrostatic

END MODULE atmosphere_mod
