MODULE standard_disk_mod
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a standard disk, with a power law radial profile in surface density
!> and in temperature.  In addition it has the realistic feature that the gas
!> spirals in, with a constant factor alpha between the radial and azimuthal
!> speeds.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE trace_mod
  USE star_mod
  USE units_mod
  USE scaling_mod
  implicit none
  private
  type, public:: standard_disk_t
    real:: sigma                                ! surface mass density
    real:: sigma_p                              ! power exponent
    real:: temp                                 ! temperature normalization
    real:: temp_p                               ! power exponent
    real:: mu                                   ! mean molecular mass
    real:: gamma                                ! gas gamma
    real:: alpha                                ! spiraling factor
    real:: particle_size                        ! particle size in code units
    real:: particle_density                     ! particle solid density in code units
  contains
    procedure:: init
    procedure:: drag
  end type
  type(standard_disk_t), public:: standard_disk
CONTAINS

!> =============================================================================
!> Read parameters for a standard disk.
!> =============================================================================
SUBROUTINE init (self)
  class(standard_disk_t):: self
  real:: sigma_cgs=1700., sigma_p=-1.5, temperature=100., temperature_p=-1.0, &
         mu=2.4, gamma=1.4, alpha=0.001
  real:: particle_density_cgs = 2.
  real:: particle_size_cgs = 1.
  namelist /disk_params/ sigma_cgs, sigma_p, temperature, temperature_p, mu, &
    gamma, alpha, particle_size_cgs, particle_density_cgs
  !.............................................................................
  rewind (io%input)
  read (io%input, disk_params)
  if (io%master) write (*, disk_params)
  self%sigma = sigma_cgs*scaling%l**2/scaling%m
  self%sigma_p = sigma_p
 !self%temp = temperature/scaling%temp
  self%temp = temperature
  self%temp_p = temperature_p
  self%gamma = gamma
  self%mu = mu
  self%alpha = alpha
  self%particle_size = particle_size_cgs/scaling%l
  self%particle_density = particle_density_cgs/scaling%d
END SUBROUTINE init

!> =============================================================================
!> Return gas velocity and drag coefficient (inverse stopping time).
!> =============================================================================
SUBROUTINE drag (self, r, v, v_gas, v_kep, c)
  class(standard_disk_t):: self
  real, dimension(3):: r, v, v_gas, v_kep
  real:: c
  real:: rm, sigma, temp, cs, h, density, v_kepler, v_g, dv2, aver
  !.............................................................................
  aver = 0.0
  rm = sqrt(r(1)**2+r(2)**2+r(3)**2)
  sigma = self%sigma*rm**self%sigma_p
  temp = self%temp*rm**self%temp_p
  cs = sqrt(self%gamma*cgs%k_b*temp/(self%mu*cgs%m_u))/scaling%u
  v_kepler = sqrt(star%g*star%m/rm)
  h = rm*cs/v_kepler
  v_kep(1) =  v_kepler*r(2)/rm
  v_kep(2) = -v_kepler*r(1)/rm
  v_kep(3) = 0.0
  v_g = -v_kepler*(cs/v_kepler)**2
  v_gas (1) = ( v_g*r(2) - self%alpha*v_kepler*r(1))/rm
  v_gas (2) = (-v_g*r(1) - self%alpha*v_kepler*r(2))/rm
  v_gas (3) =  0.0
  dv2 = (v_gas(1)-v(1))**2+(v_gas(2)-v(2))**2+(v_gas(3)-v(3))**2
  density = sigma/h*exp(-(r(3)/h)**2)
  c = (density/self%particle_density)*sqrt(cs**2+dv2)/self%particle_size
  aver = aver + dv2/v_kepler**2
END SUBROUTINE drag

END MODULE standard_disk_mod
