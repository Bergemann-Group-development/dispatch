MODULE standard_disk_mod
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a standard disk
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
         mu=2.4, gamma=1.4
  namelist /disk_params/ sigma_cgs, sigma_p, temperature, temperature_p, mu, gamma
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
END SUBROUTINE init

!> =============================================================================
!> Return gas velocity and drag coefficient (inverse stopping time).
!> =============================================================================
SUBROUTINE drag (self, position, v, v_gas, c)
  class(standard_disk_t):: self
  real, dimension(3):: position, v, v_gas
  real:: c
  real, parameter:: density_solid_cgs = 3.
  real, parameter:: particle_size_cgs = 0.1
  real:: rm, sigma, temp, cs, h, density, v_kepler, v_g, dv2, aver
  real:: density_solid, particle_size
  !.............................................................................
  density_solid = density_solid_cgs/scaling%d
  particle_size = particle_size_cgs/scaling%l
  aver = 0.0
  rm = sqrt(position(1)**2+position(2)**2+position(3)**2)
  sigma = self%sigma*rm**self%sigma_p
  print*,sigma,self%sigma,rm
  temp = self%temp*rm**self%temp_p
  cs = sqrt(self%gamma*cgs%k_b*temp/(self%mu*cgs%m_u))/scaling%u
  v_kepler = sqrt(star%g*star%m/rm)
  h = rm*cs/v_kepler
  v_g = v_kepler*(1.-(cs/v_kepler)**2)
  v_gas (1) =  v_g*position(2)/rm
  v_gas (2) = -v_g*position(1)/rm
  v_gas (3) =  0.0
  dv2 = (v_gas(1)-v(1))**2+(v_gas(2)-v(2))**2+(v_gas(3)-v(3))**2
  density = sigma/h*exp(-(position(3)/h)**2)
  c = (density/density_solid)*sqrt(cs**2+dv2)/particle_size
  print*,sigma,h,density,position(3)/h,c
  aver = aver + dv2/v_kepler**2
END SUBROUTINE drag

END MODULE standard_disk_mod
