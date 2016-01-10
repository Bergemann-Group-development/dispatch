MODULE standard_disk_mod
  USE io_mod
  USE trace_mod
  USE star_mod
  USE units_mod
  USE scaling_mod
  implicit none
  private
  type, public:: standard_disk_t
    real:: sigma
    real:: sigma_p
    real:: temp
    real:: temp_p
    real:: mu
    real:: gamma
  contains
    procedure:: init
    procedure:: structure
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
  self%sigma = sigma_cgs/scaling%l**2
  self%sigma_p = sigma_p
  self%temp = temperature/scaling%temp
  self%temp_p = temperature_p
  self%gamma = gamma
  self%mu = mu
END SUBROUTINE init

!> =============================================================================
!> Return density, temperature, and rotation velocity for a standard disk.  The
!> sub-keplerian roation speed is an approximation, correct to order of magnitude.
!> =============================================================================
SUBROUTINE structure (self, position, density, temperature, velocity)
  class(standard_disk_t):: self
  real, dimension(:,:):: position, velocity
  real, dimension(:):: density, temperature
  integer:: i, nr
  real:: rm, v_kepler, v_gas, sigma, temp, cs, h
  !.............................................................................
  nr = size(position,dim=2)
  do i=1,nr
    rm = sqrt(position(i,1)**2+position(i,2)**2+position(i,3)**2)
    v_kepler = sqrt(star%g*star%m/rm)
    sigma = self%sigma*rm**self%sigma_p
    temp = self%temp*rm**self%temp_p
    cs = sqrt(self%gamma*cgs%k_b*temp/(self%mu*cgs%m_u))/scaling%u
    h = rm*cs/v_kepler
    density(i) = sigma/h*exp(-(position(i,3)/h)**2)
    v_gas = v_kepler*(1.-(cs/v_kepler)**2)
    temperature(i) = temp
    velocity (i,1) =  v_gas*position(i,2)/rm
    velocity (i,2) = -v_gas*position(i,1)/rm
    velocity (i,3) =  0.0
  end do
END SUBROUTINE structure

!> =============================================================================
!> Return gas velocity and drag coefficient (inverse stopping time).
!> =============================================================================
SUBROUTINE drag (self, position, v, v_gas, c)
  class(standard_disk_t):: self
  real, dimension(:,:):: position, v, v_gas
  real, dimension(:):: c
  integer:: i, m
  real, parameter:: density_solid_cgs = 3.
  real, parameter:: particle_size_cgs = 0.1
  real, dimension(:), allocatable:: rm, sigma, temp, cs, h, density, v_kepler, v_g, dv2
  real:: density_solid, particle_size
  !.............................................................................
  m = size(c)
  density_solid = density_solid_cgs/scaling%d
  particle_size = particle_size_cgs/scaling%l
  allocate(rm(m),sigma(m),temp(m),cs(m),h(m),density(m),v_kepler(m),v_g(m),dv2(m))
  do i=1,m
    rm(i) = sqrt(position(i,1)**2+position(i,2)**2+position(i,3)**2)
    sigma(i) = self%sigma*rm(i)**self%sigma_p
    temp(i) = self%temp*rm(i)**self%temp_p
    cs(i) = sqrt(self%gamma*cgs%k_b*temp(i)/(self%mu*cgs%m_u))/scaling%u
    v_kepler(i) = sqrt(star%g*star%m/rm(i))
    h(i) = rm(i)*cs(i)/v_kepler(i)
    v_g(i) = v_kepler(i)*(1.-(cs(i)/v_kepler(i))**2)
    v_gas (i,1) =  v_g(i)*position(i,2)/rm(i)
    v_gas (i,2) = -v_g(i)*position(i,1)/rm(i)
    v_gas (i,3) =  0.0
    dv2 = (v_gas(i,1)-v(i,1))**2+(v_gas(i,2)-v(i,2))**2+(v_gas(i,3)-v(i,3))**2
    density(i) = sigma(i)/h(i)*exp(-(position(i,3)/h(i))**2)
    c(i) = (density(i)/density_solid)*sqrt(cs(i)**2+dv2(i))/particle_size
  end do
END SUBROUTINE drag

END MODULE standard_disk_mod
