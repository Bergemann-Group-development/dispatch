MODULE disk_model_mod
  USE io_mod
  USE trace_mod
  USE scaling_mod
  USE units_mod
  USE math_mod
  USE eos_mod
  implicit none
  private
  type, public:: disk_model_t
    real :: m_star      = 0.0
    real :: sigma_ratio = 0.0
  contains
    procedure:: init
    procedure:: r_au
    procedure:: sigma
    procedure:: sigma_grad
    procedure:: sigma_dust
    procedure:: temperature
    procedure:: temperature_grad
    procedure:: c_sound
    procedure:: c_sound_grad
    procedure:: c_sound_2_grad
    procedure:: d_midplane
    procedure:: d_midplane_grad
    procedure:: d_midplane_dust
    procedure:: h_scale
    procedure:: h_scale_grad
    procedure:: h_scale_dust
    procedure:: v_kepler
    procedure:: v_kepler_grad
    procedure:: v_head
  end type
  !-----------------------------------------------------------------------------
  ! Shared module parameters
  !-----------------------------------------------------------------------------
  real:: m_solar     = 1.0    ! stellar mass in Solar masses
  real:: temp_1      = 280.   ! temperature at 1 au (in K)
  real:: temp_exp    = -0.5   ! exponent of power-law temperature profile
  real:: sigma_1     = 1700.  ! surface density at 1 au (in g/cm^2)
  real:: sigma_exp   = -1.0   ! exponent of power-law surface density profile
  real:: sigma_ratio = 0.05   ! dust-to-gas surface density ratio
  real:: h_ratio     = 0.1    ! initial dust-to-gas scale height ratio
  type(disk_model_t), public:: disk_model
CONTAINS

!===============================================================================
!> Initialize parameters
!===============================================================================
SUBROUTINE init (self)
  class(disk_model_t):: self
  integer            :: iostat
  logical, save      :: first_time = .true.
  real(8)            :: au
  namelist /disk_model_params/ &
    m_solar, sigma_1, sigma_exp, temp_1, temp_exp, sigma_ratio, h_ratio
  !-----------------------------------------------------------------------------
  call io%lock%set ('disk_model_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, disk_model_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, disk_model_params)
    end if
    call io%namelist_warning ('disk_model_params', iostat)
    write (io%namelist, disk_model_params)
    !...........................................................................
    call io%header ('disk_model:')
    self%m_star = m_solar*cgs%m_sun/scaling%m
    self%sigma_ratio = sigma_ratio
    au = cgs%AU/scaling%l
    write (stdout,1) '   sigma:', self%sigma (au)
    write (stdout,1) '   d_mid:', self%d_midplane (au)
    write (stdout,1) '    temp:', self%temperature (au)
    write (stdout,1) ' c_sound:', self%c_sound (au)
    write (stdout,1) 'v_kepler:', self%v_kepler (au)
    write (stdout,1) '  v_head:', self%v_head (au, 0d0)
    write (stdout,1) ' h_scale:', self%h_scale (au)
    1 format(2x,a,1p,e12.3)
    first_time = .false.
  end if
  call io%lock%unset ('disk_model_t%init')
END SUBROUTINE init

!===============================================================================
!> Radius in au
!===============================================================================
FUNCTION r_au (self, r) RESULT (rr)
  class(disk_model_t):: self
  real(8):: r, rr
  !.............................................................................
  rr = r*scaling%l/cgs%AU
END FUNCTION r_au

!===============================================================================
!> Surface density
!===============================================================================
FUNCTION sigma (self, r) RESULT (s)
  class(disk_model_t):: self
  real(8):: r, s
  !.............................................................................
  s = sigma_1*self%r_au(r)**sigma_exp/(scaling%d*scaling%l)
END FUNCTION sigma

!===============================================================================
!> Radial gradient of surface density
!===============================================================================
FUNCTION sigma_grad (self, r) RESULT (s_grad)
  class(disk_model_t):: self
  real(8):: r, s_grad
  !.............................................................................
  s_grad = sigma_exp*sigma_1*self%r_au(r)**(sigma_exp-1d0)/(cgs%AU/scaling%d)
END FUNCTION sigma_grad

!===============================================================================
!> Dust surface density
!===============================================================================
FUNCTION sigma_dust (self, r) RESULT (s)
  class(disk_model_t):: self
  real(8):: r, s
  !.............................................................................
  s = sigma_ratio*sigma_1*self%r_au(r)**sigma_exp/(scaling%d*scaling%l)
END FUNCTION sigma_dust

!===============================================================================
!> Vertical scale height
!===============================================================================
FUNCTION h_scale (self, r) RESULT (h)
  class(disk_model_t):: self
  real(8):: r, h
  !.............................................................................
  h = r*self%c_sound(r)/self%v_kepler(r)/sqrt(eos%gamma)
END FUNCTION h_scale

!===============================================================================
!> Radial gradient of vertical scale height
!===============================================================================
FUNCTION h_scale_grad (self, r) RESULT (h_grad)
  class(disk_model_t):: self
  real(8):: r, v, h_grad
  !.............................................................................
  v = self%v_kepler(r)
  h_grad = (r*v*self%c_sound_grad(r) + &
    self%c_sound(r)*(v-r*self%v_kepler_grad(r))) / &
    (v**2*sqrt(eos%gamma))
END FUNCTION h_scale_grad

!===============================================================================
!> Vertical dust scale height
!===============================================================================
FUNCTION h_scale_dust (self, r) RESULT (h)
  class(disk_model_t):: self
  real(8):: r, h
  !.............................................................................
  h = h_ratio*r*self%c_sound(r)/self%v_kepler(r)/sqrt(eos%gamma)
END FUNCTION h_scale_dust

!===============================================================================
!> Temperature
!===============================================================================
FUNCTION temperature (self, r) RESULT (tt)
  class(disk_model_t):: self
  real(8):: r, tt
  !.............................................................................
  tt = temp_1*self%r_au(r)**temp_exp/scaling%temp
END FUNCTION temperature

!===============================================================================
!> Radial gradient of temperature
!===============================================================================
FUNCTION temperature_grad (self, r) RESULT (tt_grad)
  class(disk_model_t):: self
  real(8):: r, tt_grad
  !.............................................................................
  tt_grad = temp_exp*temp_1*self%r_au(r)**(temp_exp-1d0) / &
    (cgs%AU/scaling%l*scaling%temp)
END FUNCTION temperature_grad

!===============================================================================
!> Mid-plane density
!===============================================================================
FUNCTION d_midplane (self, r) RESULT (d)
  class(disk_model_t):: self
  real(8):: r, d
  !.............................................................................
  d = self%sigma(r)/(self%h_scale(r)*sqrt(2.*math%pi))
END FUNCTION d_midplane

!===============================================================================
!> Radial gradient of mid-plane density
!===============================================================================
FUNCTION d_midplane_grad (self, r) RESULT (d_grad)
  class(disk_model_t):: self
  real(8):: r, h, d_grad
  !.............................................................................
  h = self%h_scale(r)
  d_grad = (h*self%sigma_grad(r) - self%sigma(r)*self%h_scale_grad(r)) / &
    (h**2*sqrt(2.*math%pi))
END FUNCTION d_midplane_grad

!===============================================================================
!> Mid-plane dust density
!===============================================================================
FUNCTION d_midplane_dust (self, r) RESULT (d)
  class(disk_model_t):: self
  real(8):: r, d
  !.............................................................................
  d = self%sigma_dust(r)/(self%h_scale_dust(r)*sqrt(2.*math%pi))
END FUNCTION d_midplane_dust

!===============================================================================
!> Sound speed
!===============================================================================
FUNCTION c_sound(self, r) RESULT (c)
  class(disk_model_t):: self
  real(8):: r, c
  !.............................................................................
  c = sqrt(eos%gamma*self%temperature(r)/eos%mu)
END FUNCTION c_sound

!===============================================================================
!> Radial gradient of sound speed
!===============================================================================
FUNCTION c_sound_grad(self, r) RESULT (c_grad)
  class(disk_model_t):: self
  real(8):: r, c_grad
  !.............................................................................
  c_grad = self%temperature_grad(r)*self%c_sound(r)/2d0/self%temperature(r)
END FUNCTION c_sound_grad

!===============================================================================
!> Radial gradient of square of sound speed
!===============================================================================
FUNCTION c_sound_2_grad(self, r) RESULT (c_2_grad)
  class(disk_model_t):: self
  real(8):: r, c_2_grad
  !.............................................................................
  c_2_grad = self%temperature_grad(r)*self%c_sound(r)**2/self%temperature(r)
END FUNCTION c_sound_2_grad

!===============================================================================
!> Kepler velocity
!===============================================================================
FUNCTION v_kepler (self, r) RESULT (v)
  class(disk_model_t):: self
  real(8):: r, v
  !.............................................................................
  v = sqrt(scaling%grav*self%m_star/max(r,0.1d0))
END FUNCTION v_kepler

!===============================================================================
!> Radial gradient of Kepler velocity
!===============================================================================
FUNCTION v_kepler_grad (self, r) RESULT (v_grad)
  class(disk_model_t):: self
  real(8):: r, v_grad
  !.............................................................................
  v_grad = -self%v_kepler(r)/2d0/r
END FUNCTION v_kepler_grad

!===============================================================================
!> Head wind velocity
!===============================================================================
FUNCTION v_head (self, r, z) RESULT (v)
  class(disk_model_t):: self
  real(8), intent(in) :: r, z
  real(8)             :: h, dmid, d, exp_grad, d_grad, p_grad, v_k, v
  !.............................................................................
  h = self%h_scale (r)                                    ! scale height
  dmid = self%d_midplane (r)                              ! mid-plane density
  d = dmid*exp(-(z/h)**2/2d0)                             ! Gaussian profile
  exp_grad = exp(-(z/h)**2/2d0) * &                       ! Gaussian profile
    z**2*self%h_scale_grad(r)/h**3                        ! gradient
  d_grad = exp(-(z/h)**2/2d0)*self%d_midplane_grad(r) + & ! density gradient
    dmid*exp_grad
  p_grad = d*self%c_sound_2_grad(r) + &                   ! pressure gradient
    self%c_sound(r)**2*d_grad
  v_k = self%v_kepler (r)                                 ! Keplerian velocity
  v = sqrt(v_k**2+r*p_grad/d)-v_k                         ! headwind speed
END FUNCTION v_head

END MODULE disk_model_mod