!===============================================================================
!> $Id$
!> HD module for the Sedov blast wave exercise
!===============================================================================
MODULE hd_solver_mod
  USE scalar_mod
  USE vector_mod
  USE tensor_mod
  USE patch_mod
  USE operators_mod
  USE eos_mod
  implicit none
  private
  integer, parameter:: i_d=1, i_s=2, i_px=3, i_py=4, i_pz=5
  type, public, extends(patch_t):: hd_t
    type(scalar_t):: d, dddt
    type(scalar_t):: s, dsdt
    type(vector_t):: p, dpdt, rho_face
    integer:: i_d=i_d, i_s=i_s, i_px=i_px, i_py=i_py, i_pz=i_pz
    real:: nu(2)
  contains
    procedure:: init
    procedure:: update
    procedure:: log_density
    procedure:: log_pressure
    procedure:: velocity_magnitude
    procedure:: magnetic_field_magnitude
  end type
  type(operator_t):: o
CONTAINS

!===============================================================================
!> Map variables to the current memory slot and initialize the deriv operators.
!> hd_t%init is called from timestep_t%init, again called from task_list_t%init.
!===============================================================================
SUBROUTINE init (self)
  class(hd_t):: self
  logical:: first_time=.true.
  real:: nu(2)=[0.5,2.0]
  namelist /hd_params/ gamma, nu
  !.............................................................................
  call trace_begin ('hd_mod:hd_t%init', 2)
  call self%patch_t%init
  call o%init (self%gn, self%ds)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, hd_params)
    if (io%master) write (*,hd_params)
  end if
  self%nu = nu
  !-----------------------------------------------------------------------------
  ! Memory mapping; this is cheap, and may be done every time
  !-----------------------------------------------------------------------------
  self%d%f    => self%mem(:,:,:,i_d ,self%it,1)
  self%s%f    => self%mem(:,:,:,i_s ,self%it,1)
  self%p%x    => self%mem(:,:,:,i_px,self%it,1)
  self%p%y    => self%mem(:,:,:,i_py,self%it,1)
  self%p%z    => self%mem(:,:,:,i_pz,self%it,1)
  self%dddt%f => self%mem(:,:,:,i_d ,self%it,2)
  self%dsdt%f => self%mem(:,:,:,i_s ,self%it,2)
  self%dpdt%x => self%mem(:,:,:,i_px,self%it,2)
  self%dpdt%y => self%mem(:,:,:,i_py,self%it,2)
  self%dpdt%z => self%mem(:,:,:,i_pz,self%it,2)
  self%nv = 5
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Evaluates the time derivatives of the hydrodynamic variables
!===============================================================================
SUBROUTINE update (self)
  class(hd_t):: self
  type(scalar_t), save:: viscosity, c_sound, log_rho, artificial_pressure, &
                   diffusion_speed, viscous_dissipation, entropy_per_mass, gas_pressure
  type(vector_t), save:: velocity, delta_u, mass_flux, entropy_flux, &
                   mass_diffusion
  type(tensor_t), save:: reynolds_stress, strain_tensor, viscous_stress, momentum_flux
associate (rho=>self%d, s=>self%s, p=>self%p, dd_dt=>self%dddt, ds_dt=>self%dsdt, &
           ds=>self%ds, dp_dt=>self%dpdt, gas_gamma=>eos%gamma, nu=>self%nu)
  !.............................................................................
  call trace_begin ('hd_mod::hd_t%update')
  call entropy_per_mass%init    (self%gn, self%size/self%n, 'entropy')
  call gas_pressure%init        (self%gn, self%size/self%n, 'p_gas')
  call c_sound%init             (self%gn, self%size/self%n, 'c_sound')
  call delta_u%init             (self%gn, self%size/self%n, 'delta_u')
  call velocity%init            (self%gn, self%size/self%n, 'veloc')
  call log_rho%init             (self%gn, self%size/self%n, 'log_d')
  call self%rho_face%init       (self%gn, self%size/self%n, 'rho_f')
  call artificial_pressure%init (self%gn, self%size/self%n, 'p_artif')
  call viscous_dissipation%init (self%gn, self%size/self%n, 'visc_d')
  call mass_flux%init           (self%gn, self%size/self%n, 'mass_fl')
  call entropy_flux%init        (self%gn, self%size/self%n, 'entr_fl')
  call diffusion_speed%init     (self%gn, self%size/self%n, 'dif_sp')
  call mass_diffusion%init      (self%gn, self%size/self%n, 'mass_dif')
  call reynolds_stress%init     (self%gn, self%size/self%n, 'reyn_str')
  call strain_tensor%init       (self%gn, self%size/self%n, 'strain')
  call viscous_stress%init      (self%gn, self%size/self%n, 'visc')
  call momentum_flux%init       (self%gn, self%size/self%n, 'mom_fl')
  !-----------------------------------------------------------------------------
  ! Memory mapping; this is cheap, and may be done every time
  !-----------------------------------------------------------------------------
  self%d%f    => self%mem(:,:,:,i_d ,self%it,1)
  self%s%f    => self%mem(:,:,:,i_s ,self%it,1)
  self%p%x    => self%mem(:,:,:,i_px,self%it,1)
  self%p%y    => self%mem(:,:,:,i_py,self%it,1)
  self%p%z    => self%mem(:,:,:,i_pz,self%it,1)
  self%dddt%f => self%mem(:,:,:,i_d ,self%it,2)
  self%dsdt%f => self%mem(:,:,:,i_s ,self%it,2)
  self%dpdt%x => self%mem(:,:,:,i_px,self%it,2)
  self%dpdt%y => self%mem(:,:,:,i_py,self%it,2)
  self%dpdt%z => self%mem(:,:,:,i_pz,self%it,2)
  !-----------------------------------------------------------------------------
  ! Equation of state for ideal gas
  !-----------------------------------------------------------------------------
  entropy_per_mass%f        = s%f/rho%f
  if (gas_gamma==1.0) then
    gas_pressure%f          = rho%f
    c_sound%f = 1.0
  else
    if (minval(rho%f) < 0.0) then
      print*,'density negative: id, loc=', self%id, minloc(rho%f)
      call self%output
    end if
    gas_pressure%f = rho%f**gas_gamma*exp(entropy_per_mass%f*(gas_gamma-1.0))
    c_sound%f = sqrt(gas_gamma*(gas_pressure%f/rho%f))
  end if
  !-----------------------------------------------------------------------------
  ! Velocity field from linear momenta
  !-----------------------------------------------------------------------------
  log_rho%f   = log(rho%f)
  self%rho_face%x  = exp(o%xdn(log_rho%f))
  self%rho_face%y  = exp(o%ydn(log_rho%f))
  self%rho_face%z  = exp(o%zdn(log_rho%f))
  velocity%x  = p%x/self%rho_face%x
  velocity%y  = p%y/self%rho_face%y
  velocity%z  = p%z/self%rho_face%z
  self%u_max = self%fmaxval(sqrt(velocity%x**2+velocity%y**2+velocity%z**2)+c_sound%f)
  !-----------------------------------------------------------------------------
  ! Diffusion added for stability
  !-----------------------------------------------------------------------------
  diffusion_speed%f = nu(1)*c_sound%f
  delta_u%x         = (ds(1)*nu(2))*o%ddxup(velocity%x)
  delta_u%y         = (ds(3)*nu(2))*o%ddyup(velocity%y)
  delta_u%z         = (ds(3)*nu(2))*o%ddzup(velocity%z)
  artificial_pressure%f = rho%f*min(0.0,delta_u%x + delta_u%y + delta_u%z)**2
  !-----------------------------------------------------------------------------
  ! Evolution of mass density, with diffusive term proportional to -\grad\ln\rho
  !-----------------------------------------------------------------------------
  mass_diffusion%x  =  - ds(1)*self%rho_face%x*o%xdn(diffusion_speed%f)*o%ddxdn(log_rho%f)
  mass_diffusion%y  =  - ds(2)*self%rho_face%y*o%ydn(diffusion_speed%f)*o%ddydn(log_rho%f)
  mass_diffusion%z  =  - ds(3)*self%rho_face%z*o%zdn(diffusion_speed%f)*o%ddzdn(log_rho%f)
  mass_flux%x        = p%x + mass_diffusion%x
  mass_flux%y        = p%y + mass_diffusion%y
  mass_flux%z        = p%z + mass_diffusion%z
  dd_dt%f = - o%ddxup(mass_flux%x)- o%ddyup(mass_flux%y)- o%ddzup(mass_flux%z)
  self%u_max = max(self%u_max, ds(1)*self%fmaxval(abs(dd_dt%f/rho%f)))
  !-----------------------------------------------------------------------------
  ! Evolution of linear momenta = equation-of-motion
  !-----------------------------------------------------------------------------
  reynolds_stress%xx = o%xup(velocity%x)**2
  reynolds_stress%yy = o%yup(velocity%y)**2
  reynolds_stress%zz = o%zup(velocity%z)**2
  reynolds_stress%xy = o%ydn(velocity%x)*o%xdn(velocity%y)
  reynolds_stress%yz = o%zdn(velocity%y)*o%ydn(velocity%z)
  reynolds_stress%zx = o%xdn(velocity%z)*o%zdn(velocity%x)
  !-----------------------------------------------------------------------------
  ! Strain tensor {s_{ij} = \partial u_i}/{\partial x_j}
  !-----------------------------------------------------------------------------
  strain_tensor%xx = o%ddxup(velocity%x)
  strain_tensor%yy = o%ddyup(velocity%y)
  strain_tensor%zz = o%ddzup(velocity%z)
  strain_tensor%xy = 0.5*(o%ddydn(velocity%x) + o%ddxdn(velocity%y))
  strain_tensor%yz = 0.5*(o%ddzdn(velocity%y) + o%ddydn(velocity%z))
  strain_tensor%zx = 0.5*(o%ddxdn(velocity%z) + o%ddzdn(velocity%x))
  !-----------------------------------------------------------------------------
  ! (negative) viscous stress tensor in per-unit-mass variables
  !-----------------------------------------------------------------------------
  viscous_stress%xx = -ds(1)*diffusion_speed%f*strain_tensor%xx
  viscous_stress%yy = -ds(2)*diffusion_speed%f*strain_tensor%yy
  viscous_stress%zz = -ds(3)*diffusion_speed%f*strain_tensor%zz
  viscous_stress%xy = -0.5*(ds(1)+ds(2))*o%xdn(o%ydn(diffusion_speed%f))*strain_tensor%xy
  viscous_stress%yz = -0.5*(ds(2)+ds(3))*o%ydn(o%zdn(diffusion_speed%f))*strain_tensor%yz
  viscous_stress%zx = -0.5*(ds(3)+ds(1))*o%zdn(o%xdn(diffusion_speed%f))*strain_tensor%zx
  !-----------------------------------------------------------------------------
  ! Viscous dissipation per unit volume
  !-----------------------------------------------------------------------------
  viscous_dissipation%f = rho%f*diffusion_speed%f*( &
    ds(1)*strain_tensor%xx**2 + &
    ds(2)*strain_tensor%yy**2 + &
    ds(3)*strain_tensor%zz**2 + &
   (ds(1)+ds(2))*o%xup(o%yup(strain_tensor%xy))**2 + &
   (ds(2)+ds(3))*o%zup(o%xup(strain_tensor%yz))**2 + &
   (ds(3)+ds(1))*o%zup(o%xup(strain_tensor%zx))**2)
  !-----------------------------------------------------------------------------
  ! add up contributions to momentum flux
  !-----------------------------------------------------------------------------
  momentum_flux%xx = rho%f*(reynolds_stress%xx + viscous_stress%xx) + gas_pressure%f + artificial_pressure%f
  momentum_flux%yy = rho%f*(reynolds_stress%yy + viscous_stress%yy) + gas_pressure%f + artificial_pressure%f
  momentum_flux%zz = rho%f*(reynolds_stress%zz + viscous_stress%zz) + gas_pressure%f + artificial_pressure%f
  momentum_flux%xy = exp(o%xdn(o%ydn(log_rho%f)))*(reynolds_stress%xy + viscous_stress%xy)
  momentum_flux%yz = exp(o%ydn(o%zdn(log_rho%f)))*(reynolds_stress%yz + viscous_stress%yz)
  momentum_flux%zx = exp(o%zdn(o%xdn(log_rho%f)))*(reynolds_stress%zx + viscous_stress%zx)
  !-----------------------------------------------------------------------------
  ! Add diffusion of momentum consistent with diffusion of mass
  !-----------------------------------------------------------------------------
  momentum_flux%xx = momentum_flux%xx + o%xup(mass_diffusion%x)*o%xup(velocity%x)
  momentum_flux%yy = momentum_flux%yy + o%yup(mass_diffusion%y)*o%yup(velocity%y)
  momentum_flux%zz = momentum_flux%zz + o%zup(mass_diffusion%z)*o%zup(velocity%z)
  momentum_flux%xy = momentum_flux%xy + 0.5*(o%ydn(mass_diffusion%x*velocity%x)+ o%xdn(mass_diffusion%y*velocity%y))
  momentum_flux%yz = momentum_flux%yz + 0.5*(o%zdn(mass_diffusion%y*velocity%y)+ o%ydn(mass_diffusion%z*velocity%z))
  momentum_flux%zx = momentum_flux%zx + 0.5*(o%xdn(mass_diffusion%z*velocity%z)+ o%zdn(mass_diffusion%x*velocity%x))
  !-----------------------------------------------------------------------------
  ! dp_dt = -div(momentum_flux)
  !-----------------------------------------------------------------------------
  dp_dt%x = - o%ddxdn(momentum_flux%xx) - o%ddyup(momentum_flux%xy) - o%ddzup(momentum_flux%zx)
  dp_dt%y = - o%ddxup(momentum_flux%xy) - o%ddydn(momentum_flux%yy) - o%ddzup(momentum_flux%yz)
  dp_dt%z = - o%ddxup(momentum_flux%zx) - o%ddyup(momentum_flux%yz) - o%ddzdn(momentum_flux%zz)
  !-----------------------------------------------------------------------------
  ! Evolution of entropy:  {\partial s}/{\partial t} = -div(entropy_flux) + Q/T
  !-----------------------------------------------------------------------------
  entropy_flux%x  =  - ds(1)*o%xdn(diffusion_speed%f) * o%ddxdn (entropy_per_mass%f)
  entropy_flux%y  =  - ds(2)*o%ydn(diffusion_speed%f) * o%ddydn (entropy_per_mass%f)
  entropy_flux%z  =  - ds(3)*o%zdn(diffusion_speed%f) * o%ddzdn (entropy_per_mass%f)
  !-----------------------------------------------------------------------------
  ! Add diffusion of entropy consistent with diffusion of mass
  !-----------------------------------------------------------------------------
  entropy_flux%x  =  entropy_flux%x + mass_diffusion%x*o%xdn(entropy_per_mass%f)
  entropy_flux%y  =  entropy_flux%y + mass_diffusion%y*o%ydn(entropy_per_mass%f)
  entropy_flux%z  =  entropy_flux%z + mass_diffusion%z*o%zdn(entropy_per_mass%f)
  ds_dt%f = - o%ddxup(entropy_flux%x) - o%ddyup(entropy_flux%y) - o%ddzup(entropy_flux%z)
!  ds_dt%f = ds_dt%f + viscous_dissipation%f * rho%f/gas_pressure%f
  !-----------------------------------------------------------------------------
    call trace_end
end associate
END SUBROUTINE update

!===============================================================================
!> Compute the log density and other similar things -- for use in refinement
!===============================================================================
SUBROUTINE log_density (self, out)
  class(hd_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: out
  !...............................................................................
  call trace_begin ('hd_mod:hd_t%log_density')
  if (io%verbose>3) then
    print'(1x,a,4i4,3(2x,3i4),3x,a,1p,4e10.2)', &
      'log_d', self%id, self%it, self%new, i_d, self%li, self%ui, self%n, &
      'minmax:', minval(self%d%f), maxval(self%d%f), &
      minval(self%mem(self%li(1):self%ui(1), &
                      self%li(2):self%ui(2), &
                      self%li(3):self%ui(3), i_d, self%it, 1)), &
      maxval(self%mem(self%li(1):self%ui(1), &
                      self%li(2):self%ui(2), &
                      self%li(3):self%ui(3), i_d, self%it, 1))
  end if
  out(1:self%n(1), &
      1:self%n(2), &
      1:self%n(3)) = log(self%mem(self%li(1):self%ui(1), &
                                            self%li(2):self%ui(2), &
                                            self%li(3):self%ui(3), i_d, self%it, 1))
  call trace_end
END SUBROUTINE log_density
!===============================================================================
SUBROUTINE log_pressure (self, out)
  class(hd_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: out
  !...............................................................................
  out = 0.0
END SUBROUTINE log_pressure
!===============================================================================
SUBROUTINE velocity_magnitude (self, out)
  class(hd_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: out
  !...............................................................................
  out = 0.0
END SUBROUTINE velocity_magnitude
!===============================================================================
SUBROUTINE magnetic_field_magnitude (self, out)
  class(hd_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: out
  !...............................................................................
  out = 0.0
END SUBROUTINE magnetic_field_magnitude
!===============================================================================
END MODULE hd_solver_mod
