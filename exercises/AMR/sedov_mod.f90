!> =============================================================================
!> $Id$
!> HD module for the Sedov blast wave exercise
!> =============================================================================
MODULE pde_mod
  USE scalar_mod
  USE vector_mod
  USE tensor_mod
  USE patches_mod
  USE operators_mod
  implicit none
  private
  integer, parameter:: i_d=1, i_s=2, i_px=3, i_py=4, i_pz=5
  type, public, extends(patch_t):: pde_t
    type(scalar_t):: d, dddt
    type(scalar_t):: s, dsdt
    type(vector_t):: p, dpdt
    real:: nu(2)
    real:: gamma
    real:: amplitude, amplitude_rho
    real:: w_gauss
    real:: p_gauss
    real:: u_max
  contains
    procedure:: init
    procedure:: update
    procedure:: log_density
    procedure:: log_pressure
    procedure:: velocity_magnitude
    procedure:: magnetic_field_magnitude
  end type
  type(operator_t), pointer:: o
  public experiment_name
CONTAINS

!> =============================================================================
FUNCTION experiment_name()
  character(len=64):: experiment_name
  experimentname = 'sedov'
  experiment_name = experimentname
END FUNCTION
!> =============================================================================
!> Map variables to the current memory slot and initialize the deriv operators.
!> pde_t%init is called from timestep_t%init, again called from task_list_t%init.
!> =============================================================================
SUBROUTINE init(self, position, size)
  class(pde_t):: self
  real:: position(3), size(3)
  type(tensor_t):: tensor
  integer:: i, j
  real:: x, y, f
  logical:: first_time=.true.
  real:: w_gauss=3.0, p_gauss=4.0, nu(2)=[0.5,2.0], &
         amplitude=3.0, amplitude_rho=3., gamma=5./3.
  namelist /sedov_params/ w_gauss, p_gauss, amplitude, amplitude_rho, gamma, nu
  !.............................................................................
  call trace_begin ('hd_mod:pde_t%init', 2)
  if (first_time) then
    first_time=.false.
    rewind (io%input)
    read (io%input, sedov_params)
    if (io%master) write (*,sedov_params)
  end if
  self%w_gauss = w_gauss
  self%p_gauss = p_gauss
  self%amplitude = amplitude
  self%amplitude_rho = amplitude_rho
  self%gamma = gamma
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
  if (.not.associated(o)) allocate (o)
  print*,self%gn,self%size,self%n
  call o%init (self%gn, self%size/self%n)
  !-----------------------------------------------------------------------------
  ! Physical initial conditions: Gaussian top hat in entropy.  Only at time=0.0
  !-----------------------------------------------------------------------------
  self%nv = 5
  if (self%time==0.0) then
    self%d%f = 1.0
    self%s%f = 0.0
    self%p%x = 0.0
    self%p%y = 0.0
    self%p%z = 0.0
    do j=1,self%gn(2)
    do i=1,self%gn(1)
      x = (i-(self%gn(1)+1)*0.5) + self%position(1)/self%ds(1)
      y = (j-(self%gn(2)+1)*0.5) + self%position(2)/self%ds(2)
      f = exp(-(sqrt(x**2+y**2)/self%w_gauss)**self%p_gauss)
      self%d%f(i,j,:) = exp(self%amplitude_rho*f)
      self%s%f(i,j,:) = self%amplitude*f*self%d%f(i,j,:)
    end do
    end do
    if (io%verbose>1) then
      print'(i4,2x,3i4,1p,6e10.2)', &
        self%id,self%n,self%w_gauss,self%p_gauss,self%position
    else
      print'(1x,a,i4,2(3f8.3,2x),a,1p,2e10.2)', &
        'initial:', self%id, self%position, self%size, &
        '  min, max:', &
        minval(self%mem(:,:,:,1,1,1)), &
        maxval(self%mem(:,:,:,1,1,1))
    end if
  end if
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Evaluates the time derivatives of the hydrodynamic variables
!> =============================================================================
SUBROUTINE update (self)
  class(pde_t):: self
  type(scalar_t), save:: viscosity, c_sound, log_rho, artificial_pressure, &
                   diffusion_speed, viscous_dissipation, entropy_per_mass, gas_pressure
  type(vector_t), save:: velocity, delta_u, mass_flux, entropy_flux, &
                   mass_diffusion, rho_face
  type(tensor_t), save:: reynolds_stress, strain_tensor, viscous_stress, momentum_flux
associate (rho=>self%d, s=>self%s, p=>self%p, dd_dt=>self%dddt, ds_dt=>self%dsdt, &
           ds=>self%ds, dp_dt=>self%dpdt, gas_gamma=>self%gamma, nu=>self%nu)
  !.............................................................................
  call trace_begin ('sedov_mod::pde_t%update')
  call entropy_per_mass%init    (self%gn, self%size/self%n, 'entropy')
  call gas_pressure%init        (self%gn, self%size/self%n, 'p_gas')
  call c_sound%init             (self%gn, self%size/self%n, 'c_sound')
  call delta_u%init             (self%gn, self%size/self%n, 'delta_u')
  call velocity%init            (self%gn, self%size/self%n, 'veloc')
  call log_rho%init             (self%gn, self%size/self%n, 'log_d')
  call rho_face%init            (self%gn, self%size/self%n, 'rho_f')
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
  rho_face%x  = exp(o%xdn(log_rho%f))
  rho_face%y  = exp(o%ydn(log_rho%f))
  rho_face%z  = exp(o%zdn(log_rho%f))
  velocity%x  = p%x/rho_face%x
  velocity%y  = p%y/rho_face%y
  velocity%z  = p%y/rho_face%z
  self%u_max = maxval(sqrt(velocity%x**2+velocity%y**2+velocity%z**2)+c_sound%f)
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
  mass_diffusion%x  =  - ds(1)*rho_face%x*o%xdn(diffusion_speed%f)*o%ddxdn(log_rho%f)
  mass_diffusion%y  =  - ds(2)*rho_face%y*o%ydn(diffusion_speed%f)*o%ddydn(log_rho%f)
  mass_diffusion%z  =  - ds(3)*rho_face%z*o%zdn(diffusion_speed%f)*o%ddzdn(log_rho%f)
  mass_flux%x        = p%x + mass_diffusion%x
  mass_flux%y        = p%y + mass_diffusion%y
  mass_flux%z        = p%z + mass_diffusion%z
  dd_dt%f = - o%ddxup(mass_flux%x)- o%ddyup(mass_flux%y)- o%ddzup(mass_flux%z)
  self%u_max = max(self%u_max, ds(1)/maxval(dd_dt%f/rho%f))
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

!> =============================================================================
!> Compute the log density and other similar things -- for use in refinement
!> =============================================================================
SUBROUTINE log_density (self, out)
  class(pde_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: out
  !...............................................................................
  call trace_begin ('hd_mod:pde_t%log_density')
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
  out(self%li(1):self%ui(1), &
      self%li(2):self%ui(2), &
      self%li(3):self%ui(3)) = log(self%mem(self%li(1):self%ui(1), &
                                            self%li(2):self%ui(2), &
                                            self%li(3):self%ui(3), i_d, self%it, 1))
  call trace_end
END SUBROUTINE log_density
!> =============================================================================
SUBROUTINE log_pressure (self, out)
  class(pde_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: out
  !...............................................................................
  out = 0.0
END SUBROUTINE log_pressure
!> =============================================================================
SUBROUTINE velocity_magnitude (self, out)
  class(pde_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: out
  !...............................................................................
  out = 0.0
END SUBROUTINE velocity_magnitude
!> =============================================================================
SUBROUTINE magnetic_field_magnitude (self, out)
  class(pde_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: out
  !...............................................................................
  out = 0.0
END SUBROUTINE magnetic_field_magnitude
!> =============================================================================
END MODULE pde_mod
