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
    real:: gamma = 5./3.
    real:: amplitude
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
  experiment_name = 'sedov'
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
  real:: x, y
  real:: w_gauss=3.0, p_gauss=4.0, amplitude=2.0, nu(2)=[0.5,2.0]
  namelist /sedov_params/ w_gauss, p_gauss, amplitude, nu
  !.............................................................................
  call trace_begin ('hd_mod:pde_t%init', 2)
  rewind (io%input)
  read (io%input, sedov_params)
  if (io%master) write (*,sedov_params)
  self%w_gauss = w_gauss
  self%p_gauss = p_gauss
  self%amplitude = amplitude
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
  allocate (o)
  print*,self%gn,self%size,self%n
  call o%init (self%gn, self%size/self%n)
  !-----------------------------------------------------------------------------
  ! Physical initial conditions: Gaussian top hat in entropy.  Only at time=0.0
  !-----------------------------------------------------------------------------
  self%nv = 5
  if (self%time==0.0) then
    self%d%f = 1.0
    print*,'MK1', self%it, self%new
    self%p%x = 0.0
    self%p%y = 0.0
    self%p%z = 0.0
    do j=1,self%gn(2)
    do i=1,self%gn(1)
      x = (i-(self%gn(1)+1)*0.5) + self%position(1)/self%ds(1)
      y = (j-(self%gn(2)+1)*0.5) + self%position(2)/self%ds(2)
      self%s%f(i,j,:) = self%amplitude*exp(-(x**2+y**2)/self%w_gauss**2)**self%p_gauss
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
  type(scalar_t):: viscosity, c_sound, log_density, artificial_pressure, &
                   viscous_dissipation, entropy_per_mass, gas_pressure
  type(vector_t):: velocity, delta_u, mass_flux, entropy_flux, &
                   diffusion_speed, mass_diffusion
  type(tensor_t):: reynolds_stress, strain_tensor, viscous_stress, momentum_flux
associate (rho=>self%d, s=>self%s, p=>self%p, dd_dt=>self%dddt, ds_dt=>self%dsdt, &
           ds=>self%ds, dp_dt=>self%dpdt, gas_gamma=>self%gamma, nu=>self%nu)
  !.............................................................................
  call trace_begin ('sedov_mod::pde_t%update')
  call entropy_per_mass%init    (self%gn, self%size/self%n, 'entropy')
  call gas_pressure%init        (self%gn, self%size/self%n, 'p_gas')
  call c_sound%init             (self%gn, self%size/self%n, 'c_sound')
  call delta_u%init             (self%gn, self%size/self%n, 'delta_u')
  call velocity%init            (self%gn, self%size/self%n, 'veloc')
  call log_density%init         (self%gn, self%size/self%n, 'log_d')
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
  ! Equation of state for ideal gas
  !-----------------------------------------------------------------------------
  entropy_per_mass%f        = s%f/rho%f
  if (gas_gamma==1.0) then
    gas_pressure%f          = rho%f
    c_sound%f = 1.0
  else
    gas_pressure%f = rho%f**gas_gamma*exp(entropy_per_mass%f*(gas_gamma-1.0))
    c_sound%f = gas_gamma*(gas_pressure%f/rho%f)
    c_sound%f = sqrt(c_sound%f)
  end if
  !-----------------------------------------------------------------------------
  ! Velocity field from linear momenta
  !-----------------------------------------------------------------------------
  log_density%f     = log(rho%f)
  velocity%x        = p%x/exp(o%xdn(log_density%f))
  velocity%y        = p%y/exp(o%ydn(log_density%f))
  velocity%z        = p%z/exp(o%zdn(log_density%f))
  self%u_max = sqrt(maxval(velocity%x**2+velocity%y**2+velocity%z**2)) &
             + maxval(c_sound%f)
  !-----------------------------------------------------------------------------
  ! Diffusion added for stability
  !-----------------------------------------------------------------------------
  diffusion_speed%x = (ds(1)*nu(1))*c_sound%f
  diffusion_speed%y = (ds(2)*nu(1))*c_sound%f
  diffusion_speed%z = (ds(3)*nu(1))*c_sound%f
  delta_u%x         = (ds(1)*nu(2))*o%ddxup(velocity%x)
  delta_u%y         = (ds(1)*nu(2))*o%ddyup(velocity%y)
  delta_u%z         = (ds(1)*nu(2))*o%ddzup(velocity%z)
  artificial_pressure%f = rho%f*min(0.0,delta_u%x + delta_u%y + delta_u%z)**2
  !-----------------------------------------------------------------------------
  ! Evolution of mass density
  !-----------------------------------------------------------------------------
  mass_diffusion%x  =  - o%xdn(diffusion_speed%x) * o%ddxdn (log_density%f)
  mass_diffusion%y  =  - o%ydn(diffusion_speed%y) * o%ddydn (log_density%f)
  mass_diffusion%z  =  - o%zdn(diffusion_speed%z) * o%ddzdn (log_density%f)
  mass_flux%x        = p%x + mass_diffusion%x
  mass_flux%y        = p%y + mass_diffusion%y
  mass_flux%z        = p%z + mass_diffusion%z
  call minus_div_vector (dd_dt, mass_flux)
  !-----------------------------------------------------------------------------
  ! Evolution of linear momenta = equation-of-motion
  !-----------------------------------------------------------------------------
  reynolds_stress%xx = o%xup(velocity%x)*o%xup(p%x)
  reynolds_stress%yy = o%yup(velocity%y)*o%yup(p%y)
  reynolds_stress%zz = o%zup(velocity%z)*o%zup(p%z)
  reynolds_stress%xy = o%ydn(velocity%x)*o%xdn(p%y)
  reynolds_stress%yz = o%zdn(velocity%y)*o%ydn(p%z)
  reynolds_stress%zx = o%xdn(velocity%z)*o%zdn(p%x)

  strain_tensor%xx = o%ddxup(velocity%x)
  strain_tensor%yy = o%ddyup(velocity%y)
  strain_tensor%zz = o%ddzup(velocity%z)
  strain_tensor%xy = 0.5*(o%ddydn(velocity%x) + o%ddxdn(velocity%y))
  strain_tensor%yz = 0.5*(o%ddzdn(velocity%y) + o%ddydn(velocity%z))
  strain_tensor%zx = 0.5*(o%ddxdn(velocity%z) + o%ddzdn(velocity%x))

  viscous_stress%xx = log_density%f*o%xup(diffusion_speed%x)*strain_tensor%xx
  viscous_stress%yy = log_density%f*o%yup(diffusion_speed%y)*strain_tensor%yy
  viscous_stress%zz = log_density%f*o%zup(diffusion_speed%z)*strain_tensor%zz
  viscous_stress%xy = o%ydn(o%xdn(log_density%f)*diffusion_speed%x)*strain_tensor%xy
  viscous_stress%yz = o%zdn(o%ydn(log_density%f)*diffusion_speed%y)*strain_tensor%yz
  viscous_stress%zx = o%xdn(o%zdn(log_density%f)*diffusion_speed%z)*strain_tensor%zx

  viscous_stress%xx = viscous_stress%xx + o%xup(mass_diffusion%x*velocity%x)
  viscous_stress%yy = viscous_stress%yy + o%yup(mass_diffusion%y*velocity%y)
  viscous_stress%zz = viscous_stress%zz + o%zup(mass_diffusion%z*velocity%z)
  viscous_stress%xy = viscous_stress%xy + o%xdn(mass_diffusion%y)*o%ydn(velocity%x)
  viscous_stress%yz = viscous_stress%yz + o%ydn(mass_diffusion%z)*o%zdn(velocity%y)
  viscous_stress%zx = viscous_stress%zx + o%zdn(mass_diffusion%x)*o%xdn(velocity%z)

  viscous_stress%xx = viscous_stress%xx + artificial_pressure%f
  viscous_stress%yy = viscous_stress%yy + artificial_pressure%f
  viscous_stress%zz = viscous_stress%zz + artificial_pressure%f

  viscous_dissipation%f = &
    strain_tensor%xx*viscous_stress%xx + &
    strain_tensor%yy*viscous_stress%yy + &
    strain_tensor%zz*viscous_stress%zz + &
    2*(o%xup(o%yup(strain_tensor%xy*viscous_stress%xy)) &
     + o%yup(o%zup(strain_tensor%yz*viscous_stress%yz)) &
     + o%zup(o%xup(strain_tensor%zx*viscous_stress%zx)))

  viscous_stress%xx = viscous_stress%xx + gas_pressure%f
  viscous_stress%yy = viscous_stress%yy + gas_pressure%f
  viscous_stress%zz = viscous_stress%zz + gas_pressure%f

  momentum_flux%xx = reynolds_stress%xx + viscous_stress%xx
  momentum_flux%yy = reynolds_stress%yy + viscous_stress%yy
  momentum_flux%zz = reynolds_stress%zz + viscous_stress%zz
  momentum_flux%xy = reynolds_stress%xy + viscous_stress%xy
  momentum_flux%yz = reynolds_stress%yz + viscous_stress%yz
  momentum_flux%zx = reynolds_stress%zx + viscous_stress%zx

  call minus_div_tensor (dp_dt, momentum_flux)
  !-----------------------------------------------------------------------------
  ! Evolution of entropy
  !-----------------------------------------------------------------------------
  entropy_flux%x  =  - o%xdn(diffusion_speed%x) * o%ddxdn (entropy_per_mass%f)
  entropy_flux%y  =  - o%ydn(diffusion_speed%y) * o%ddydn (entropy_per_mass%f)
  entropy_flux%z  =  - o%zdn(diffusion_speed%z) * o%ddzdn (entropy_per_mass%f)
  call minus_div_vector (ds_dt, entropy_flux)
  ds_dt%f = ds_dt%f + viscous_dissipation%f * gas_pressure%f/rho%f
  !-----------------------------------------------------------------------------
    call trace_end
end associate
END SUBROUTINE

!> =============================================================================
!> This should produce down-shifted gradients, for face centered fluxes
!> =============================================================================
FUNCTION minus_grad (in1, in2) RESULT (out)
  type(scalar_t):: in2
  type(vector_t):: in1, out
  out%x = - o%xdn(in1%x)*o%ddxdn(in2%f)
  out%y = - o%ydn(in1%y)*o%ddydn(in2%f)
  out%z = - o%zdn(in1%z)*o%ddzdn(in2%f)
END FUNCTION
!> =============================================================================
!> The in argument is downshfted, so the div should have up-derivatives
!> =============================================================================
SUBROUTINE minus_div_vector (out,in)
  type(scalar_t):: out
  type(vector_t):: in
  out%f = - o%ddxup(in%x) - o%ddyup(in%y) - o%ddzup(in%z)
END SUBROUTINE
!> =============================================================================
!> in1%yx is the transport of y-mom by u_x, and is taken ddxup() of for non-
!> diagonal elements, and ddxdn() of for diagonal ones.
!> in2 is velocity, in1 is momentum.   The tensor is symmetric, so only
!> 6 elements are significant.  The off-diagonal elements have o%ddydn(out%xy).
!> FIXME: In IDL the density factor is interpolated separately, which is better
!> =============================================================================
SUBROUTINE outer_product (out, in1, in2)
  type(tensor_t):: out
  type(vector_t):: in1, in2
  print*,'zz'
END SUBROUTINE
!> =============================================================================
!> These elements should be cell centered on the diagonal, and edge centered
!> elsewhere
!> =============================================================================
SUBROUTINE symmetrized_derivatives (out, in1)
  type(tensor_t):: out
  type(vector_t):: in1
END SUBROUTINE
!> =============================================================================
!> The in1 argument is cell centered, the in2 argument is face-centered, and
!> the in3 argument is edge-centered, as should the output be, except in the
!> diagonal, which is cell centered
!> =============================================================================
SUBROUTINE tensor_factor (out, in1, in2, in3)
  type(tensor_t):: out
  type(scalar_t):: in1
  type(vector_t):: in2
  type(tensor_t):: in3
END SUBROUTINE
!> =============================================================================
!> This is to compute the cell-centered sum of squares for dissipation
!> =============================================================================
FUNCTION scalar_tensor_product (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(tensor_t):: in1, in2
END FUNCTION
!> =============================================================================
!> This adds a centered component to the diagonal
!> =============================================================================
SUBROUTINE tensor_add_diagonal (out, in1)
  type(tensor_t):: out
  type(scalar_t):: in1
  out%xx = out%xx + in1%f
  out%yy = out%yy + in1%f
  out%zz = out%zz + in1%f
END SUBROUTINE
!> =============================================================================
!> The first argument is a face centered flux, the 2nd argument is a face-
!> centered velocity, and the output should have a cell centered diagonal and
!> edge-centered off-diagonals
!> =============================================================================
FUNCTION edge_multiply (in1, in2) RESULT (out)
  type(tensor_t):: out
  type(vector_t):: in1
  type(vector_t):: in2
  out%xx = o%xup(in1%x*in2%x)
  out%yy = o%yup(in1%y*in2%y)
  out%zz = o%zup(in1%z*in2%z)
  out%xy = o%xdn(in1%y)*o%ydn(in2%x)
  out%yz = o%ydn(in1%z)*o%zdn(in2%y)
  out%zx = o%zdn(in1%x)*o%xdn(in2%z)
  out%yx = out%xy
  out%zy = out%yz
  out%xz = out%zx
END FUNCTION
!> =============================================================================
SUBROUTINE minus_div_tensor (out,in)
  type(vector_t):: out
  type(tensor_t):: in
  out%x = - o%ddxup(in%xx) - o%ddyup(in%xy) - o%ddzup(in%xz)
  out%y = - o%ddxup(in%yx) - o%ddyup(in%yy) - o%ddzup(in%yz)
  out%z = - o%ddxup(in%zx) - o%ddyup(in%zy) - o%ddzup(in%zz)
END SUBROUTINE

!> =============================================================================
!> Compute the log density and other similar things -- for use in refinement
!> =============================================================================
SUBROUTINE log_density (self, out)
  class(pde_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: out
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
  out = log(self%mem(self%li(1):self%ui(1), &
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
