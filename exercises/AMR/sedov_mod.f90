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
    real:: amplitude = 2.3
    real:: w_gauss = 3.0
    real:: p_gauss = 4.0
    real:: u_max = 1.0
  contains
    procedure:: init
    procedure:: update
    procedure:: log_density
    procedure:: log_pressure
    procedure:: velocity_magnitude
    procedure:: magnetic_field_magnitude
  end type
  type(operator_t), pointer:: d
CONTAINS

!> =============================================================================
!> Map variables to the current time slot in memory, and initialize the deriv
!> operators.
!> =============================================================================
SUBROUTINE init(self, position, size)
  class(pde_t):: self
  real:: position(3), size(3)
  type(tensor_t):: tensor
  integer:: n(3)
  integer:: i, j
  real:: x, y
  !.............................................................................
  call trace_begin ('hd_mod:pde_t%init', 2)
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
  allocate (d)
  call d%init (self%gn, self%size/self%n)
  !-----------------------------------------------------------------------------
  ! Physical initial conditions: Gaussian top hat in entropy.  Only at time=0.0
  !-----------------------------------------------------------------------------
  self%nv = 5
  if (self%time==0.0) then
    call self%patch_t%init (position, size)
    self%d = 1.0
    self%p = 0.0
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
SUBROUTINE update (self)
  class(pde_t):: self
  type(scalar_t):: viscosity, p_gas, c_sound, log_density, artificial_pressure, &
                   viscous_dissipation, entropy_per_mass, diffusion_coefficient, &
                   gas_pressure, mass_diffusion
  type(vector_t):: velocity, delta_u, mass_flux, entropy_flux, &
                   diffusion_speed
  type(tensor_t):: reynolds_stress, strain_tensor, viscous_stress, momentum_flux
associate (d=>self%d, s=>self%s, p=>self%p, dd_dt=>self%dddt, ds_dt=>self%dsdt, &
           ds=>self%ds, dp_dt=>self%dpdt, gas_gamma=>self%gamma, nu=>self%nu)
  !.............................................................................
  if (gas_gamma==1.0) then
    p_gas = d
    c_sound = 1.0
  else
    p_gas = d**gas_gamma*exp(entropy_per_mass*(gas_gamma-1.0))
    c_sound = gas_gamma*(p_gas/d)
    c_sound = sqrt(c_sound)
  end if
  !-----------------------------------------------------------------------------
  log_density         = log(d)
  velocity            = velocity_interpolation (p, log_density)
  entropy_per_mass    = s/d
  !-----------------------------------------------------------------------------
  diffusion_speed     = c_sound
  diffusion_speed     = (ds*nu(1))*diffusion_speed
  delta_u             = (ds*nu(2))*velocity
  !-----------------------------------------------------------------------------
  reynolds_stress     = outer_product (velocity, velocity)
  strain_tensor       = symmetrized_derivatives (velocity)
  viscous_stress      = tensor_factor (viscosity, strain_tensor)
  viscous_stress      = viscous_stress + edge_multiply (mass_diffusion, velocity)
  artificial_pressure = square_magnitude (delta_u)
  call                tensor_add_diagonal (viscous_stress, artificial_pressure)
  viscous_dissipation = scalar_tensor_product (strain_tensor, viscous_stress)
  call                tensor_add_diagonal (viscous_stress, gas_pressure)
  momentum_flux       = reynolds_stress + viscous_stress
  !-----------------------------------------------------------------------------
  mass_flux           = p + minus_grad (diffusion_coefficient, log_density)
  entropy_flux        =     minus_grad (diffusion_coefficient, entropy_per_mass)
  !-----------------------------------------------------------------------------
  call minus_div_vector (dd_dt, mass_flux)
  call minus_div_vector (ds_dt, entropy_flux)
  call minus_div_tensor (dp_dt, momentum_flux)
end associate
END SUBROUTINE

!> =============================================================================
!> Functions and subroutines in this block support the hydrodynamics equations
!> in update.
!> =============================================================================
FUNCTION scalar_tensor_product (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(tensor_t):: in1, in2
  out%f = in1%xx * in2%xx + in1%xy * in2%xy + in1%xz * in2%xz &
        + in1%yx * in2%yx + in1%yy * in2%yy + in1%yz * in2%yz &
        + in1%zx * in2%zx + in1%zy * in2%zy + in1%zz * in2%zz
END FUNCTION
!> =============================================================================
SUBROUTINE tensor_add_diagonal (out, in1)
  type(tensor_t):: out
  type(scalar_t):: in1
END SUBROUTINE
!> =============================================================================
FUNCTION outer_product (in1, in2) RESULT (out)
  type(tensor_t):: out
  type(vector_t):: in1, in2
END FUNCTION
!> =============================================================================
FUNCTION tensor_factor (in1, in2) RESULT (out)
  type(tensor_t):: out
  type(scalar_t):: in1
  type(tensor_t):: in2
  out%xx = in1%f*out%xx; out%xy = in1%f*out%xy; out%xz = in1%f*out%xz
  out%yx = in1%f*out%yx; out%yy = in1%f*out%yy; out%yz = in1%f*out%yz
  out%zx = in1%f*out%zx; out%zy = in1%f*out%zy; out%zz = in1%f*out%zz
END FUNCTION
!> =============================================================================
FUNCTION edge_multiply (in1, in2) RESULT (out)
  type(tensor_t):: out
  type(scalar_t):: in1
  type(vector_t):: in2
END FUNCTION
!> =============================================================================
FUNCTION symmetrized_derivatives (in1) RESULT (out)
  type(tensor_t):: out
  type(vector_t):: in1
  out%xx = d%ddxup(in1%x)
  out%yy = d%ddyup(in1%y)
  out%zz = d%ddzup(in1%z)
  out%xy = 0.5*(d%ddydn(in1%x) + d%ddxdn(in1%y))
  out%yz = 0.5*(d%ddzdn(in1%y) + d%ddydn(in1%z))
  out%zx = 0.5*(d%ddxdn(in1%z) + d%ddzdn(in1%x))
  out%yx = out%xy
  out%zy = out%yz
  out%xz = out%zx
END FUNCTION
!> =============================================================================
FUNCTION minus_grad (c, in) RESULT (out)
  type(scalar_t):: c, in
  type(vector_t):: out
  out%x = - d%ddxdn(in%f)
  out%y = - d%ddydn(in%f)
  out%z = - d%ddzdn(in%f)
END FUNCTION
!> =============================================================================
SUBROUTINE minus_div_vector (out,in)
  type(scalar_t):: out
  type(vector_t):: in
  out%f = - d%ddxup(in%x) - d%ddyup(in%y) - d%ddzup(in%z)
END SUBROUTINE
!> =============================================================================
SUBROUTINE minus_div_tensor (out,in)
  type(vector_t):: out
  type(tensor_t):: in
  out%x = - d%ddxup(in%xx) - d%ddyup(in%yx) - d%ddzup(in%zx)
  out%y = - d%ddxup(in%xy) - d%ddyup(in%yy) - d%ddzup(in%zy)
  out%z = - d%ddxup(in%xz) - d%ddyup(in%yz) - d%ddzup(in%zz)
END SUBROUTINE
!> =============================================================================
FUNCTION velocity_interpolation (in1, in2) RESULT (out)
  type(vector_t):: in1
  type(scalar_t):: in2
  type(vector_t):: out
  out%x = in1%x/exp(d%xdn(alog(in2%f)))
  out%y = in1%y/exp(d%xdn(alog(in2%f)))
  out%z = in1%z/exp(d%xdn(alog(in2%f)))
END FUNCTION
!> =============================================================================
FUNCTION square_magnitude (in1) RESULT (out)
  type(vector_t):: in1
  type(scalar_t):: out
  out%f = in1%x**2 + in1%y**2 + in1%z**2
END FUNCTION

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
