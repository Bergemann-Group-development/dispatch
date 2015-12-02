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
  type(operator_t), pointer:: d
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
  real:: w_gauss=3.0, p_gauss=4.0, amplitude=2.0
  namelist /sedov_params/ w_gauss, p_gauss, amplitude
  !.............................................................................
  call trace_begin ('hd_mod:pde_t%init', 2)
  rewind (io%input)
  read (io%input, sedov_params)
  if (io%master) write (*,sedov_params)
  self%w_gauss = w_gauss
  self%p_gauss = p_gauss
  self%amplitude = amplitude
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
  print*,self%gn,self%size,self%n
  call d%init (self%gn, self%size/self%n)
  !-----------------------------------------------------------------------------
  ! Physical initial conditions: Gaussian top hat in entropy.  Only at time=0.0
  !-----------------------------------------------------------------------------
  self%nv = 5
  if (self%time==0.0) then
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
!> Evaluates the time derivatives of the hydrodynamic variables
!> =============================================================================
SUBROUTINE update (self)
  class(pde_t):: self
  type(scalar_t):: viscosity, p_gas, c_sound, log_density, artificial_pressure, &
                   viscous_dissipation, entropy_per_mass, gas_pressure
  type(vector_t):: velocity, delta_u, mass_flux, entropy_flux, &
                   diffusion_speed, mass_diffusion
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
  entropy_per_mass    = s/d
  velocity            = velocity_interpolation (p, log_density)                 ! OK
  delta_u             = square_magnitude (velocity)                             ! OK
  self%u_max          = sqrt(maxval(delta_u%f))         ! for Courant condition
  !-----------------------------------------------------------------------------
  diffusion_speed     = c_sound
  diffusion_speed     = (ds*nu(1))*diffusion_speed
  delta_u             = (ds*nu(2))*velocity
  !-----------------------------------------------------------------------------
  mass_diffusion      =       minus_grad (diffusion_speed, log_density)         ! OK
  mass_flux           = p + mass_diffusion
  call                minus_div_vector (dd_dt, mass_flux)                       ! OK
  !-----------------------------------------------------------------------------
  entropy_flux        = p*s + minus_grad (diffusion_speed, entropy_per_mass)    ! OK
  call                minus_div_vector (ds_dt, entropy_flux)                    ! OK
  !-----------------------------------------------------------------------------
  reynolds_stress     = outer_product (p, velocity)                             ! OK
  strain_tensor       = symmetrized_derivatives (velocity)                      ! OK
  viscous_stress      = tensor_factor (log_density, diffusion_speed, strain_tensor) ! OK
  viscous_stress      = viscous_stress + edge_multiply (mass_diffusion, velocity)   ! ??
  artificial_pressure = square_magnitude (delta_u)                              ! OK
  call                tensor_add_diagonal (viscous_stress, artificial_pressure) ! OK
  viscous_dissipation = scalar_tensor_product (strain_tensor, viscous_stress)   ! OK
  call                tensor_add_diagonal (viscous_stress, gas_pressure)        ! OK
  momentum_flux       = reynolds_stress + viscous_stress
  call                minus_div_tensor (dp_dt, momentum_flux)                   ! OK
  !-----------------------------------------------------------------------------
end associate
END SUBROUTINE

!> =============================================================================
!> Functions and subroutines in this block support the hydrodynamics equations
!> in update.
!> =============================================================================
FUNCTION velocity_interpolation (in1, in2) RESULT (out)
  type(vector_t):: in1
  type(scalar_t):: in2
  type(vector_t):: out
  out%x = in1%x/exp(d%xdn(alog(in2%f)))
  out%y = in1%y/exp(d%ydn(alog(in2%f)))
  out%z = in1%z/exp(d%zdn(alog(in2%f)))
END FUNCTION
!> =============================================================================
!> This should produce centered square values
!> =============================================================================
FUNCTION square_magnitude (in1) RESULT (out)
  type(vector_t):: in1
  type(scalar_t):: out
  out%f = d%xup(in1%x**2) + d%yup(in1%y**2) + d%zup(in1%z**2)
END FUNCTION
!> =============================================================================
!> This should produce down-shifted gradients, for face centered fluxes
!> =============================================================================
FUNCTION minus_grad (in1, in2) RESULT (out)
  type(scalar_t):: in2
  type(vector_t):: in1, out
  out%x = - d%xdn(in1%x)*d%ddxdn(in2%f)
  out%y = - d%ydn(in1%y)*d%ddydn(in2%f)
  out%z = - d%zdn(in1%z)*d%ddzdn(in2%f)
END FUNCTION
!> =============================================================================
!> The in argument is downshfted, so the div should have up-derivatives
!> =============================================================================
SUBROUTINE minus_div_vector (out,in)
  type(scalar_t):: out
  type(vector_t):: in
  out%f = - d%ddxup(in%x) - d%ddyup(in%y) - d%ddzup(in%z)
END SUBROUTINE
!> =============================================================================
!> in1%yx is the transport of y-mom by u_x, and is taken ddxup() of for non-
!> diagonal elements, and ddxdn() of for diagonal ones.
!> in2 is velocity, in1 is momentum.   The tensor is symmetric, so only
!> 6 elements are significant.  The off-diagonal elements have d%ddydn(out%xy).
!> FIXME: In IDL the density factor is interpolated separately, which is better
!> =============================================================================
FUNCTION outer_product (in1, in2) RESULT (out)
  type(tensor_t):: out
  type(vector_t):: in1, in2
  out%xx = d%xup(in2%x)*d%xup(in1%x)
  out%yy = d%yup(in2%y)*d%yup(in1%y)
  out%zz = d%zup(in2%z)*d%zup(in1%z)
  out%xy = d%ydn(in2%x)*d%xdn(in1%y)
  out%yz = d%zdn(in2%y)*d%ydn(in1%z)
  out%zx = d%xdn(in2%z)*d%zdn(in1%x)
  out%yx = out%xy
  out%zy = out%yz
  out%xz = out%zx
END FUNCTION
!> =============================================================================
!> These elements should be cell centered on the diagonal, and edge centered
!> elsewhere
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
!> The in1 argument is cell centered, the in2 argument is face-centered, and
!> the in3 argument is edge-centered, as should the output be, except in the
!> diagonal, which is cell centered
!> =============================================================================
FUNCTION tensor_factor (in1, in2, in3) RESULT (out)
  type(tensor_t):: out
  type(scalar_t):: in1
  type(vector_t):: in2
  type(tensor_t):: in3
  out%xx = in1%f*d%xup(in2%x)*in3%xx
  out%yy = in1%f*d%yup(in2%y)*in3%yy
  out%zz = in1%f*d%zup(in2%z)*in3%zz
  out%xy = d%ydn(d%xdn(in1%f)*in2%x)*in3%xy
  out%yz = d%zdn(d%ydn(in1%f)*in2%y)*in3%yz
  out%zx = d%xdn(d%zdn(in1%f)*in2%z)*in3%zx
  out%yx = out%xy
  out%zy = out%yz
  out%xz = out%zx
END FUNCTION
!> =============================================================================
!> This is to compute the cell-centered sum of squares for dissipation
!> =============================================================================
FUNCTION scalar_tensor_product (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(tensor_t):: in1, in2
  out%f = in1%xx*in2%xx + in1%yy*in2%yy + in1%zz*in2%zz + 2*( &
        + d%xup(d%yup(in1%yx*in2%yx)) &
        + d%yup(d%zup(in1%zy*in2%zy)) &
        + d%zup(d%xup(in1%xz*in2%xz)))
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
  out%xx = d%xup(in1%x*in2%x)
  out%yy = d%yup(in1%y*in2%y)
  out%zz = d%zup(in1%z*in2%z)
  out%xy = d%xdn(in1%y)*d%ydn(in2%x)
  out%yz = d%ydn(in1%z)*d%zdn(in2%y)
  out%zx = d%zdn(in1%x)*d%xdn(in2%z)
  out%yx = out%xy
  out%zy = out%yz
  out%xz = out%zx
END FUNCTION
!> =============================================================================
SUBROUTINE minus_div_tensor (out,in)
  type(vector_t):: out
  type(tensor_t):: in
  out%x = - d%ddxup(in%xx) - d%ddyup(in%xy) - d%ddzup(in%xz)
  out%y = - d%ddxup(in%yx) - d%ddyup(in%yy) - d%ddzup(in%yz)
  out%z = - d%ddxup(in%zx) - d%ddyup(in%zy) - d%ddzup(in%zz)
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
