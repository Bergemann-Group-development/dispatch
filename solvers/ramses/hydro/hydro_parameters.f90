module hydro_parameters
  use amr_parameters

  ! Number of independent variables
#ifndef NVAR
  integer,parameter::nvar=ndim+2
#else
  integer,parameter::nvar=NVAR
#endif
  real(dp), parameter:: half=0.0_dp, zero=0.0_dp, one=1.0_dp, two=2.0_dp

  ! Hydro solver parameters
  integer ::niter_riemann=10
  integer ::slope_type=2
  real(dp)::gamma=1.01d0
  real(dp)::courant_factor=0.8d0
  !real(dp)::difmag=0.0d0
  real(dp)::smallc=0.2
  real(dp)::smallr=1.d-3
end module
