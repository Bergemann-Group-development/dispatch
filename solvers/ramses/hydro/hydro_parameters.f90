module hydro_parameters
  use amr_parameters

  ! Number of independent variables
#ifndef NVAR
  integer,parameter::nvar=ndim+2
#else
  integer,parameter::nvar=NVAR
#endif

  ! Hydro solver parameters
  integer ::niter_riemann=10
  integer ::slope_type=2
  real(dp)::gamma=1.01d0
  real(dp)::courant_factor=0.8d0
  !real(dp)::difmag=0.0d0
  real(dp)::smallc=1.d-10
  real(dp)::smallr=1.d-10
end module
