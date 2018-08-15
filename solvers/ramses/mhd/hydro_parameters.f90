module hydro_parameters
  use amr_parameters

  ! Number of independent variables
#ifndef NVAR
  integer,parameter::nvar=8
#else
  integer,parameter::nvar=NVAR
#endif

  ! Max Number of Metals
  integer, parameter :: nmetal = nvar-8

  ! Size of hydro kernel
  integer,parameter::iu1=-1
  integer,parameter::iu2=+4
  integer,parameter::ju1=(1-ndim/2)-1*(ndim/2)
  integer,parameter::ju2=(1-ndim/2)+4*(ndim/2)
  integer,parameter::ku1=(1-ndim/3)-1*(ndim/3)
  integer,parameter::ku2=(1-ndim/3)+4*(ndim/3)
  integer,parameter::if1=1
  integer,parameter::if2=3
  integer,parameter::jf1=1
  integer,parameter::jf2=(1-ndim/2)+3*(ndim/2)
  integer,parameter::kf1=1
  integer,parameter::kf2=(1-ndim/3)+3*(ndim/3)

  ! Imposed boundary condition variables
  real(dp),dimension(1:MAXBOUND,1:nvar+3)::boundary_var
  real(dp),dimension(1:MAXBOUND)::d_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::p_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::u_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::v_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::w_bound=0.0d0
  real(dp),dimension(1:nmetal,1:MAXBOUND)::metal_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::A_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::B_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::C_bound=0.0d0

  ! Refinement parameters for hydro
  real(dp),dimension(MAXLEVEL)::err_grad_d=-1.0  ! Density gradient
  real(dp),dimension(MAXLEVEL)::err_grad_s=-1.0  ! Passive scalar gradient
  real(dp),dimension(MAXLEVEL)::err_grad_u= 0.0  ! Velocity gradient
  real(dp),dimension(MAXLEVEL)::err_grad_p=-1.0  ! Pressure gradient
  real(dp),dimension(MAXLEVEL)::err_grad_A=-1.0  ! Bx gradient
  real(dp),dimension(MAXLEVEL)::err_grad_B=-1.0  ! By gradient
  real(dp),dimension(MAXLEVEL)::err_grad_C=-1.0  ! Bz gradient
  real(dp),dimension(MAXLEVEL)::err_grad_B2=-1.0 ! B L2 norm gradient
  real(dp),dimension(MAXLEVEL)::gfloor_u=-1.0    ! global velocity floor
  real(dp),dimension(MAXLEVEL)::gfloor_b=-1.0    ! global magnetic field floor
  real(dp)::floor_d=1.d-10   ! Density floor
  real(dp)::floor_s=1.d-10   ! Passsive scalar floor
  real(dp)::floor_u=1.d-10   ! Velocity floor
  real(dp)::floor_p=1.d-10   ! Pressure floor
  real(dp), dimension(1:nmetal) ::floor_metal=1.d-10  ! Metal floor
  real(dp)::floor_A=1.d-10   ! Bx floor
  real(dp)::floor_B=1.d-10   ! By floor
  real(dp)::floor_C=1.d-10   ! Bz floor
  real(dp)::floor_b2=1.d-10  ! B L2 norm floor
  real(dp)::mass_sph=0.0D0   ! mass_sph
  real(dp)::rho_sph=0.0D0    ! rho_sph -- alternative to mass_sph
  real(dp),dimension(1:MAXLEVEL)::jeans_refine=-1.0
  integer::levelmax_u=99
  integer::levelmax_d=99
  integer::levelmax_p=99
  integer::levelmax_b=99
  integer::h_verbose=0

  ! Initial conditions hydro variables
  real(dp),dimension(1:MAXREGION)::d_region=0.
  real(dp),dimension(1:MAXREGION)::u_region=0.
  real(dp),dimension(1:MAXREGION)::v_region=0.
  real(dp),dimension(1:MAXREGION)::w_region=0.
  real(dp),dimension(1:MAXREGION)::p_region=0.
  real(dp),dimension(1:MAXREGION)::A_region=0.
  real(dp),dimension(1:MAXREGION)::B_region=0.
  real(dp),dimension(1:MAXREGION)::C_region=0.
  real(dp),dimension(1:nmetal,1:MAXREGION)::metal_region=0.

  real(dp),dimension(1:nmetal) :: t_metaldecay=0.                               ! Half-life of metals in years
  character(LEN=40),dimension(1:nmetal) :: metal_table                          ! Yield table filename
  
  ! Hydro solver parameters
  integer ::niter_riemann=10
  integer ::courant_type=0                                                      ! see godunov_utils
  logical:: isothermal=.false.
  real(dp) ::slope_type=1
  real(dp)::limit_slope=0
  real(dp)::gamma=1.4d0
  real(dp)::courant_factor=0.5d0
  real(dp)::dlnr_max=-1.0_dp                                                    ! disabled by default
  real(dp)::dlnp_max=-1.0_dp                                                    ! disabled by default
  real(dp)::smallc=1.d-10
  real(dp)::smallr=1.d-10
  real(dp)::eta_mag=0.0d0
  real(dp)::cmax=0.0,cmax2=0.0                                                  ! current limiter
  real(dp)::cbmax=-1.0                                                          ! max value of cb where we should shift to safe solver
  real(dp)::camax=-1.0                                                          ! max value of cb where we should shift to safe solver
  real(dp)::t_damp=0.0                                                          ! damping time, artificial damping / friction
  character(LEN=10)::scheme='muscl'
  character(LEN=10)::riemann='llf'
  character(LEN=10)::riemann2d='llf'
  logical::read_conservative=.false.
  logical::write_conservative=.false.
  logical::add_passive_scalars=.false.                                          ! if set to true, allow adding passive scalars on restart, if they do not exist
  logical::symmetrize=.false.

  ! Interpolation parameters
  integer ::interpol_var=0
  integer ::interpol_type=1

  ! Passive variables index
  integer::imetal=9
  integer::idelay=9
  integer::ixion=9
  integer::ichem=9
  integer, parameter ::ichem_param=9

  ! Global translated distance
  real(kind=dp) :: r_translate(3,1:MAXLEVEL)=0.0_dp, v_translate(3)=0.0_dp

end module hydro_parameters
