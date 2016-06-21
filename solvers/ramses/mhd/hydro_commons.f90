module statevar
  use amr_parameters
  use hydro_parameters
  implicit none

  ! 
  ! Enumeration of state variables (second index of uold and unew)
  !

  integer, parameter :: irho=1          ! mass density                             
  integer, parameter :: ipx=2           ! momentum density vector                          
#if NDIM>1
  integer, parameter :: ipy=3
#endif
#if NDIM>2
  integer, parameter :: ipz=4
#endif
  integer, parameter :: ie_tot=ndim+2   ! total energy density
  integer, parameter :: ibxl=ndim+3     ! left magnetic field vector
#if NDIM>1
  integer, parameter :: ibyl=ndim+4
#endif
#if NDIM>2
  integer, parameter :: ibzl=ndim+5                 
#endif
  ! metals are 9 to nvar
  integer, parameter :: ibrx = nvar+1   ! right magnetic field vector
#if NDIM>1
  integer, parameter :: ibry = nvar+2
#endif
#if NDIM>2
  integer, parameter :: ibrz = nvar+3              
#endif
end module

module hydro_commons
  use amr_parameters
  use hydro_parameters
  real(dp),allocatable,dimension(:,:)::uold,unew,usav ! State vector and its update
  real(dp),allocatable,dimension(:)  ::divu,enew ! Non conservative variables
  real(dp)::mass_tot=0.0D0,mass_tot_0=0.0D0,dmax,dmin,daver
  real(dp)::Grho,Grho0,Grho_time,t_Grho               ! Newtons constant and ramp-up params
  type courant_t
    real(dp) :: speed
    real(dp) :: cb, cv, cg, cr, cp
    real(dp) :: x(ndim)
    real(dp) :: urms, rho, umag, bmag, etot
    real(8)  :: volume
    integer  :: level
  end type
  type(courant_t) :: courant_number
  integer::n_switch
  integer:: n_camax, n_cbmax
  integer(kind=8), save:: n_negative_e=0_8, n_camax_tot, n_cbmax_tot, not_iso_tot
  real(dp):: dlnr(1:MAXLEVEL)=0.0_dp
  real(dp):: dlnp(1:MAXLEVEL)=0.0_dp
end module hydro_commons

module const
  use amr_parameters
  real(dp),parameter::bigreal = 1.0e+30
  real(dp),parameter::zero = 0.0
  real(dp),parameter::one = 1.0
  real(dp),parameter::two = 2.0
  real(dp),parameter::three = 3.0
  real(dp),parameter::four = 4.0
  real(dp),parameter::two3rd = 0.6666666666666667
  real(dp),parameter::half = 0.5
  real(dp),parameter::third = 0.33333333333333333
  real(dp),parameter::forth = 0.25
  real(dp),parameter::sixth = 0.16666666666666667
end module const

