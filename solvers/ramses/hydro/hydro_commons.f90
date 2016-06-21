module hydro_commons
  use amr_parameters
  use hydro_parameters
  real(dp),allocatable,dimension(:,:,:,:,:)::uold,unew ! State vector and its update
  real(dp),allocatable,dimension(:,:,:,:,:)::divu,enew ! Non conservative variables
  real(dp)::mass_tot=0.0D0,mass_tot_0=0.0D0
  real(dp)::Grho
  type courant_type
    real(dp) :: cb, cv, cg
    real(dp) :: x(ndim)
    real(dp) :: urms
    real(8)  :: volume
  end type
  type(courant_type) :: courant_number
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

