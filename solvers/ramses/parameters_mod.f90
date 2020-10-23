
module parameters_mod
#define NPRE 4
  ! Define real types
  integer,parameter::sp=kind(1.0E0)
#if NPRE==4
  integer,parameter::dp=kind(1.0E0) ! real*4
#else
  integer,parameter::dp=kind(1.0D0) ! real*8
#endif
  logical:: isothermal=.false.
  real(8):: gamma=1.4d0
  real(dp):: smallc=1.d-10
  real(dp):: smallr=1.d-10
  ! from module const
  real(dp), parameter:: zero=0.0_dp
  real(dp), parameter:: one=1.0_dp
  real(dp), parameter:: two=2.0_dp
  real(dp), parameter:: half=0.5_dp
  real(dp), parameter:: forth=0.25_dp
end module parameters_mod
