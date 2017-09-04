MODULE kinds_mod
  use iso_fortran_env, only : real32, real64
  implicit none
  public

#ifdef DOUBLE
  integer, parameter :: KindScalarVar = real64
#else
  integer, parameter :: KindScalarVar = real32
#endif DOUBLE
#ifdef RT_DOUBLE
  integer, parameter :: KindRTVar = real64
#else
  integer, parameter :: KindRTVar = real32
#endif RT_DOUBLE
contains

END MODULE kinds_mod
