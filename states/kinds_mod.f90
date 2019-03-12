MODULE kinds_mod
  use iso_fortran_env, only : real32, real64
  implicit none
  private
#ifdef DOUBLE
  integer, public, parameter :: KindScalarVar = real64
#else
  integer, public, parameter :: KindScalarVar = real32
#endif DOUBLE
#ifdef RT_DOUBLE
  integer, public, parameter :: KindRTVar = real64
#else
  integer, public, parameter :: KindRTVar = real32
#endif RT_DOUBLE
  type, public:: kinds_t
  contains
    procedure, nopass:: is_set
  end type
  type(kinds_t), public:: kinds
CONTAINS

!===============================================================================
!> Utility function to handle optional logical switches in parameter lists
!===============================================================================
LOGICAL FUNCTION is_set (switch)
  logical, optional:: switch
  if (present(switch)) then
    is_set = switch
  else
    is_set = .false.
  end if
END  FUNCTION

END MODULE kinds_mod
