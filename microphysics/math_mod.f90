!===============================================================================
!> Fundamental constants in CGS and SI units
!===============================================================================
MODULE math_mod
  implicit none
  private
  !---------------------------------------------------------------------
  ! Data type holding fundamental constants of nature
  !---------------------------------------------------------------------
  type, public:: math_t
    character(len=16):: name='not set'
    real(kind=8):: pi  =         acos(-1.0_8)
    real(kind=8):: pi2 = 2.0_8 * acos(-1.0_8)
    real(kind=8):: pi4 = 4.0_8 * acos(-1.0_8)
    real(kind=8):: e = exp(1d0)
    real(kind=8):: ln10 = log(10d0)
  contains
    procedure:: init
  end type
  type(math_t), public:: math
CONTAINS

SUBROUTINE init (self)
  class(math_t):: self
  self%name    = 'math'  
END SUBROUTINE init

END MODULE math_mod
