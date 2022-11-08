!===============================================================================
!> Define code units, in terms of (here) CGS units.  Demonstrates how to pass
!> on the three standard scaling values (l,t,m), and benefit from the standard
!> metadata produced and passed on to Python by self%set (...)
!===============================================================================
MODULE scaling_mod
  USE units_mod
  USE io_mod
  USE trace_mod
  USE standard_scaling_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Extend the units_t data typem which handles scaling metadata output
  !-----------------------------------------------------------------------------
  type, public, extends(standard_scaling_t):: scaling_t
  contains
    procedure:: init
    procedure:: output
  end type
  type(scaling_t), public:: scaling
CONTAINS

!===============================================================================
!> Default scaling module, optionally overloaded with local version.
!===============================================================================
SUBROUTINE init (self)
  class(scaling_t):: self
  real(8), save:: l=1d8, d=1d-7, v=1d6
  character(4), save:: system='CGS'
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  call trace_begin('scaling_t%init')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    call self%units_t%init ()
    call self%set (system, l=l, t=l/v, m=d*l**3)
  end if
  call io%lock%unset
END SUBROUTINE init

!===============================================================================
!> overwrites standard_scaling output
!===============================================================================
SUBROUTINE output (self)
  class(scaling_t):: self
  !----------------------------------------------------------------------------
  ! we don't need this extra output as the code units are constant
END SUBROUTINE output

END MODULE scaling_mod
