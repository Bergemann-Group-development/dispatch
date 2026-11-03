!===============================================================================
!> Code units of the run: CGS with fixed length, time and mass scales.
!===============================================================================
MODULE scaling_mod
  USE units_mod
  USE io_mod
  USE trace_mod
  USE standard_scaling_mod
  implicit none
  private
  !===============================================================================
  !> standard_scaling_t with fixed scales and no output of its own.
  !===============================================================================
  type, public, extends(standard_scaling_t):: scaling_t
  contains
    procedure:: init
    procedure:: output
  end type
  type(scaling_t), public:: scaling
CONTAINS

!===============================================================================
!> Sets the code units once per process: CGS with a length of 1e8 cm, a velocity
!> of 1e6 cm s^-1 and a density of 1e-7 g cm^-3, which make the time unit 100 s
!> and the mass unit 1e17 g.
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
!> Writes nothing: the code units are constant, so the scaling output of
!> standard_scaling_t is not needed.
!===============================================================================
SUBROUTINE output (self)
  class(scaling_t):: self
  !----------------------------------------------------------------------------
END SUBROUTINE output

END MODULE scaling_mod
