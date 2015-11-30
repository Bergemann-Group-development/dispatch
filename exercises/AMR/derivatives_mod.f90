!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE derivatives_mod
  USE io_mod
  USE trace_mod
  USE central_mod
  USE upwind_mod
  implicit none
  private
  type, public, extends(upwind_t):: derivatives_t
  contains
    procedure:: init
  end type
  logical:: first_time = .true.
CONTAINS

!> =============================================================================
SUBROUTINE init (self, dims, ds)
  class(derivatives_t):: self
  integer:: dims(3)
  real:: ds(3)
  logical:: do_upwind = .false.
  real:: factor = 1.5
  namelist /deriv_params/ do_upwind, factor
  !.............................................................................
  call trace_begin('derivatives_t%init', 1)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, deriv_params)
    if (io%master) write (*, deriv_params)
  end if
  self%do_upwind = do_upwind
  self%factor = factor
  self%n = dims
  self%ds = ds
  call trace_end
END SUBROUTINE init
!> =============================================================================
END MODULE derivatives_mod
