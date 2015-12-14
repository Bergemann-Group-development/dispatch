!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE rt_mod
  USE io_mod
  USE trace_mod
  USE patches_mod
  implicit none
  private
  type, public:: rt_t
    integer :: li(3), ui(3), n(3), gn(3)
    real(kind=8) :: ds(3)
    procedure(integral), pointer :: solver => null()
  contains
    procedure:: init
    procedure:: integral
    procedure:: feautrier
  end type
  logical :: first_time = .true.
  type(rt_t), public:: rt
CONTAINS
!> =============================================================================
SUBROUTINE init (self, solver_name)
  class(rt_t)       :: self
  character(len=12) :: solver_name
  !.............................................................................
  select case(solver_name)
  case('integral')
    self%solver => integral
  case('feautrier')
    self%solver => feautrier
  end select
END SUBROUTINE
!> =============================================================================
SUBROUTINE solver (self, intensity)
  class(rt_t):: self
  real, pointer, dimension(:,:,:):: intensity
END SUBROUTINE
!> =============================================================================
SUBROUTINE integral (self)
  class(rt_t):: self
END SUBROUTINE
!> =============================================================================
SUBROUTINE feautrier (self)
  class(rt_t):: self
END SUBROUTINE
!> =============================================================================
END MODULE rt_mod
