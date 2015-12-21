!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE rt_mod
  USE io_mod
  USE trace_mod
  USE patches_mod
  implicit none
  private
  interface
    subroutine solve (m, n, dtau, s, ifw1, irv1, q, ifwN, irvN)
      integer:: m, n
      real, dimension(m,n), intent(in)  :: dtau, s
      real, dimension(m,n), intent(out) :: q
      real, dimension(m)  , intent(in)  :: ifw1, irvN
      real, dimension(m)  , intent(out) :: irv1, ifwN
    end subroutine
  end interface
  external integral, feautrier, feautrier8
  logical :: first_time = .true.
  type, public:: rt_t
    procedure(solve), nopass, public, pointer:: solver => null()
  contains
    procedure:: init
  end type
CONTAINS
!> =============================================================================
SUBROUTINE init (self, solver_name)
  class(rt_t):: self
  character(len=12) :: solver_name
  !.............................................................................
  select case(solver_name)
  case('integral')
    self%solver => integral
  case('integral2')
    self%solver => integral
  case('feautrier')
    self%solver => feautrier
  case('feautrier8')
    self%solver => feautrier8
  end select
END SUBROUTINE
!> =============================================================================
END MODULE rt_mod
