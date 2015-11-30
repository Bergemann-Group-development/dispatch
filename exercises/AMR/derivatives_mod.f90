!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE derivatives_mod
  implicit none
  private
  type, abstract, public:: derivatives_t
    real:: ds(3)=1.0, dt, factor
    integer:: n(3)
  contains
    procedure:: init
    procedure(deriv), deferred:: dx, dy, dz
  end type
  abstract interface
    function deriv (self, f)
      import derivatives_t
      class(derivatives_t):: self
      real, dimension(self%n(1), self%n(2), self%n(3)):: f, deriv
    end function
  end interface
CONTAINS
!> =============================================================================
SUBROUTINE init (self, ds)
  class(derivatives_t):: self
  real:: ds(3)
  self%ds = ds
END SUBROUTINE init
!> =============================================================================
END MODULE derivatives_mod

