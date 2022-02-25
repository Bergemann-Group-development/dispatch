!===============================================================================
!> $Id$
!> A so-called "abstract type definition", in which parameters may be given,
!> but procedure definitions are "deferred" until another (non-abstract) derived
!> type (object) is defined as an extension of this one
!===============================================================================
MODULE derivatives_mod
  implicit none
  private
  type, abstract, public:: derivatives_t
    real:: ds(3)=1.0
    integer:: n(3)=1
    logical:: do_upwind
  contains
    procedure(deriv), deferred:: dx
    procedure(deriv), deferred:: dy
    procedure(deriv), deferred:: dz
  end type
  abstract interface
    function deriv (self, f)
      import derivatives_t
      class(derivatives_t):: self
      real, dimension(self%n(1), self%n(2), self%n(3)):: f, deriv
    end function
  end interface
END MODULE derivatives_mod

