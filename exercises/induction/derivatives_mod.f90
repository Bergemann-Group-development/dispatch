!> =============================================================================
!> $Id$
!> =============================================================================
MODULE derivatives_mod
  USE mpi_send_mod
  USE mpi_coords_mod
  implicit none
  private
  type, public:: derivatives_t
    real:: ds(3)=1.0
    integer:: n(3)
  contains
    procedure:: dx=>void, dy=>void, dz=>void
  end type
CONTAINS

FUNCTION void (self, f)
  class(derivatives_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: f, void
END FUNCTION

END MODULE derivatives_mod

