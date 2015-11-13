!> =============================================================================
!> $Id$
!> PDE for simple advection
!> =============================================================================
MODULE pde_mod
  USE stagger_mod
  implicit none
  type pde_t
  contains
    procedure:: dfdt
  end type
CONTAINS

!> =============================================================================
!> The PDE gives the time derivative of the variables.  Anything not set is zero
!> from the allocation of the memory.  Staggered or centered derivatives is set
!> by manipulating the mesh, task%m.
!> =============================================================================
SUBROUTINE dfdt (task)
  class(mhd_t):: task
  type(vector_t):: fd
  !.............................................................................
  fd%x = xdn(task%m,v%d)*v%ux
  fd%y = ydn(task%m,v%d)*v%uy
  task%d%d = - ddxup(task%m,fd%x) - ddyup(task%m fd%y)
END SUBROUTINE

END MODULE experiment_mod
