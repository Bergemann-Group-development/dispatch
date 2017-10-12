!> =============================================================================
MODULE solver_mod
  USE kinds_mod
  USE gpatch_mod
  
  type, public, extends(gpatch_t):: solver_t
  contains
    procedure:: gas_pressure
    procedure:: gas_velocity_scalar
    procedure:: gas_velocity_vector
    generic:: gas_velocity => gas_velocity_scalar, gas_velocity_vector
  end type

CONTAINS

!> =============================================================================
FUNCTION gas_pressure (self) RESULT (pg)
  class (solver_t):: self
  real(kind=KindScalarVar), dimension(self%gn(1),self%gn(2),self%gn(3)):: pg
  !-----------------------------------------------------------------------------
END FUNCTION gas_pressure

!> =============================================================================
FUNCTION gas_velocity_vector (self) RESULT (v)
  class (solver_t):: self
  real(kind=KindScalarVar), dimension(self%gn(1),self%gn(2),self%gn(3),3):: v
  !-----------------------------------------------------------------------------
END FUNCTION gas_velocity_vector

!> =============================================================================
FUNCTION gas_velocity_scalar (self, idir) RESULT (v)
  class (solver_t):: self
  integer:: idir
  real(kind=KindScalarVar), dimension(self%gn(1),self%gn(2),self%gn(3)):: v
  !-----------------------------------------------------------------------------
END FUNCTION gas_velocity_scalar


END MODULE solver_mod
