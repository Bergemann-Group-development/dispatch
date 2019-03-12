!===============================================================================
!> This data type is used only to present a common name for all RT solvers
!===============================================================================
MODULE rt_solver_mod
  USE rt_mod
  implicit none
  private
  type, public, extends(rt_t):: rt_solver_t
  end type
CONTAINS
END MODULE rt_solver_mod
