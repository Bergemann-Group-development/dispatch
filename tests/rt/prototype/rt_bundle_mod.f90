!> =====================================================================
!> =====================================================================
MODULE rt_bundle_mod
  USE io_mod
  USE trace_mod
  USE rt_solver_mod
  implicit none
  private
  ! --------------------------------------------------------------------
  type, public:: rt_bundle_t
  contains
    procedure :: init
    procedure :: update
    procedure :: solve
  end type
  type(rt_bundle_t), public:: rt_bundle
! ----------------------------------------------------------------------
CONTAINS

!> =====================================================================
!> =====================================================================
SUBROUTINE init (self)
  class (rt_bundle_t) :: self
  ! --------------------------------------------------------------------
  call trace%begin ('rt_bundle_t%init')
  call rt_solver%init
  call trace%end
END SUBROUTINE init

!> =====================================================================
!> =====================================================================
SUBROUTINE update (self)
  class (rt_bundle_t) :: self
  ! --------------------------------------------------------------------
  call trace%begin ('rt_bundle_t%update')
  call trace%end
END SUBROUTINE update

!> =====================================================================
!> =====================================================================
SUBROUTINE solve (self)
  class (rt_bundle_t) :: self
  ! --------------------------------------------------------------------
  call trace%begin ('rt_bundle_t%solve')
  call trace%end
END SUBROUTINE solve

END MODULE rt_bundle_mod
