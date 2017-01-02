!> =====================================================================
!> =====================================================================
MODULE rt_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE rt_bundle_mod
  implicit none
  private
  ! --------------------------------------------------------------------
  type, public, extends(task_t):: rt_t
    type(rt_bundle_t)                 :: bundle
    real, dimension(:,:,:),   pointer :: rk
    real, dimension(:,:,:,:), pointer :: src
    real, dimension(:,:,:),   pointer :: dedt
  contains
    procedure :: init
    procedure :: update
    procedure :: solve
    procedure :: interpolate
  end type
! ----------------------------------------------------------------------
CONTAINS

!> =====================================================================
!> =====================================================================
SUBROUTINE init (self)
  class (rt_t) :: self
  ! --------------------------------------------------------------------
  call trace_begin ('rt_t%init')
  call trace_end
END SUBROUTINE init

!> =====================================================================
!> =====================================================================
SUBROUTINE update (self)
  class (rt_t) :: self
  ! --------------------------------------------------------------------
  call trace_begin ('rt_t%update')
  call trace_end
END SUBROUTINE update

!> =====================================================================
!> =====================================================================
SUBROUTINE solve (self)
  class (rt_t) :: self
  ! --------------------------------------------------------------------
  call trace_begin ('rt_t%solve')
  call trace_end
END SUBROUTINE solve

!> =====================================================================
!> =====================================================================
SUBROUTINE interpolate (self, time)
  class (rt_t) :: self
  real(8)      :: time
  ! --------------------------------------------------------------------
  call trace_begin ('rt_t%interpolate')
  call trace_end
END SUBROUTINE interpolate

END MODULE rt_mod
