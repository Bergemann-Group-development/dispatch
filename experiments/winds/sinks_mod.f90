!===============================================================================
!> Sink module stub -- needs only to provide mass and position
!===============================================================================
MODULE sinks_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE scaling_mod
  USE units_mod
  implicit none
  private
  type, public, extends(task_t):: sinks_t
    real(8):: mass=0d0
  contains
    procedure:: init
    procedure:: update
  end type
  real:: mass_msun=30
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self)
  class(sinks_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%init')
  self%mass = mass_msun*cgs%m_sun/scaling%m
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE update (self)
  class(sinks_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%update')
  call trace%end ()
END SUBROUTINE update

END MODULE sinks_mod
