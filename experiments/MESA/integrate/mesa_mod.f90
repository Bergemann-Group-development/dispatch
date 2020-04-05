!===============================================================================
!> MESA module
!===============================================================================
MODULE mesa_mod
  USE io_mod
  USE trace_mod
  USE sinks_mod
  USE gpatch_mod
  implicit none
  private
  type, public, extends(sinks_t):: mesa_t
  contains
    procedure:: accrete
  end type
  type(mesa_t), public:: mesa
CONTAINS

!===============================================================================
!> Overload the accretion procedure in sinks_t, calling it first and then doing
!> MESA related work
!===============================================================================
SUBROUTINE accrete (self, patch)
  class(mesa_t):: self
  class(gpatch_t), pointer:: patch
  !-----------------------------------------------------------------------------
  call trace%begin ('mesa_t%accrete')
  call self%sinks_t%accrete (patch)
  ! ... do whatever ...
  call trace%end ()
END SUBROUTINE accrete

END MODULE mesa_mod
