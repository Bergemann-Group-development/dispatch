!===============================================================================
!> Experiment module that uses the BIFROST solver_t data type
!===============================================================================
MODULE experiment_mod
  USE io_unit_mod
  USE trace_mod
  USE timestep_mod
  USE solver_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
    procedure:: update
    procedure:: output
  end type
  type(experiment_t), public:: experiment
  character(len=80), save:: id = &
  '$Id$'
CONTAINS

!===============================================================================
!> Interface to let the external code initialize a patch
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%init')
  call self%solver_t%init
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Update the external code with one time step.  Note that this requires the
!> bifrost_update to set the interface_t%u_max value, which is used by the Courant
!> condition (patch_t%courant_condition); cf. code/dispatch_mod.f90
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  real, dimension(:,:,:,:), pointer:: mem
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%update')
  call self%solver_t%update
  call trace%end()
END SUBROUTINE update

!===============================================================================
!> Give the external code a chance to do I/O.  Note that this may not be 
!> necessary, if calling self%patch_t%output is enough.
!===============================================================================
SUBROUTINE output (self, experiment_name)
  class(experiment_t):: self
  character(len=64), optional:: experiment_name
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%output')
  !call self%solver_t%output
  call trace%end()
END SUBROUTINE output

END MODULE experiment_mod
