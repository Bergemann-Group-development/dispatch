!> =============================================================================
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE mhd_mod
  USE initial_mod
  USE timestep_mod
  USE atmosphere_mod
  USE units_mod
  USE scaling_mod
  implicit none
  private
  type, public, extends(mhd_t):: experiment_t
    type(initial_t):: initial
    type(atmosphere_t):: atmosphere
  contains
    procedure init
    procedure update
  end type
  logical, save:: first_time=.true.
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init')
  call scaling%init
  !----------------------------------------------------------------------------
  ! Allocate memory and set IC values
  !----------------------------------------------------------------------------
  self%nv = 8
  call self%initial%init
  call self%patch_t%init
  call self%set(bits%static)
  call self%atmosphere%init
  self%mem(:,:,:,:,1,1) = self%initial%selected(self%mesh)
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  type(timestep_t):: timestep
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call self%mhd_t%update
  call self%courant_condition
  call timestep%update (self%id, self%iit, self%dt, self%mem)
  call trace_end
END SUBROUTINE update

END MODULE
