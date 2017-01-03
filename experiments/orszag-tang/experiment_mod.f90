!> =============================================================================
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE mhd_mod
  USE initial_mod
  USE timestep_mod
  implicit none
  private
  type, public, extends(mhd_t):: experiment_t
    type(initial_t):: initial
  contains
    procedure init
    procedure update
  end type
  logical, save:: first_time=.true.
CONTAINS

!> =============================================================================
!> Setup initial states in the patches
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init')

  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  self%nv = 8
  call self%mhd_t%init
  call self%set (bits%static)
  !
  !----------------------------------------------------------------------------
  ! Initial conditions
  !----------------------------------------------------------------------------
  call self%initial%init()
  self%mhd = .true.
  call self%initial%init_mhd_vars(self%mesh, self%gamma, self%mem(:,:,:,:,1,1))
  !
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Solver update; the default should be built into the solver
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  type(timestep_t):: timestep
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  !
  call self%mhd_t%update
  call self%courant_condition
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  !
  call trace_end
END SUBROUTINE update

END MODULE
