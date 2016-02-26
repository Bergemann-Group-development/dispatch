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
    integer:: i_d=1
    real:: G
    type(initial_t):: initial
  contains
    procedure init
    procedure update
  end type
  logical, save:: first_time=.true.
CONTAINS

!> =============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init')

  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  self%nv = 8
  call self%patch_t%init
  call self%set(bits%static)
  call self%mhd_t%init
  
  !----------------------------------------------------------------------------
  ! Initial conditions
  !----------------------------------------------------------------------------
  call self%initial%init()
  self%mhd = .true.
  self%mem(:,:,:,:,1,1) = self%initial%init_mhd_vars(self%mesh,self%gamma)

  call trace_end
END SUBROUTINE init

!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  type(timestep_t):: timestep
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  
  call self%mhd_t%update
  call self%courant_condition
  call timestep%update (self%id, self%iit, self%dt, self%mem)
  call self%patch_t%update

  call trace_end
END SUBROUTINE update

END MODULE
