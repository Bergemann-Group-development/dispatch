!> =============================================================================
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE scalar_mod
  USE vector_mod
  USE mhd_mod
  USE forcing_mod
  USE initial_mod
  USE timestep_mod
  implicit none
  private
  real, parameter:: pi2=8.0*atan(1.0)
  type, public, extends(mhd_t):: experiment_t
    integer:: i_d=1
    real:: G
    type(force_t):: force
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
  real, save:: B0(3)=[0.0,0.0,1.0], d0=1.0
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init')
  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  self%nv = 8
  call self%patch_t%init
  call self%set(bits%static)
  !----------------------------------------------------------------------------
  ! Initial conditions
  !----------------------------------------------------------------------------
  call self%initial%init
  self%mem(:,:,:,:,1,1) = self%initial%selected(self%mesh)
  call self%force%init
  call trace_end
END SUBROUTINE init

!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  type(timestep_t):: timestep
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call self%mhd_t%update
  block
    real, dimension(:,:,:), pointer:: px
    px => self%mem(:,:,:,3,self%it,2)
  call trace('before force')
    call self%check_error(px,"before force")
  self%dpdt = self%dpdt + self%force%selected (self%time, self%mesh)
    call self%check_error(px,"after force")
    px => self%mem(:,:,:,3,self%it,1)
  call trace('after force')
  call self%courant_condition
  call timestep%update (self%id, self%iit, self%dt, self%mem)
    call self%check_error(px,"after timestep")
  end block
  call self%patch_t%update
  call trace_end
END SUBROUTINE update

END MODULE
