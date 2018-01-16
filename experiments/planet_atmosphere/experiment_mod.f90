!> =============================================================================
!> Setup the initial conditions and experiment update
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
  USE planet_mod
  USE mesh_mod
  USE boundary_mod
  USE force_mod
  implicit none
  private
  type, public, extends(mhd_t):: experiment_t
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
  real, dimension(:,:,:,:), pointer:: ff
  !----------------------------------------------------------------------------
  if (self%mem_allocated) return
  call trace_begin('experiment_t%init')
  call scaling%init
  !----------------------------------------------------------------------------
  ! Allocate memory and set IC values
  !----------------------------------------------------------------------------
  self%nv = 8
  self%msplit = 3
  call self%mhd_t%init
  call self%set (bits%static)
  call self%initial%init (self%kind)
  ff => self%mem(:,:,:,:,self%it,1)
  call self%initial%condition (self%mesh, ff, self%idx)
  call force%init (self%kind)
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
  call boundary%condition (self%mem(:,:,:,1,self%it,2), self%mesh, [0d0,0d0,0d0])
  call boundary%condition (self%mem(:,:,:,2,self%it,2), self%mesh, [0d0,0d0,0d0])
  call boundary%condition (self%mem(:,:,:,3,self%it,2), self%mesh, [-0.5d0,0d0,0d0])
  call boundary%condition (self%mem(:,:,:,4,self%it,2), self%mesh, [0d0,-0.5d0,0d0])
  call boundary%condition (self%mem(:,:,:,5,self%it,2), self%mesh, [0d0,0d0,-0.5d0])
  call self%courant_condition
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call trace_end
END SUBROUTINE update

END MODULE
