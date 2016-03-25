!> =============================================================================
!> Setup the initial conditions and experiment update
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE pebble_patch_mod
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
  type, public, extends(pebble_patch_t):: experiment_t
    type(initial_t):: initial
  contains
    procedure init
    procedure update
  end type
  logical, save:: first_time=.true.
CONTAINS

!> =============================================================================
!> The initialization must be done in the right order.  The rotation speed
!> omega must be known before the initial%init runs, for example.
!>
!> dispatch
!>   setup::init
!>     scaling::init                            ! code units
!>     planet::init                             ! planet size, orbit
!>     atmosphere::init                         ! construct hydrostatic state
!>     rubiks_mod::init                         ! makes all patches
!>       experiment_mod::init                   ! calls other inits
!>         mhd_mod::init                        ! reads gamma
!>         pebbles_mod::init                    ! reads omega
!>         initial_mod::init                    ! selects ICs
!>         initial_mod::selected                ! computes ICs
!>         force_mod::init                      ! initialize forcing
!>         boundary_mod::condition              ! apply BCs to ICs
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  if (self%mem_allocated) return
  call trace_begin('experiment_t%init')
  call scaling%init
  !----------------------------------------------------------------------------
  ! Allocate memory and set IC values
  !----------------------------------------------------------------------------
  self%nv = 8
  self%msplit = 3
  call self%mhd_t%init                        ! reads gamma
  call self%pebble_patch_t%init               ! reads omega
  call self%set(bits%static)
  call self%initial%init
  call self%initial%selected(self%mem, self%mesh, self%omega)
  call force%init
  call boundary%condition (self%mem(:,:,:,3,self%it,1), self%mesh, [-0.5d0,0d0,0d0])
  call boundary%condition (self%mem(:,:,:,4,self%it,1), self%mesh, [0d0,-0.5d0,0d0])
  call boundary%condition (self%mem(:,:,:,5,self%it,1), self%mesh, [0d0,0d0,-0.5d0])
  call self%output
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Choice of update steps for this experiment
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !.............................................................................
  type(timestep_t):: timestep
  !-----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call self%check_density ('before update')
  !-----------------------------------------------------------------------------
  ! Option for transporting pebbles with frozen MHD values, in this case, make
  ! sure to wrap around in the time buffer before freezing, so the value needed
  ! by the pebbles have been computed and stored everywhere.
  !-----------------------------------------------------------------------------
  if (self%only_pebbles) then
    if (self%istep<self%nt+1) then
      call self%mhd_t%update
    end if
    if (self%is_set (bits%has_BC)) then
      call boundary%condition (self%mem(:,:,:,3,self%it,1), self%mesh, [-0.5d0,0d0,0d0])
      call boundary%condition (self%mem(:,:,:,4,self%it,1), self%mesh, [0d0,-0.5d0,0d0])
      call boundary%condition (self%mem(:,:,:,5,self%it,1), self%mesh, [0d0,0d0,-0.5d0])
    end if
    self%mem(:,:,:,:,self%it,2) = 0.0
  else
    call self%mhd_t%update
    if (self%is_set (bits%has_BC)) then
      call boundary%condition (self%mem(:,:,:,1,self%it,2), self%mesh, [0d0,0d0,0d0])
      call boundary%condition (self%mem(:,:,:,2,self%it,2), self%mesh, [0d0,0d0,0d0])
      call boundary%condition (self%mem(:,:,:,3,self%it,2), self%mesh, [-0.5d0,0d0,0d0])
      call boundary%condition (self%mem(:,:,:,4,self%it,2), self%mesh, [0d0,-0.5d0,0d0])
      call boundary%condition (self%mem(:,:,:,5,self%it,2), self%mesh, [0d0,0d0,-0.5d0])
      call boundary%condition (self%mem(:,:,:,3,self%it,1), self%mesh, [-0.5d0,0d0,0d0])
      call boundary%condition (self%mem(:,:,:,4,self%it,1), self%mesh, [0d0,-0.5d0,0d0])
      call boundary%condition (self%mem(:,:,:,5,self%it,1), self%mesh, [0d0,0d0,-0.5d0])
    end if
  end if
  call self%courant_condition
  call self%pebble_patch_t%update_pebbles
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call self%check_density ('after timestep')
  call trace_end
END SUBROUTINE update

END MODULE
