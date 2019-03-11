!> =============================================================================
!> This is a bare minimum interface to the new extras module, showing how to
!> interface with a pre-existing force_mod%selected() function.
!>
!> In this version, the force_per_unit_mass array is allocated in %init, and its
!> values are set by a call to the previous standard force_t interface, in order
!> for the new arrangement to be as non-intrusive as possible, relative to pre-
!> existing experiments.
!>
!> In the future, external forces should NOT use this interface.  Instead, each
!> type of force should use its own data type (primarily to hold parameters
!> specific to that force), while the resulting force should be addded to
!> (not overwrite) the generic patch_t%force_per_unit_mass or _volume arrays.
!>
!> To comply with that policy, this interface also adds the force, and thus
!> assumes that the generic array is set to zero elsewhere, before each timestep.
!> The proper place to do that is in the extras_t data type, which is where the
!> user chooses which forces to include and combine.
!> =============================================================================
MODULE forces_mod
  USE patch_mod
  USE force_mod
  USE link_mod
  implicit none
  private
  type, public:: forces_t
    type(force_t):: force
    class(patch_t), pointer:: patch
    real, dimension(:,:,:), pointer:: Ux=>null(), Uy=>null(), Uz=>null()
  contains
    procedure:: init
    procedure:: pre_update
  end type
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self, link)
  class(forces_t):: self
  class(link_t):: link
  !.............................................................................
  self%patch => task2patch (link%task)
  associate (patch => self%patch)
  call self%force%init (patch%kind, patch%id, patch%mesh)
  !-----------------------------------------------------------------------------
  ! Make sure the LHS is allocated
  !-----------------------------------------------------------------------------
  if (.not.allocated(patch%force_per_unit_mass)) then
    associate (m => patch%gn)
    allocate (patch%force_per_unit_mass(m(1),m(2),m(3),3))
    end associate
  end if
  end associate
END SUBROUTINE init

!> =============================================================================
!> Use the pre-existing force%selected () function to set or add the force
!> =============================================================================
SUBROUTINE pre_update (self)
  class(forces_t):: self
  !.............................................................................
  real, dimension(:,:,:), pointer:: d
  real, dimension(:,:,:,:), pointer:: p
  !-----------------------------------------------------------------------------
  ! Use existing function call
  !-----------------------------------------------------------------------------
  associate (patch => self%patch)
  d => patch%mem(:,:,:,patch%idx%d,patch%it,1)
  p => patch%mem(:,:,:,patch%idx%px:patch%idx%pz,patch%it,1)
  patch%force_per_unit_mass = &
  patch%force_per_unit_mass + &
  self%force%selected(patch%time, d, p, self%Ux, self%Uy, self%Uz, patch%mesh)
  end associate
END SUBROUTINE pre_update

END MODULE forces_mod
