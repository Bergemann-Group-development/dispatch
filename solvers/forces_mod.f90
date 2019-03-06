!> =============================================================================
!> This is a bare minimum interface to the new extras module, showing how to
!> interface with a pre-existing force_mod%selected() function.
!> =============================================================================
MODULE forces_mod
  USE patch_mod
  USE force_mod
  USE kinds_mod
  implicit none
  private
  type, public:: forces_t
    type(force_t):: force
    real, dimension(:,:,:), pointer:: Ux=>null(), Uy=>null(), Uz=>null()
  contains
    procedure:: init
    procedure:: pre_update
  end type
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self, patch)
  class(forces_t):: self
  class(patch_t):: patch
  !.............................................................................
  call self%force%init (patch%kind, patch%id, patch%mesh)
END SUBROUTINE init

!> =============================================================================
!> Use the pre-existing force%selected () function to set or add the force
!> =============================================================================
SUBROUTINE pre_update (self, patch, add)
  class(forces_t):: self
  class(patch_t):: patch
  logical, optional:: add
  !.............................................................................
  real, dimension(:,:,:), pointer:: d
  real, dimension(:,:,:,:), pointer:: p
  !-----------------------------------------------------------------------------
  ! Make sure the LHS is allocated
  !-----------------------------------------------------------------------------
  if (.not.allocated(patch%force_per_unit_mass)) then
    associate (m => patch%gn)
    allocate (patch%force_per_unit_mass(m(1),m(2),m(3),3))
    end associate
  end if
  !-----------------------------------------------------------------------------
  ! Allow several forces, as controlled by calls in extras_mod
  !-----------------------------------------------------------------------------
  d => patch%mem(:,:,:,patch%idx%d,patch%it,1)
  p => patch%mem(:,:,:,patch%idx%px:patch%idx%pz,patch%it,1)
  if (kinds%is_set(add)) then
    patch%force_per_unit_mass = &
    patch%force_per_unit_mass + &
    self%force%selected(patch%time, d, p, self%Ux, self%Uy, self%Uz, patch%mesh)
  !-----------------------------------------------------------------------------
  ! Use existing function call
  !-----------------------------------------------------------------------------
  else
    patch%force_per_unit_mass = &
    self%force%selected(patch%time, d, p, self%Ux, self%Uy, self%Uz, patch%mesh)
  end if
END SUBROUTINE pre_update

END MODULE forces_mod

