!===============================================================================
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
!> (not overwrite) the generic scratch_t%force_per_unit_mass or _volume arrays.
!>
!> To comply with that policy, this interface also adds the force, and thus
!> assumes that the generic array is set to zero elsewhere, before each timestep.
!> The proper place to do that is in the extras_t data type, which is where the
!> user chooses which forces to include and combine.
!===============================================================================
MODULE forces_mod
  USE io_mod
  USE trace_mod
  USE gpatch_mod
  USE force_mod
  USE link_mod
  USE index_mod
  USE scratch_mod
  implicit none
  private
  type, public:: forces_t
    type(force_t):: force
    class(gpatch_t), pointer:: patch => null()
    real, dimension(:,:,:), pointer:: Ux=>null(), Uy=>null(), Uz=>null()
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: dealloc
    procedure:: pre_update
    procedure:: post_update
  end type
CONTAINS

!===============================================================================
!> Initialize patch link
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(forces_t):: self
  class(gpatch_t), pointer:: patch
  !.............................................................................
  call trace%begin ('forces_t%pre_init')
  self%patch => patch
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize forces_t interface to force_t
!===============================================================================
SUBROUTINE init (self)
  class(forces_t):: self
  !.............................................................................
  call trace%begin ('forces_t%init')
  call self%force%init (self%patch%kind, self%patch%id, self%patch%mesh)
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Deallocate
!===============================================================================
SUBROUTINE dealloc (self)
  class(forces_t):: self
  !.............................................................................
  call trace%begin ('forces_t%dealloc')
  call self%force%dealloc
  call trace%end()
END SUBROUTINE dealloc

!===============================================================================
!> Use the pre-existing force%selected () function to set or add the force
!===============================================================================
SUBROUTINE pre_update (self)
  class(forces_t):: self
  !.............................................................................
  real, dimension(:,:,:)  , pointer, contiguous :: d
  real, dimension(:,:,:,:), pointer :: p
  class(gpatch_t), pointer:: patch
  type(index_t), pointer:: idx
  integer:: it, m(4)
  logical, save:: first_time=.true.
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('forces_t%pre_update', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Use existing function call, if associated
  !-----------------------------------------------------------------------------
  if (associated(self%force%selected)) then
    patch => self%patch
    m(1:3) = patch%mesh%gn
    m(4) = 3
    call scratch%allocate (scratch%force_per_unit_mass, m)
    idx => patch%idx
    it = patch%it
    d => patch%mem(:,:,:,idx%d,it,1)
    p => patch%mem(:,:,:,idx%px:idx%pz,it,1)
    scratch%force_per_unit_mass = &
      self%force%selected(patch%time, d, p, self%Ux, self%Uy, self%Uz, patch%mesh)
  else if (first_time) then
    !$omp critical (forces_cr)
    write (io%output,*) 'WARNING: no force selected'
    first_time = .false.
    !$omp end critical (forces_cr)
  end if
  call trace%end (itimer)
END SUBROUTINE pre_update

!===============================================================================
!> Stub
!===============================================================================
SUBROUTINE post_update (self)
  class(forces_t)                               :: self
  !.............................................................................
END SUBROUTINE post_update

END MODULE forces_mod
