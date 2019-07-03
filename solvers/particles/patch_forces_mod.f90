!===============================================================================
!> Data type for storing the fime slices of forces from a patch onto a sink 
!> particle.  The forces for all patches (of relevance) are stored in a list of
!> patch_force_t instances.  Each task should hold a pointer directly to the
!> patch_force_t corresponding patch_force_t instance, to be able to easily
!> update the force history.  To compute the total force from all patches on
!> a sinkparticle, one loops over the list of patches, using an interpolation
!> procedure patch_force_t%interpolate to compute the contribution to be added
!> from each patch.
!===============================================================================
MODULE patch_forces_mod
  USE patch_force_mod
  USE dll_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  type, public, extends(dll_t):: patch_forces_t
  contains
    procedure:: init
    procedure:: force
  end type
  integer, save:: nt=5
CONTAINS

!===============================================================================
!> Initialize storage for one  patch
!===============================================================================
SUBROUTINE init (self, name)
  class(patch_forces_t):: self
  character(len=*), optional:: name
  !.............................................................................
END SUBROUTINE init

!===============================================================================
!> Sum up forces from all patches in the list
!===============================================================================
FUNCTION force (self, time)
  class (patch_forces_t):: self
  real(8):: time
  real:: force(3)
  class (dll_node_t), pointer:: node
  !-----------------------------------------------------------------------------
  force = 0.0
  node => self%head
  do while (associated(node))
    select type (node)
    class is (patch_force_t)
    force = force + node%force (time)
    end select
    node => node%next
  end do
END FUNCTION force

END MODULE patch_forces_mod
