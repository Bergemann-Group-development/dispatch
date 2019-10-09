!===============================================================================
!> Data type for storing the time slices of forces from a patch onto a sink 
!> particle.  The forces for all patches (of relevance) are stored in a list of
!> patch_force_t instances.  To compute the total force from all patches on
!> a sink particle, one loops over the list of patches, using an interpolation
!> procedure patch_force_t%force to compute the contribution to be added from
!> each patch.
!>
!> Note that this data type inherits standard doubly linked list operators, such
!> as %append and %remove (overloaded below)
!===============================================================================
MODULE patch_forces_mod
  USE patch_force_mod
  USE dll_mod
  USE trace_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  type, public, extends(dll_t):: patch_forces_t
    integer:: verbose=0
  contains
    procedure:: init
    procedure:: remove
    procedure:: force
  end type
  type(patch_forces_t), public:: patch_forces
CONTAINS

!===============================================================================
!> Initialize the force list; namelist reading could be added here, if relevant
!===============================================================================
SUBROUTINE init (self, name)
  class(patch_forces_t):: self
  character(len=*), optional:: name
  !.............................................................................
  call self%dll_t%init (name)
END SUBROUTINE init

!===============================================================================
!> Remove one patch force from the list, deallocating its allocated memory
!===============================================================================
SUBROUTINE remove (self, old)
  class(patch_forces_t):: self
  class(dll_node_t), pointer:: old
  !.............................................................................
  call trace%begin ('patch_forces_t%remove')
  select type (old)
  class is (patch_force_t)
  call old%dealloc
  end select
  call self%dll_t%remove (old)
  call trace%end()
END SUBROUTINE remove

!===============================================================================
!> Sum up forces from all patches in the list
!===============================================================================
FUNCTION force (self, time)
  class (patch_forces_t):: self
  real(8):: time
  real:: force(3)
  class (dll_node_t), pointer:: node
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_forces_t%force')
  force = 0.0
  node => self%head
  do while (associated(node))
    select type (node)
    class is (patch_force_t)
    force = force + node%force (time)
if (patch_forces%verbose > 1) &
print *, 'patch_forces_t%force: id, force =', node%id, force
    end select
    node => node%next
  end do
  call trace%end()
END FUNCTION force

END MODULE patch_forces_mod
