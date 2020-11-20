!===============================================================================
!> Module that implements halo loading of tasks, using either guard zone
!> loading or boundary conditions.  This relies on going through an init
!> procedure when creating new nbor lists, which allocates either a guard
!> zone loader or a boundary condition data type, depending on the type of
!> boundary.
!>
!> One simple way to do this initially, without disturbing the existing setup,
!> is to first startup as now, with cartesian_mod creating tasks and with an
!> nbor list attached to each task link.
!>
!> Then, in a 2nd step, one loops through all task links, and all nbors on
!> the existing nbor list, replacing the nbor_t data type with a guard_t
!> data type, and finally adding a boundary_t data type for each face that
!> corresponds to a physical boundary.
!===============================================================================
MODULE halo_mod
  USE link_mod
  USE guard_mod
  USE boundary_mod
  !USE face_mod
  implicit none
  private
  !.............................................................................
  type, public, extends(link_t):: halo_t
    !class(face_t), pointer:: face
  contains
    procedure:: init
    procedure:: update
  end type 
CONTAINS

SUBROUTINE init (self)
  class(halo_t):: self
  class(nbor_t), pointer:: nbor
  !class(face_t), pointer:: face
  !-----------------------------------------------------------------------------
  nbor => self%nbor
  do while (associated(nbor))
    nbor => nbor%next
  end do
END SUBROUTINE init

SUBROUTINE update (self)
  class(halo_t):: self
  class(nbor_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  nbor => self%nbor
  do while (associated(nbor))
    call nbor%update
    nbor => nbor%next
  end do
END SUBROUTINE update

END MODULE halo_mod
