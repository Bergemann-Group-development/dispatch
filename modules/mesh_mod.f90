!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for mesh
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE mesh_mod
  implicit none
  !-----------------------------------------------------------------------------
  ! Mesh and memory for HD added
  !-----------------------------------------------------------------------------
  type:: mesh_t
    real(8):: p, s
    real(4), dimension(:), pointer:: r
    integer:: lb, ub, lbi, ubi, lbo, ubo
    contains
    procedure:: init
    procedure:: print
  end type
PRIVATE
PUBLIC mesh_t
CONTAINS
SUBROUTINE init (self)
  class(mesh_t):: self
END SUBROUTINE
SUBROUTINE print (self)
  class(mesh_t):: self
END SUBROUTINE
END MODULE
