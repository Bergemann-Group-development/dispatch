!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for task
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE patch_mod
  USE task_mod
  USE mesh_mod
  implicit none
  type, extends(tasks_t):: patch_t
    real, dimension(:,:,:,:,:,:), pointer:: mem
    type(mesh_t), dimension(:), pointer:: m
    contains
    procedure:: init
    procedure:: update
    procedure:: allocate
    procedure:: deallocate
    procedure:: print
  end type
CONTAINS

SUBROUTINE init (self, id)
  class(patch_t):: self
  integer:: id
END SUBROUTINE

SUBROUTINE update (self)
  class(patch_t):: self
END SUBROUTINE

SUBROUTINE allocate (self)
  class(patch_t):: self
END SUBROUTINE

SUBROUTINE deallocate (self)
  class(patch_t):: self
END SUBROUTINE

SUBROUTINE print (self)
  class(patch_t):: self
END SUBROUTINE

END MODULE
