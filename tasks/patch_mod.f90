!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for patches, which adds pointers to memory and mesh, and
!> number of dimensions and variables.  The number of dimensions is hardwired
!> to 3, since we can just choose the actual dimensions to get 2-D and 1-D.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE patch_mod
  USE task_mod
  USE mesh_mod
  implicit none
  integer, parameter:: ndim = 3
  type, extends(tasks_t):: patch_t
    real, dimension(:,:,:,:,:,:), pointer:: mem
    type(mesh_t), dimension(:), pointer:: m
    integer:: ndim = ndim           ! Number of dimensions
    integer:: nvar                  ! Number of variables
  end type
PRIVATE
PUBLIC patch_t, ntime

END MODULE
