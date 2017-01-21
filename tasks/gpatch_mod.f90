!> =============================================================================
!> The gpatch_t data type extends the patch_t with a generic guard zone dnload
!> procedure, which covers the need for all simple, Cartesian mesh based 
!> solvers, and eliminates the need for a dnload procedure in each one of them.
!> =============================================================================
MODULE gpatch_mod
  USE patch_mod
  USE download_mod
  implicit none
  private
  type, public, extends(patch_t):: gpatch_t
  contains
    procedure:: dnload
  end type
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE dnload (self)
  class(gpatch_t):: self
  call download%download_link (self%link, all_cells=.false.)
END SUBROUTINE

END MODULE gpatch_mod
