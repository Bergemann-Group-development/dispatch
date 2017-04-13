!> =============================================================================
!> The gpatch_t data type extends the patch_t with a generic guard zone dnload
!> procedure, which covers the need for all simple, Cartesian mesh based 
!> solvers, and eliminates the need for a dnload procedure in each one of them.
!> =============================================================================
MODULE gpatch_mod
  USE io_mod
  USE patch_mod
  USE download_mod
  USE initial_mod
  USE kinds_mod
  implicit none
  private
  type, public, extends(patch_t):: gpatch_t
    type(initial_t):: initial
  contains
    procedure:: init
    procedure:: dnload
  end type
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self)
  class(gpatch_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: ff
  logical:: ok
  !-----------------------------------------------------------------------------
  call self%patch_t%init
  call self%initial%init
  ff => self%mem(:,:,:,:,1,1)
  call self%input (ok)
  if (.not.ok) call self%initial%selected (self%mesh, ff)
END SUBROUTINE

!> =============================================================================
!> =============================================================================
SUBROUTINE dnload (self)
  class(gpatch_t):: self
  call download%download_link (self%link, all_cells=.false.)
END SUBROUTINE

END MODULE gpatch_mod
