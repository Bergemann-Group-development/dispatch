!> =============================================================================
!> Default index slots for variables.  May be overruled by placing a local
!> index_mod.f90 file in the solver directory, or the experiment directory.
!> =============================================================================
MODULE index_mod
  implicit none
  private
  type, public:: index_t
    integer:: d=1, s=5, e=5, px=2, py=3, pz=4, bx=6, by=7, bz=8, phi=-1, dphi = -1
  end type
END MODULE index_mod
