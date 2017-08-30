!> =============================================================================
!> This index file has slot indices that correspond to the original mhd_mod (and
!> hence to the original HLLC solver in the RAMSES code).
!> =============================================================================
MODULE index_mod
  implicit none
  private
  type, public:: index_t
    integer:: d=1, e=5, px=2, py=3, pz=4, bx=6, by=7, bz=8, phi=9, tt=10, s=5
  end type
END MODULE index_mod
