!> =============================================================================
!> This index file has slot indices that correspond to the original mhd_mod
!> =============================================================================
MODULE index_mod
  implicit none
  private
  type, public:: index_t
    integer:: d=1, s=5, e=5, px=2, py=3, pz=4, bx=6, by=7, bz=8, phi=9, tt=10, et=-1, dphi=-1
  end type
END MODULE index_mod
