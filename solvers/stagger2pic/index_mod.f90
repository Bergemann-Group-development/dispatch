!> =============================================================================
!> This index file has slot indices that correspond to the original mhd_mod
!> =============================================================================
MODULE index_mod
  implicit none
  private
  type, public:: index_t
    integer:: d=1, s=2, e=2, px=3, py=4, pz=5, bx=6, by=7, bz=8, phi=9, tt=10
  end type
END MODULE index_mod
