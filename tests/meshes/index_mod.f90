!> =============================================================================
!> This index file has slot indices that correspond to the original mhd_mod
!> =============================================================================
MODULE index_mod
  implicit none
  private
  type, public:: index_t
    integer:: d=1, e=2, et=3, phi=4, px=5, py=6, pz=7, bx=8, by=9, bz=10
  end type
END MODULE index_mod
