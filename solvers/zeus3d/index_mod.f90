!> =============================================================================
!> This index file has slot indices that correspond to the original mhd_mod
!> =============================================================================
MODULE index_mod
  implicit none
  private
  type, public:: index_t
    integer:: d=1, e=2, et=3, px=4, py=5, pz=6, bx=7, by=8, bz=9, phi=10, tt=11, s=-1, dphi=-1
  end type
END MODULE index_mod
