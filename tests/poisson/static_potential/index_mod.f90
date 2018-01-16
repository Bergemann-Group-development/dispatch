!> =============================================================================
!> This index file has slot indices that correspond to the original mhd_mod
!> =============================================================================
MODULE index_mod
  implicit none
  private
  type, public:: index_t
    integer:: d=1, phi=2, px=-1, py=-1, pz=-1, s=-1, e=-1, bx=-1, by=-1, bz=-1, dphi=-1
  end type
END MODULE index_mod
