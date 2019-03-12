!===============================================================================
!> This index file has slot indices that correspond to the original mhd_mod (and
!> hence to the original HLLC solver in the RAMSES code).
!===============================================================================
MODULE index_mod
  USE io_unit_mod
  implicit none
  private
  type, public:: index_t
    integer:: d=1, e=5, px=2, py=3, pz=4, bx=-1, by=-1, bz=-1, phi=6, tt=-1, s=-1, et=-1, dphi=-1
  contains
    procedure:: output
  end type
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE output (self)
  class(index_t):: self
  integer:: fx=1, fy=2, fz=3
  namelist /idx_nml/ fx, fy, fz
  !-----------------------------------------------------------------------------
  write (io_unit%data, idx_nml)
END SUBROUTINE output

END MODULE index_mod
