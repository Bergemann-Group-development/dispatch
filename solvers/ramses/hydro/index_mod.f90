!> =============================================================================
!> This index file has slot indices that correspond to the original mhd_mod (and
!> hence to the original HLLC solver in the RAMSES code).
!> =============================================================================
MODULE index_mod
  USE io_unit_mod
  implicit none
  private
  type, public:: index_t
    integer:: d=1, e=5, px=2, py=3, pz=4, bx=-1, by=-1, bz=-1, phi=6, tt=-1, s=-1, et=-1, dphi=-1
  contains
    procedure:: write
  end type
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE write (self)
  class(index_t):: self
  integer:: d, e, s, px, py, pz, bx, by, bz, phi
  namelist /idx_nml/ d, e, s, px, py, pz, bx, by, bz, phi
  d   = self%d   - 1
  e   = self%e   - 1
  s   = self%s   - 1
  px  = self%px  - 1
  py  = self%py  - 1
  pz  = self%pz  - 1
  bx  = self%bx  - 1
  by  = self%by  - 1
  bz  = self%bz  - 1
  phi = self%phi - 1
    write (io_unit%data, idx_nml)
END SUBROUTINE write

END MODULE index_mod
