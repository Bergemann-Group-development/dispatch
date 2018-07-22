!> =============================================================================
!> Default index slots for variables.  May be overruled by placing a local
!> index_mod.f90 file in the solver directory, or the experiment directory.
!> =============================================================================
MODULE index_mod
  implicit none
  private
  type, public:: index_t
    integer:: d=1, s=2, e=2, px=3, py=4, pz=5, bx=6, by=7, bz=8, phi=9, tt=10, dphi=11
  contains
    procedure:: write
  end type
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE write (self)
  class(index_t):: self
  integer:: d, e, s, px, py, pz, bx, by, bz, phi
  namelist /index_nml/ d, e, s, px, py, pz, bx, by, bz, phi
  d   = self%d   - 1
  e   = self%e   - 1
  s   = self%s   - 1
  px  = self%px  - 1
  py  = self%py  - 1
  pz  = self%pz  - 1
  bx  = self%bx  - 1
  by  = self%by  - 1
  bz  = self%bz  - 1
  tt  = self%tt  - 1
  phi = self%phi - 1
  write (io_unit%data, index_nml)
END SUBROUTINE write

END MODULE index_mod
