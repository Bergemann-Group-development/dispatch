!> =============================================================================
!> Boundary conditions for a 2-D box
!> =============================================================================
MODULE pboundary_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE patch_mod
  USE bits_mod
  USE initial_mod, only: func_d0, func_phi0
  implicit none
  private

  type, public:: pboundary_t
  contains
    procedure:: conditions
  end type
  type(pboundary_t), public:: pboundaries

CONTAINS

!> =============================================================================
!> Impose boundary conditions.
!> All boundaries are fixed and given by the exact solution.
!> =============================================================================
SUBROUTINE conditions (self, patch, it)
  class(pboundary_t) :: self
  class(patch_t), target :: patch
  class(mesh_t), pointer :: m1, m2, m3
  integer:: it, ix, iy
  real:: x, y
  !-----------------------------------------------------------------------------
  call trace%begin('pboundary_t%conditions')

  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)

  if (patch%boundaries%is_set(bits%xl)) then
    do iy=m2%lb,m2%ub
      y = m2%centre_nat + m2%r(iy)
      do ix=m1%lb,m1%lo
        x = m1%centre_nat + m1%r(ix)
        patch%mem(ix,iy,m3%lb:m3%ub,patch%idx%d,it,1)   = func_d0(x,y)
        patch%mem(ix,iy,m3%lb:m3%ub,patch%idx%phi,it,1) = func_phi0(x,y)
      end do
    end do
  end if

  if (patch%boundaries%is_set(bits%xu)) then
    do iy=m2%lb,m2%ub
      y = m2%centre_nat + m2%r(iy)
      do ix=m1%uo,m1%ub
        x = m1%centre_nat + m1%r(ix)
        patch%mem(ix,iy,m3%lb:m3%ub,patch%idx%d,it,1)   = func_d0(x,y)
        patch%mem(ix,iy,m3%lb:m3%ub,patch%idx%phi,it,1) = func_phi0(x,y)
      end do
    end do

  end if

  if (patch%boundaries%is_set(bits%yl)) then
    do iy=m2%lb,m2%lo
      y = m2%centre_nat + m2%r(iy)
      do ix=m1%lb,m1%ub
        x = m1%centre_nat + m1%r(ix)
        patch%mem(ix,iy,m3%lb:m3%ub,patch%idx%d,it,1)   = func_d0(x,y)
        patch%mem(ix,iy,m3%lb:m3%ub,patch%idx%phi,it,1) = func_phi0(x,y)
      end do
    end do
  end if

  if (patch%boundaries%is_set(bits%yu)) then
    do iy=m2%uo,m2%ub
      y = m2%centre_nat + m2%r(iy)
      do ix=m1%lb,m1%ub
        x = m1%centre_nat + m1%r(ix)
        patch%mem(ix,iy,m3%lb:m3%ub,patch%idx%d,it,1)   = func_d0(x,y)
        patch%mem(ix,iy,m3%lb:m3%ub,patch%idx%phi,it,1) = func_phi0(x,y)
      end do
    end do
  end if

  call trace%end
END SUBROUTINE conditions

END MODULE pboundary_mod
