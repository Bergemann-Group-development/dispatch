!===============================================================================
!> Boundary conditions for a 2-D box
!===============================================================================
MODULE pboundary_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE patch_mod
  USE bits_mod
  implicit none
  private

  type, public:: pboundary_t
  contains
    procedure:: conditions
  end type
  type(pboundary_t), public:: pboundaries

CONTAINS

!===============================================================================
!> Impose boundary conditions.
!> All conditions are assumed "outflow" (i.e., zero gradient).
!===============================================================================
SUBROUTINE conditions (self, patch, it)
  class(pboundary_t) :: self
  class(patch_t), target :: patch
  class(mesh_t), pointer :: mx, my, mz
  integer :: ix, iy, iz, it
  !-----------------------------------------------------------------------------
  call trace%begin('boundary_t%conditions')

  mx => patch%mesh(1)
  my => patch%mesh(2)
  mz => patch%mesh(3)

  ! this experiment is periodic in the x- and z-direction

  ! y-boundaries
  if (patch%n(2) > 1) then
    ! lower y
    if (patch%boundaries%is_set(bits%yl)) then
      patch%mem(mx%lb:mx%ub,my%li,mz%lb:mz%ub,patch%idx%py,it,1) = 0.0
      do iy=my%lb,my%lo
        patch%mem(mx%lb:mx%ub,iy,mz%lb:mz%ub,:,it,1) = patch%mem(mx%lb:mx%ub,my%li+1,mz%lb:mz%ub,:,it,1)
        ! reflective in the y-component of the momentum
        patch%mem(mx%lb:mx%ub,iy,mz%lb:mz%ub,patch%idx%py,it,1) = -patch%mem(mx%lb:mx%ub,my%li+1,mz%lb:mz%ub,patch%idx%py,it,1)
      end do
    end if
    ! upper y
    if (patch%boundaries%is_set(bits%yu)) then
      patch%mem(mx%lb:mx%ub,my%ui,mz%lb:mz%ub,patch%idx%py,it,1) = 0.0
      do iy=my%uo,my%ub
        patch%mem(mx%lb:mx%ub,iy,mz%lb:mz%ub,:,it,1) = patch%mem(mx%lb:mx%ub,my%ui-1,mz%lb:mz%ub,:,it,1)
        ! reflective in the y-component of the momentum
        patch%mem(mx%lb:mx%ub,iy,mz%lb:mz%ub,patch%idx%py,it,1) = -patch%mem(mx%lb:mx%ub,my%ui-1,mz%lb:mz%ub,patch%idx%py,it,1)
      end do
    end if
  end if

  nullify (mx, my, mz)
  call trace%end()
END SUBROUTINE conditions

END MODULE pboundary_mod
