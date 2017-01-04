!> =============================================================================
!> Boundary conditions for a 2-D box
!> =============================================================================
MODULE boundary_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE patch_mod
  implicit none
  private

  type, public:: boundary_t
  contains
    procedure:: conditions
  end type
  type(boundary_t), public:: boundary

CONTAINS

!> =============================================================================
!> Impose boundary conditions.
!> All conditions are assumed "outflow" (i.e., zero gradient).
!> =============================================================================
SUBROUTINE conditions (self, patch, it)
  class(boundary_t) :: self
  class(patch_t), target :: patch
  class(mesh_t), pointer :: mx, my, mz
  integer :: ix, iy, iz, it
  !-----------------------------------------------------------------------------
  call trace_begin('boundary_t%conditions')

  mx => patch%mesh(1)
  my => patch%mesh(2)
  mz => patch%mesh(3)

  ! x-boundaries
  if (patch%n(1) > 1) then
    do iz=mz%lb,mz%ub
      do iy=my%lb,my%ub
        ! lower x
        do ix=mx%lb,mx%lo
          patch%mem(ix,iy,iz,:,it,1) = patch%mem(mx%li+1,iy,iz,:,it,1)
        end do
        ! upper x
        do ix=mx%uo,mx%ub
          patch%mem(ix,iy,iz,:,it,1) = patch%mem(mx%ui-1,iy,iz,:,it,1)
        end do
      end do
    end do
  end if

  ! y-boundaries
  ! NOTE: I am intentionally changing the do-loop structure a little in order to
  ! later compare the speed of the the different implementations.
  if (patch%n(2) > 1) then
    do iz=mz%lb,mz%ub
      ! lower y
      do iy=my%lb,my%lo
        patch%mem(mx%lb:mx%ub,iy,iz,:,it,1) = patch%mem(mx%lb:mx%ub,my%li+1,iz,:,it,1)
      end do
      ! upper y
      do iy=my%uo,my%ub
        patch%mem(mx%lb:mx%ub,iy,iz,:,it,1) = patch%mem(mx%lb:mx%ub,my%ui-1,iz,:,it,1)
      end do
    end do
  end if

  ! z-boundaries
  if (patch%n(3) > 1) then
    ! lower z
    do iz=mz%lb,mz%lo
      patch%mem(mx%lb:mx%ub,my%lb:my%ub,iz,:,it,1) = patch%mem(mx%lb:mx%ub,my%lb:my%ub,mz%li+1,:,it,1)
    end do
    ! upper z
    do iz=mz%uo,mz%ub
      patch%mem(mx%lb:mx%ub,my%lb:my%ub,iz,:,it,1) = patch%mem(mx%lb:mx%ub,my%lb:my%ub,mz%ui-1,:,it,1)
    end do
  end if

  nullify (mx, my, mz)
  call trace_end
END SUBROUTINE conditions

END MODULE boundary_mod
