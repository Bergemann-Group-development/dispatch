!> =============================================================================
!> Boundary conditions for a planet surface
!> =============================================================================
MODULE boundary_mod
  USE io_mod
  USE trace_mod
  USE planet_mod
  USE mesh_mod
  implicit none
  private
  type, public:: boundary_t
  contains
    procedure condition
  end type
  type(boundary_t), public:: boundary
CONTAINS

!> =============================================================================
!> Impose boundary conditions, taking staggering into account
!> =============================================================================
SUBROUTINE condition (self, mem, m, h)
  class(boundary_t)      :: self
  real, dimension(:,:,:) :: mem
  class(mesh_t), pointer, dimension(:):: m
  real(8)                :: h(3)
  !............................................................................
  integer      :: ix, iy, iz, iv
  real(8)      :: x, y, z, r
  !-----------------------------------------------------------------------------
  call trace_begin('boundary_t%condition')
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + r3(iz) + h(3)*m(3)%d - planet%position(3)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + r2(iy) + h(2)*m(2)%d - planet%position(2)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + r1(ix) + h(1)*m(1)%d - planet%position(1)
        r = sqrt(x**2+y**2+z**2)
        if (r < planet%radius+m(1)%d*0.1) then
          mem(ix,iy,iz) = 0.0
        end if
      end do
    end do
  end do
  end associate
  call trace_end
END SUBROUTINE condition

END MODULE boundary_mod
