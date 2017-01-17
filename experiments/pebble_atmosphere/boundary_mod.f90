!> =============================================================================
!> Boundary conditions for a planet surface
!> =============================================================================
MODULE boundary_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE planet_mod
  USE mesh_mod
  implicit none
  private
  type, public:: boundary_t
  contains
    procedure:: condition
    procedure:: fixed0
    procedure:: fixed3
    generic:: fixed => fixed0, fixed3
  end type
  type(boundary_t), public:: boundary
CONTAINS

!> =============================================================================
!> Compute stagger offsets, taking care with buggy compilers, which may need 
!> help from specific mesh information.
!> =============================================================================
SUBROUTINE set_stagger (m, h, iv)
  class(mesh_t), pointer :: m(:)
  real(8)                :: h(3)
  integer                :: iv, i
  !-----------------------------------------------------------------------------
  select type (m)
  type is (Cartesian_mesh)
    do i=1,3
      h(i) = m(i)%h(iv)
    end do
  class default
    call mpi%abort ('boundary_mod::set_stagger: staggering h undefined')
  end select
END SUBROUTINE set_stagger

!> =============================================================================
!> Impose fixed boundary condition; scalar values
!> =============================================================================
SUBROUTINE fixed0 (self, mem, v, m, iv)
  class(boundary_t)      :: self
  real, dimension(:,:,:) :: mem
  real                   :: v
  class(mesh_t), pointer :: m(:)
  integer                :: iv
  !............................................................................
  integer                :: ix, iy, iz
  real(8)                :: x, y, z, r, h(3)
  !-----------------------------------------------------------------------------
  call trace_begin('boundary_t%fixed0')
  call set_stagger (m, h, iv)
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + r3(iz) + h(3)*m(3)%d - planet%position(3)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + r2(iy) + h(2)*m(2)%d - planet%position(2)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + r1(ix) + h(1)*m(1)%d - planet%position(1)
        r = sqrt(x**2+y**2+z**2)
        if (r < planet%radius-m(1)%d*0.1) then
          mem(ix,iy,iz) = v
        end if
      end do
    end do
  end do
  end associate
  call trace_end
END SUBROUTINE fixed0

!> =============================================================================
!> Impose fixed boundary condition; array values
!> =============================================================================
SUBROUTINE fixed3 (self, mem, v, m, iv)
  class(boundary_t)      :: self
  real, dimension(:,:,:) :: mem, v
  class(mesh_t), pointer :: m(:)
  integer                :: iv
  !............................................................................
  integer                :: ix, iy, iz
  real(8)                :: x, y, z, r, h(3)
  !-----------------------------------------------------------------------------
  call trace_begin('boundary_t%fixed3')
  call set_stagger (m, h, iv)
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + r3(iz) + h(3)*m(3)%d - planet%position(3)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + r2(iy) + h(2)*m(2)%d - planet%position(2)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + r1(ix) + h(1)*m(1)%d - planet%position(1)
        r = sqrt(x**2+y**2+z**2)
        if (r < planet%radius-m(1)%d*0.1) then
          mem(ix,iy,iz) = v(ix,iy,iz)
        end if
      end do
    end do
  end do
  end associate
  call trace_end
END SUBROUTINE fixed3

!> =============================================================================
!> Skeleton for imposing boundary conditions, taking staggering into account
!> =============================================================================
SUBROUTINE condition (self, mem, m, iv)
  class(boundary_t)      :: self
  real, dimension(:,:,:) :: mem
  class(mesh_t), pointer :: m(:)
  integer                :: iv
  !............................................................................
  integer                :: ix, iy, iz
  real(8)                :: x, y, z, r, h(3)
  !-----------------------------------------------------------------------------
  call trace_begin('boundary_t%condition')
  call set_stagger (m, h, iv)
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
