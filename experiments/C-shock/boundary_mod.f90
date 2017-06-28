!> =============================================================================
!> Boundary conditions for a planet surface
!> =============================================================================
MODULE boundary_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE index_mod
  implicit none
  private
  type, public:: boundary_t
    integer:: status = 0
    real:: P_l, rho_l, s_l, vx_l, vy_l, vz_l, Bx_l, By_l, Bz_l, &
           P_r, rho_r, s_r, vx_r, vy_r, vz_r, Bx_r, By_r, Bz_r
  contains
    procedure:: condition
  end type
  type(boundary_t), public:: boundary
  type(index_t):: idx
CONTAINS

!> =============================================================================
!> Impose boundary conditions, taking staggering into account
!> =============================================================================
SUBROUTINE condition (self, iv, mem, m, lower)
  class(boundary_t)      :: self
  real, dimension(:,:,:) :: mem
  class(mesh_t), pointer, dimension(:) :: m
  integer:: iv
  logical:: lower
  !............................................................................
  integer      :: ix, iy, iz, lb, lo, ub, uo
  !-----------------------------------------------------------------------------
  call trace_begin('boundary_t%condition')

  lb = m(1)%lb
  lo = m(1)%lo
  uo = m(1)%uo
  ub = m(1)%ub
  do iz=m(3)%lb,m(3)%ub
    do iy=m(2)%lb,m(2)%ub
      if (lower) then
        ! inflow conditions
        if (iv .eq. idx%d ) mem(lb:lo,iy,iz) = self%rho_l
        if (iv .eq. idx%s ) mem(lb:lo,iy,iz) = self%s_l
        if (iv .eq. idx%px) mem(lb:lo,iy,iz) = self%rho_l*self%vx_l
        if (iv .eq. idx%py) mem(lb:lo,iy,iz) = self%rho_l*self%vy_l
        if (iv .eq. idx%pz) mem(lb:lo,iy,iz) = self%rho_l*self%vz_l
        if (iv .eq. idx%bx) mem(lb:lo,iy,iz) = self%Bx_l
        if (iv .eq. idx%by) mem(lb:lo,iy,iz) = self%By_l
        if (iv .eq. idx%bz) mem(lb:lo,iy,iz) = self%Bz_l
      else
        ! outflow conditions
        mem(uo:ub,iy,iz) = mem(m(1)%ui,iy,iz)
      end if
    end do
  end do

  call trace_end
END SUBROUTINE condition

END MODULE boundary_mod
