!===============================================================================
!> Boundary conditions for a planet surface
!===============================================================================
MODULE pboundary_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE kinds_mod
  implicit none
  private

  type, public:: pboundary_t
    ! variables for storing constant (e.g. inflow) boundary states
    real:: P_l, rho_l, e1_l, vx_l, vy_l, vz_l, Bx_l, By_l, Bz_l, &
           P_r, rho_r, e1_r, vx_r, vy_r, vz_r, Bx_r, By_r, Bz_r
  contains
    procedure:: condition
  end type

  type(pboundary_t), public:: pboundary

CONTAINS

!===============================================================================
!> Impose boundary conditions, taking staggering into account
!===============================================================================
SUBROUTINE condition (self, mem, m, lower)
  class(pboundary_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:) :: mem
  class(mesh_t), pointer, dimension(:) :: m
  logical:: lower
  !............................................................................
  integer:: i, j, k
  !-----------------------------------------------------------------------------
  call trace_begin('boundary_t%condition')

  do k=m(3)%lb,m(3)%ub
    do j=m(2)%lb,m(2)%ub
      if (lower) then
        ! outflow conditions/zero gradient
        do i=m(1)%lb,m(1)%lo
          mem(i,j,k,:) = mem(m(1)%li,j,k,:)
        end do
      else
        ! also outflow conditions/zero gradient
        do i=m(1)%uo,m(1)%gn
          mem(i,j,k,:) = mem(m(1)%ui,j,k,:)
        end do
      end if
    end do
  end do

  call trace_end
END SUBROUTINE condition

END MODULE pboundary_mod
