!===============================================================================
!> Boundary conditions for a Shu-Osher shock tube experiment
!===============================================================================
MODULE pboundary_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE kinds_mod
  USE index_mod
  implicit none
  private

  type, public:: pboundary_t
    ! variables for storing constant (e.g. inflow) boundary states
    real::   e_l = 10.33333
    real:: rho_l =  3.857143
    real::  px_l =  10.141852233 !2.629369 * 3.857143
    real:: py_l, pz_l, Bx_l, By_l, Bz_l, &
           P_r, rho_r, e1_r, px_r, py_r, pz_r, Bx_r, By_r, Bz_r
  contains
    procedure:: condition
  end type

  type(pboundary_t), public:: pboundary

CONTAINS

!===============================================================================
!> Impose boundary conditions, taking staggering into account
!===============================================================================
SUBROUTINE condition (self, mem, m, idx, lower)
  class(pboundary_t):: self
  class(index_t):: idx
  real(kind=KindScalarVar), dimension(:,:,:,:) :: mem
  class(mesh_t), pointer, dimension(:) :: m
  logical:: lower
  !............................................................................
  integer:: i, j, k
  !-----------------------------------------------------------------------------
  call trace_begin('boundary_t%condition')

      if (lower) then
        ! inflow conditions -- constant values
        do i=m(1)%lb,m(1)%lo
          mem(i,:,:,idx%d) = self%rho_l
          mem(i,:,:,idx%e) = self%e_l
          mem(i,:,:,idx%px) = self%px_l
          mem(i,:,:,idx%py) = 0.0
          mem(i,:,:,idx%pz) = 0.0
          !mem(i,:,:,:) = mem(m(1)%li,:,:,:)
        end do
      else
        ! also outflow conditions/zero gradient
        do i=m(1)%uo,m(1)%gn
          mem(i,:,:,:) = mem(m(1)%ui,:,:,:)
        end do
      end if
  call trace_end
END SUBROUTINE condition

END MODULE pboundary_mod
