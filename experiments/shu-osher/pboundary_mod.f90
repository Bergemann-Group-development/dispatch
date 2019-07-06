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
    real::  vx_l =  10.141852233 !2.629369 * 3.857143
    real::  p_l  = 10.33333
    real:: gamma=1.4
    real:: vy_l, vz_l, Bx_l, By_l, Bz_l, &
           P_r, rho_r, e1_r, px_r, py_r, pz_r, Bx_r, By_r, Bz_r
  contains
    procedure:: condition
  end type

  type(pboundary_t), public:: pboundary

CONTAINS

!===============================================================================
!> Impose boundary conditions, taking staggering into account
!===============================================================================
SUBROUTINE condition (self, mem, m, idx, lower,solver)
  class(pboundary_t):: self
  character(len=64):: solver
  class(index_t):: idx
  real(kind=KindScalarVar), dimension(:,:,:,:) :: mem
  class(mesh_t), pointer, dimension(:) :: m
  logical:: lower
  !............................................................................
  integer:: i, j, k
  real:: rhov2
  !-----------------------------------------------------------------------------
  call trace_begin('boundary_t%condition')
      if (lower) then
        ! inflow conditions -- constant values
        do i=m(1)%lb,m(1)%lo
          mem(i,:,:,idx%d) = self%rho_l
          mem(i,:,:,idx%px) = self%vx_l * self%rho_l
          mem(i,:,:,idx%py) = 0.0
          mem(i,:,:,idx%pz) = 0.0
          
          if (solver(1:14) == 'zeus_mhd_patch' .or. solver(1:10) == 'stagger2e_') then
            mem(i,:,:,idx%e) = self%P_l / (self%gamma - 1.0)
          else if (solver(1:9) == 'stagger2_') then
            mem(i,:,:,idx%s) = (log(self%P_l) - log(self%rho_l) * self%gamma) * self%rho_l / (self%gamma - 1.0)
          else if (solver(1:6) == 'ramses')  then
            mem(i,:,:,idx%e) = self%P_l / (self%gamma - 1.0)
          end if
        end do  
        if (solver(1:7) /= 'stagger') then
          do i=m(1)%lb,m(1)%lo
            if (solver(1:14) == 'zeus_mhd_patch') then
              rhov2 = (0.25 * (mem(i,1,1,idx%px) + mem(i+1,1,1,idx%px))**2) / mem(i,1,1,idx%d) 
            else if (solver(1:6) == 'ramses') then
              rhov2 = self%vx_l**2 * self%rho_l
            end if
              mem(i,:,:,idx%e) = mem(i,:,:,idx%e) + 0.5 * rhov2
          end do
        end if
      else
        ! also outflow conditions/zero gradient
        do i=m(1)%uo,m(1)%gn
          mem(i,:,:,:) = mem(m(1)%ui,:,:,:)
        end do
      end if
  call trace_end
END SUBROUTINE condition

END MODULE pboundary_mod
