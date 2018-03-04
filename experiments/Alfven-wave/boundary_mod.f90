!> =============================================================================
!> Boundary conditions for a planet surface
!> =============================================================================
MODULE boundary_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE mhd_mod
  USE index_mod
  implicit none
  private
  type, public:: boundary_t
    integer:: status = 0
    real:: P_l, rho_l, s_l, vx_l, vy_l, vz_l, Bx_l, By_l, Bz_l, et_l, &
           P_r, rho_r, s_r, vx_r, vy_r, vz_r, Bx_r, By_r, Bz_r, et_r
  contains
    procedure:: condition
  end type
  type(boundary_t), public:: boundary
  type(index_t):: idx
CONTAINS

!> =============================================================================
!> Impose boundary conditions, taking staggering into account
!> =============================================================================
SUBROUTINE condition (self, iv, mem, m, lower,idim)
  class(boundary_t)      :: self
  real, dimension(:,:,:) :: mem
  class(mesh_t), pointer, dimension(:) :: m
  integer:: iv,jay1,jay2,jay3,idim
  logical:: lower
  !............................................................................
  integer      :: i1, i2, i3, lb, lo, ub, uo
  !-----------------------------------------------------------------------------
  call trace_begin('boundary_t%condition')

  ! periodic boundary conditions

  if (idim==3) then
    if (lower) then
      mem(:,:,m(3)%lb) = mem(:,:,m(3)%ui-2)
      mem(:,:,m(3)%lo) = mem(:,:,m(3)%ui-1)
!      mem(:,:,m(3)%li) = mem(:,:,m(3)%ui  )
!      do i3=m(3)%lb,m(3)%li
!        jay3 = modulo(i3 - m(3)%li, m(3)%n) + m(3)%li
!        mem(:,:,i3) = mem(:,:,jay3)
!      end do
    else
      mem(:,:,m(3)%uo) = mem(:,:,m(3)%li+1)
      mem(:,:,m(3)%ub) = mem(:,:,m(3)%li+2)
!      do i3=m(3)%uo,m(3)%ub
!        jay3 = modulo(i3 - m(3)%li, m(3)%n) + m(3)%li
!        mem(:,:,i3) = mem(:,:,jay3)
!      end do
    endif
  else if (idim==2) then
    if (lower) then
      mem(:,m(2)%lb,m(3)%li:m(3)%ui) = mem(:,m(2)%ui-2,m(3)%li:m(3)%ui)
      mem(:,m(2)%lo,m(3)%li:m(3)%ui) = mem(:,m(2)%ui-1,m(3)%li:m(3)%ui)
!      mem(:,m(2)%li,:) = mem(:,m(2)%ui  ,:)
!      do i2=m(2)%lb,m(2)%lo
!        jay2 = modulo(i2 - m(2)%li, m(2)%n) + m(2)%li
!        mem(:,i2,:) = mem(:,jay2,:)
!      end do
    else
      mem(:,m(2)%uo,m(3)%li:m(3)%ui) = mem(:,m(2)%li+1,m(3)%li:m(3)%ui)
      mem(:,m(2)%ub,m(3)%li:m(3)%ui) = mem(:,m(2)%li+2,m(3)%li:m(3)%ui)
!      do i2=m(2)%uo,m(2)%ub
!        jay2 = modulo(i2 - m(2)%li, m(2)%n) + m(2)%li
!        mem(:,i2,:) = mem(:,jay2,:)
!      end do
    endif
  else if (idim==1) then
    if (lower) then
      mem(m(1)%lb,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui) = mem(m(1)%ui-2,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui)
      mem(m(1)%lo,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui) = mem(m(1)%ui-1,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui)
!      mem(m(1)%li,:,:) = mem(m(1)%ui  ,:,:)
!      do i1=m(1)%lb,m(1)%lo
!        jay1 = modulo(i1 - m(1)%li, m(1)%n) + m(1)%li
!        mem(i1,:,:) = mem(jay1,:,:)
!      end do
    else
      mem(m(1)%uo,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui) = mem(m(1)%li+1,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui)
      mem(m(1)%ub,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui) = mem(m(1)%li+2,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui)
!      do i1=m(1)%uo,m(1)%ub
!        jay1 = modulo(i1 - m(1)%li, m(1)%n) + m(1)%li
!        mem(i1,:,:) = mem(jay1,:,:)
!      end do
    endif
  endif

!  do iz=m(3)%lb,m(3)%ub
!    do ix=m(1)%lb,m(1)%ub
!      if (lower) then
!        ! outflow conditions
!        mem(ix,m(2)%lb:m(2)%lo,iz) = mem(ix,m(2)%ui,iz)        
!      else
!        ! outflow conditions
!        mem(ix,m(2)%uo:m(2)%ub,iz) = mem(ix,m(2)%ui,iz)
!      end if
!    end do
!  end do
!  do ix=m(1)%lb,m(1)%ub
!    do iy=m(2)%lb,m(2)%ub
!      if (lower) then
!        ! outflow conditions
!        mem(ix,iy,m(3)%lb:m(3)%lo) = mem(ix,iy,m(3)%ui)        
!      else
!        ! outflow conditions
!        mem(ix,iy,m(3)%uo:m(3)%ub) = mem(ix,iy,m(3)%ui)
!      end if
!    end do
!  end do

  call trace_end
END SUBROUTINE condition

END MODULE boundary_mod
