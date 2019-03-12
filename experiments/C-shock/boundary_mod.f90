!===============================================================================
!> Boundary conditions for a planet surface
!===============================================================================
MODULE boundary_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE index_mod
  implicit none
  private
  type, public:: boundary_t
    integer:: status = 0
    real:: P_l, rho_l, s_l, e_l, vx_l, vy_l, vz_l, Bx_l, By_l, Bz_l, &
           P_r, rho_r, s_r, e_r, vx_r, vy_r, vz_r, Bx_r, By_r, Bz_r
  contains
    procedure:: condition
  end type
  type(boundary_t), public:: boundary
CONTAINS

!===============================================================================
!> Impose boundary conditions, taking staggering into account
!===============================================================================


SUBROUTINE condition (self, iv, mem, m, lower,idim,no_mans_land)
  class(boundary_t)      :: self
  real, dimension(:,:,:) :: mem
  class(mesh_t), pointer, dimension(:) :: m
  integer:: iv,jay1,jay2,jay3,idim
  logical:: lower, no_mans_land
  !............................................................................
  integer      :: i1, i2, i3, lb, lo, ub, uo
  !-----------------------------------------------------------------------------
  call trace_begin('boundary_t%condition')



  if (idim==1) then
    if (lower) then
      if (no_mans_land) then
        ! inflow conditions
        if (iv .eq. idx%d ) mem(m(1)%lb:m(1)%lo,:,:) = self%rho_l
        if (iv .eq. idx%e ) mem(m(1)%lb:m(1)%lo,:,:) = self%e_l
        if (iv .eq. idx%s ) mem(m(1)%lb:m(1)%lo,:,:) = self%s_l
        if (iv .eq. idx%px) mem(m(1)%lb:m(1)%lo,:,:) = self%rho_l*self%vx_l
        if (iv .eq. idx%py) mem(m(1)%lb:m(1)%lo,:,:) = self%rho_l*self%vy_l
        if (iv .eq. idx%pz) mem(m(1)%lb:m(1)%lo,:,:) = self%rho_l*self%vz_l
        if (iv .eq. idx%bx) mem(m(1)%lb:m(1)%lo,:,:) = self%Bx_l
        if (iv .eq. idx%by) mem(m(1)%lb:m(1)%lo,:,:) = self%By_l
        if (iv .eq. idx%bz) mem(m(1)%lb:m(1)%lo,:,:) = self%Bz_l
      else
        ! inflow conditions
        if (iv .eq. idx%d ) mem(m(1)%lb:m(1)%lo,:,:) = self%rho_l
        if (iv .eq. idx%e ) mem(m(1)%lb:m(1)%lo,:,:) = self%e_l
        if (iv .eq. idx%s ) mem(m(1)%lb:m(1)%lo,:,:) = self%s_l
        if (iv .eq. idx%px) mem(m(1)%lb:m(1)%lo,:,:) = self%rho_l*self%vx_l
        if (iv .eq. idx%py) mem(m(1)%lb:m(1)%lo,:,:) = self%rho_l*self%vy_l
        if (iv .eq. idx%pz) mem(m(1)%lb:m(1)%lo,:,:) = self%rho_l*self%vz_l
        if (iv .eq. idx%bx) mem(m(1)%lb:m(1)%lo,:,:) = self%Bx_l
        if (iv .eq. idx%by) mem(m(1)%lb:m(1)%lo,:,:) = self%By_l
        if (iv .eq. idx%bz) mem(m(1)%lb:m(1)%lo,:,:) = self%Bz_l
      endif      
 !     if (no_mans_land) then
 !       mem(m(1)%lo,:,:) = mem(m(1)%li,:,:)
 !       mem(m(1)%lb,:,:) = mem(m(1)%li,:,:)
 !     else
 !       mem(m(1)%lo,:,:) = mem(m(1)%li,:,:)
 !       mem(m(1)%lb,:,:) = mem(m(1)%li,:,:)
 !     endif
    else
      if (no_mans_land) then
        mem(m(1)%uo,:,:) = mem(m(1)%ui,:,:)
        mem(m(1)%ub,:,:) = mem(m(1)%ui,:,:)
      else
        mem(m(1)%uo,:,:) = mem(m(1)%ui,:,:)
        mem(m(1)%ub,:,:) = mem(m(1)%ui,:,:)
      endif
!      do i3=m(3)%uo,m(3)%ub
!        jay3 = modulo(i3 - m(3)%li, m(3)%n) + m(3)%li
!        mem(:,:,i3) = mem(:,:,jay3)
!      end do
    endif
  else if (idim==2) then
    if (lower) then
      if (no_mans_land) then
        mem(m(1)%li:m(1)%ui,m(2)%lb,:) = mem(m(1)%li:m(1)%ui,m(2)%ui-1,:)
        mem(m(1)%li:m(1)%ui,m(2)%lo,:) = mem(m(1)%li:m(1)%ui,m(2)%ui  ,:)
      else
        mem(m(1)%li:m(1)%ui,m(2)%lb,:) = mem(m(1)%li:m(1)%ui,m(2)%ui-2,:)
        mem(m(1)%li:m(1)%ui,m(2)%lo,:) = mem(m(1)%li:m(1)%ui,m(2)%ui-1,:)
      endif

    else
      if (no_mans_land) then
        mem(m(1)%li:m(1)%ui,m(2)%uo,:) = mem(m(1)%li:m(1)%ui,m(2)%li  ,:)
        mem(m(1)%li:m(1)%ui,m(2)%ub,:) = mem(m(1)%li:m(1)%ui,m(2)%li+1,:)
      else
        mem(m(1)%li:m(1)%ui,m(2)%uo,:) = mem(m(1)%li:m(1)%ui,m(2)%li+1,:)
        mem(m(1)%li:m(1)%ui,m(2)%ub,:) = mem(m(1)%li:m(1)%ui,m(2)%li+2,:)
      endif
!      do i2=m(2)%uo,m(2)%ub
!        jay2 = modulo(i2 - m(2)%li, m(2)%n) + m(2)%li
!        mem(:,i2,:) = mem(:,jay2,:)
!      end do
    endif
  else if (idim==3) then
    if (lower) then
      if (no_mans_land) then
        mem(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%lb) = mem(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%ui-1)
        mem(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%lo) = mem(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%ui  )
      else
        mem(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%lb) = mem(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%ui-2)
        mem(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%lo) = mem(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%ui-1)
      endif
!      mem(m(1)%li,:,:) = mem(m(1)%ui  ,:,:)
!      do i1=m(1)%lb,m(1)%lo
!        jay1 = modulo(i1 - m(1)%li, m(1)%n) + m(1)%li
!        mem(i1,:,:) = mem(jay1,:,:)
!      end do
    else
      if (no_mans_land) then
        mem(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%uo) = mem(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li  )
        mem(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%ub) = mem(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li+1)
      else
        mem(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%uo) = mem(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li+1)
        mem(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%ub) = mem(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li+2)
      endif 
    endif
  endif

  call trace_end
END SUBROUTINE condition
END MODULE boundary_mod
