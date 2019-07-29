!===============================================================================
!> Boundary conditions for radiative heat transfer
!===============================================================================
MODULE rt_boundaries_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE bits_mod
  USE mesh_mod
  USE eos_mod
  USE units_mod
  USE scaling_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! The top temperature should be set from experiment_mod, so it is the same
  ! as the one used there.  For other stars than the Sun, the values of m_sun
  ! and r_sun should be changed to those of the star / planet
  !-----------------------------------------------------------------------------
  type, public :: rt_boundaries_t
    class(mesh_t), pointer:: mesh(:)
    real :: grav, tt_top=4500
  contains
    procedure:: init
    procedure:: condition                        ! modify outer BCs if necessary
  end type
!-------------------------------------------------------------------------------
  type(rt_boundaries_t), public :: rt_boundaries
CONTAINS

!===============================================================================
!> Initialize the boundary data type
!===============================================================================
SUBROUTINE init (self, mesh)
  class(rt_boundaries_t):: self
  class(mesh_t), pointer, optional:: mesh(:)
  !.............................................................................
  self%mesh => mesh
  self%grav = scaling%grav*(cgs%m_sun/scaling%m)/(cgs%r_sun/scaling%l)**2
END SUBROUTINE init

!===============================================================================
!> Vertical boundary conditions; here we explicitly assume the vertical direction
!> to be downwards, so the lower boundary index-wise is the top physical boundary,
!> where I=0 => Q=-S.
!===============================================================================
SUBROUTINE condition (self, mem, src, rk, tt, mu)
  class(rt_boundaries_t) :: self
  real, dimension(:,:,:,:),pointer:: mem, src, rk
  real, dimension(:,:,:), pointer:: tt
  real :: mu(3)
  logical:: lower
  !.............................................................................
  class(mesh_t), pointer:: m1, m2, m3
  integer:: ib, i3, i3b, i2, i1, n(4)
  real:: h_scale, tau
  !-----------------------------------------------------------------------------
  n=shape(mem)
  if (self%mesh(3)%lower_boundary .and. mu(3) > 0.0) then
    i3 = self%mesh(3)%lo
    do ib=1,n(4)
      do i2=self%mesh(2)%lb,self%mesh(2)%ub
      do i1=self%mesh(1)%lb,self%mesh(1)%ub
        ! -- Estimate the boundary Q from an estimate of the total optical depth
        h_scale = tt(i1,i2,i3+1)/(scaling%temp*self%grav)
        tau = h_scale*rk(i1,i2,i3+1,ib)/abs(mu(3))
        mem(i1,i2,i3,ib) = -src(i1,i2,i3,ib)*exp(-tau)
      end do
      end do
    end do
  end if
  if (self%mesh(3)%upper_boundary .and. mu(3) < 0.0) then
    mem(:,:,self%mesh(3)%uo:self%mesh(3)%ub,:) = 0.0
  end if 
END SUBROUTINE condition

END MODULE rt_boundaries_mod
