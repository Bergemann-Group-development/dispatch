!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the advection problem in 3-D, for now with the veloocity field constant
!> in space
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE pde_mod
  USE task_mod
  implicit none
  private
  type, public, extends(task_t):: advection_t
    real, pointer, dimension(:,:,:) :: rho, drhodt
    real :: u(3)
  contains
    procedure:: init
    procedure:: pde
  end type
  type(advection_t):: advect
CONTAINS
!> =============================================================================
!> Make sure the memory is allocated, then map variables to the current time slot
!> The 'it' slot must be the one we are to update, with new values going into the
!> 'new' time slot
!> =============================================================================
SUBROUTINE init (self)
  type(advection_t):: self
  integer:: it
  !.............................................................................
  call self%task_t%init
  self%rho    => task%mem(:,:,:,i_d ,task%it,1)
  self%drhodt => task%mem(:,:,:,i_d ,task%it,2)
END SUBROUTINE init
!> =============================================================================
!> Compute the time derivatives, using the symbolic names mapped into task%mem
!> above
!> =============================================================================
SUBROUTINE pde (self)
  USE derivatives_mod
  type(advection_t):: self
  real, dimension(:,:,:,:), allocatable :: rhou
  integer :: idim
  !.............................................................................
  associate (rho=>self%rho, drhodt=>self%drhodt, u=>self%u, &
             nx=>self%gn(1), ny=>self%gn(2), nz=>self%gn(3))
  allocate(rhou(nx,ny,nz,3))
  do idim=1,3
    rhou(:,:,:,idim) = u(idim)*rho
  enddo
  drhodt = -(d%dx(rhou(:,:,:,1)) + d%dy(rhou(:,:,:,2)) + d%dz(rhou(:,:,:,3)))
  deallocate(rhou)
  end associate
END SUBROUTINE pde
!> =============================================================================
END MODULE pde_mod
