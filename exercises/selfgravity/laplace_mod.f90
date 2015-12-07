!> =============================================================================
!> Different iteration methods for solving the Poisson equation
!> =============================================================================
MODULE laplace_mod
  USE io_mod
  USE trace_mod
  USE patches_mod
  implicit none
  private
  type, public:: laplace_t
    integer :: li(3), ui(3), n(3), gn(3)
    real    :: ds(3)
  contains
    procedure:: init
    procedure:: jacobi_2nd
  end type
CONTAINS
!> =============================================================================
SUBROUTINE init (self, patch)
  class(laplace_t):: self
  class(patch_t)  :: patch
  !.............................................................................
  self%n  = patch%n
  self%gn = patch%gn
  self%li = patch%li
  self%ui = patch%ui
  self%ds = patch%ds
END SUBROUTINE init
!=======================================================================
! 3pt 2nd order Jacobi iteration, for filtering high-k
!
! On output phi is destroyed in the error calculation while
! phi_new contains the new guess for phi
!=======================================================================
SUBROUTINE jacobi_2nd (self, phi, phi_new, FourPiGrho, error)
  class(laplace_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)) :: phi
  real, dimension(self%gn(1),self%gn(2),self%gn(3)), intent(out) :: &
                     phi_new
  real, dimension(self%gn(1),self%gn(2),self%gn(3)), intent(in)  :: &
                     FourPiGrho
  real    :: error
  !
  integer :: i,j,k
  real    :: a,b
!.......................................................................
  !
  ! Make one iterative update for
  !
  !   (2/dx + 2/dy + 2/dz) phi_new_ijk =
  !                                 (phi_i+1,j,k + phi_i-1,j,k)/dx 
  !                               + (phi_i,j+1,k + phi_i,j-1,k)/dy 
  !                               + (phi_i,j,k+1 + phi_i,j,k-1)/dz
  !                               - FourPiGrho 
  !
  ! In the interior of the domain, assuming a source 4piGrho, and that
  ! the boundaries are fixed and therefore copied
  !
associate(n=>self%n,ds=>self%ds,l=>self%li,u=>self%ui)
  ! Get 1 / [(2/dx + 2/dy + 2/dz)] prefactor
  a = 0.0
  if (n(1)>1) a = a + 2. / ds(1)
  if (n(2)>1) a = a + 2. / ds(2)
  if (n(3)>1) a = a + 2. / ds(3)
  a = 1. / a
  !.............................................................................
  ! Init iteration and copy boundary conditions
  phi_new = -a * FourPiGrho
  if (n(1)>1) then
    phi_new(     1:l(1)-1,:,:) = phi(     1:l(1)-1,:,:)
    phi_new(u(1)+1:n(1),  :,:) = phi(u(1)+1:n(1),  :,:)
  endif
  if (n(2)>1) then
    phi_new(:,     1:l(2)-1,:) = phi(:,     1:l(2)-1,:)
    phi_new(:,u(2)+1:n(2),  :) = phi(:,u(2)+1:n(2),  :)
  endif
  if (n(3)>1) then
    phi_new(:,     1:l(3)-1,:) = phi(:,     1:l(3)-1,:)
    phi_new(:,u(3)+1:n(3),  :) = phi(:,u(3)+1:n(3),  :)
  endif
  !
  ! Add up contributions
  if (n(1)>1) then
    b = a / ds(1)
    do k=l(3),u(3)
    do j=l(2),u(2)
    do i=l(1),u(1)
      phi_new(i,j,k) = phi_new(i,j,k) + b * (phi(i+1,j,k)+phi(i-1,j,k))
    enddo
    enddo
    enddo
  endif
  if (n(2)>1) then
    b = a / ds(2)
    do k=l(3),u(3)
    do j=l(2),u(2)
    do i=l(1),u(1)
      phi_new(i,j,k) = phi_new(i,j,k) + b * (phi(i,j+1,k)+phi(i,j-1,k))
    enddo
    enddo
    enddo
  endif
  if (n(3)>1) then
    b = a / ds(3)
    do k=l(3),u(3)
    do j=l(2),u(2)
    do i=l(1),u(1)
      phi_new(i,j,k) = phi_new(i,j,k) + b * (phi(i,j,k+1)+phi(i,j,k-1))
    enddo
    enddo
    enddo
  endif
  ! Calculate relative error using old array as scratch
  error =  0.0
  phi = FourPiGrho
  !
  ! Add up contributions from Phi's
  if (n(1)>1) then
    b = 1. / ds(1)
    do k=l(3),u(3)
    do j=l(2),u(2)
    do i=l(1),u(1)
      phi(i,j,k) = phi(i,j,k) - b * (phi_new(i+1,j,k)-2.*phi_new(i,j,k)+phi_new(i-1,j,k))
    enddo
    enddo
    enddo
  endif
  if (n(2)>1) then
    b = 1. / ds(2)
    do k=l(3),u(3)
    do j=l(2),u(2)
    do i=l(1),u(1)
      phi(i,j,k) = phi(i,j,k) - b * (phi_new(i,j+1,k)-2.*phi_new(i,j,k)+phi_new(i,j-1,k))
    enddo
    enddo
    enddo
  endif
  if (n(3)>1) then
    b = 1. / ds(3)
    do k=l(3),u(3)
    do j=l(2),u(2)
    do i=l(1),u(1)
      phi(i,j,k) = phi(i,j,k) - b * (phi_new(i,j,k+1)-2.*phi_new(i,j,k)+phi_new(i,j,k-1))
    enddo
    enddo
    enddo
  endif  
  error = maxval(abs(phi))
  error = error / max(maxval(abs(phi_new)),maxval(FourPiGrho))
end associate
END SUBROUTINE jacobi_2nd
!=======================================================================
END MODULE laplace_mod
