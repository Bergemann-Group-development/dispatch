!===============================================================================
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
    procedure(jacobi_2nd), pointer :: solver => null()
  contains
    procedure:: init
    procedure:: jacobi_2nd
    procedure:: gauss_seidel_2nd
    procedure:: error_2nd
  end type
CONTAINS
!> =============================================================================
SUBROUTINE init (self, patch, solver_name)
  implicit none
  class(laplace_t)  :: self
  class(patch_t)    :: patch
  character(len=12) :: solver_name
  !.............................................................................
  self%n  = patch%n
  self%gn = patch%gn
  self%li = patch%li
  self%ui = patch%ui
  self%ds = patch%ds
  select case(solver_name)
  case('jacobi')
    self%solver => jacobi_2nd
  case('gauss-seidel')
    self%solver => gauss_seidel_2nd
  case default
    write(*,*) 'Poisson solver '//trim(solver_name)//' is unknown'
    write(*,*) 'Only valid solvers are: jacobi and gauss-seidel'
  end select
END SUBROUTINE init
!===============================================================================
! 3pt 2nd order Jacobi iteration, for filtering high-k
!
! On output phi is destroyed in the error calculation while
! phi_new contains the new guess for phi
!===============================================================================
SUBROUTINE jacobi_2nd (self, phi, phi_new, FourPiGrho, error)
  implicit none
  class(laplace_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)), target :: phi
  real, dimension(self%gn(1),self%gn(2),self%gn(3)), intent(out), target :: &
                     phi_new
  real, dimension(self%gn(1),self%gn(2),self%gn(3)), intent(in) :: &
                     FourPiGrho
  real    :: error
  !
  integer :: i,j,k
  real    :: a,b
  !.............................................................................
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
  !.............................................................................
associate(n=>self%n,ds=>self%ds,l=>self%li,u=>self%ui)
  ! Calculate 1 / [(2/dx + 2/dy + 2/dz)] prefactor according to dimensionality
  a = 0.0
  if (n(1)>1) a = a + 2. / ds(1)**2
  if (n(2)>1) a = a + 2. / ds(2)**2
  if (n(3)>1) a = a + 2. / ds(3)**2
  a = 1. / a
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
    b = a / ds(1)**2
    do k=l(3),u(3)
    do j=l(2),u(2)
    do i=l(1),u(1)
      phi_new(i,j,k) = phi_new(i,j,k) + b * (phi(i+1,j,k)+phi(i-1,j,k))
    enddo
    enddo
    enddo
  endif
  if (n(2)>1) then
    b = a / ds(2)**2
    do k=l(3),u(3)
    do j=l(2),u(2)
    do i=l(1),u(1)
      phi_new(i,j,k) = phi_new(i,j,k) + b * (phi(i,j+1,k)+phi(i,j-1,k))
    enddo
    enddo
    enddo
  endif
  if (n(3)>1) then
    b = a / ds(3)**2
    do k=l(3),u(3)
    do j=l(2),u(2)
    do i=l(1),u(1)
      phi_new(i,j,k) = phi_new(i,j,k) + b * (phi(i,j,k+1)+phi(i,j,k-1))
    enddo
    enddo
    enddo
  endif
  ! Calculate relative error using old array as scratch
  error = error_2nd(self, phi_new, phi, FourPiGrho)
end associate
END SUBROUTINE jacobi_2nd
!===============================================================================
! 3pt 2nd order Red-Black Gauss-Seidel iteration, for filtering high-k
!
! On output phi is destroyed in the error calculation while
! phi_new contains the new guess for phi
!===============================================================================
SUBROUTINE gauss_seidel_2nd (self, phi, phi_new, FourPiGrho, error)
  implicit none
  class(laplace_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)), target :: phi
  real, dimension(self%gn(1),self%gn(2),self%gn(3)), intent(out), target :: &
                     phi_new
  real, dimension(self%gn(1),self%gn(2),self%gn(3)), intent(in)  :: &
                     FourPiGrho
  real    :: error
  !
  integer :: i,j,k, sweep
  real    :: a,b
  real, dimension(:,:,:), pointer :: source
  !.............................................................................
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
  !.............................................................................
associate(n=>self%n,ds=>self%ds,l=>self%li,u=>self%ui)
  ! Calculate 1 / [(2/dx + 2/dy + 2/dz)] prefactor according to dimensionality
  a = 0.0
  if (n(1)>1) a = a + 2. / ds(1)**2
  if (n(2)>1) a = a + 2. / ds(2)**2
  if (n(3)>1) a = a + 2. / ds(3)**2
  a = 1. / a
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
  do sweep=0,1
    if (sweep==0) then
      source => phi
    else
      source => phi_new
    endif
    ! Add up contributions
    if (n(1)>1) then
      b = a / ds(1)**2
      do k=l(3),u(3)
      do j=l(2),u(2)
      do i=l(1),u(1)
        if (modulo(i+j+k,2)==sweep) cycle
        phi_new(i,j,k) = phi_new(i,j,k) + b * (source(i+1,j,k)+source(i-1,j,k))
      enddo
      enddo
      enddo
    endif
    if (n(2)>1) then
      b = a / ds(2)**2
      do k=l(3),u(3)
      do j=l(2),u(2)
      do i=l(1),u(1)
        if (modulo(i+j+k,2)==sweep) cycle
        phi_new(i,j,k) = phi_new(i,j,k) + b * (source(i,j+1,k)+source(i,j-1,k))
      enddo
      enddo
      enddo
    endif
    if (n(3)>1) then
      b = a / ds(3)**2
      do k=l(3),u(3)
      do j=l(2),u(2)
      do i=l(1),u(1)
        if (modulo(i+j+k,2)==sweep) cycle
        phi_new(i,j,k) = phi_new(i,j,k) + b * (source(i,j,k+1)+source(i,j,k-1))
      enddo
      enddo
      enddo
    endif
  enddo
  ! Calculate relative error using old array as scratch
  error = error_2nd(self, phi_new, phi, FourPiGrho)
end associate
END SUBROUTINE gauss_seidel_2nd
!===============================================================================
! 3pt 2nd order Error calculation -- scratch is used as temp scratch space
!===============================================================================
FUNCTION error_2nd (self, phi, scratch, FourPiGrho) result(error)
  implicit none
  class(laplace_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)) :: scratch
  real, dimension(self%gn(1),self%gn(2),self%gn(3)), intent(in)  :: &
                     Phi, FourPiGrho
  real    :: error
  !
  integer :: i,j,k
  real    :: b
  !.............................................................................
  ! Calculate error as
  !
  !   error = FourPiGrho 
  !         - (phi_i+1,j,k - 2 phi_i,j,k + phi_i-1,j,k)/dx 
  !         - (phi_i,j+1,k - 2 phi_i,j,k + phi_i,j-1,k)/dy 
  !         - (phi_i,j,k+1 - 2 phi_i,j,k + phi_i,j,k-1)/dz
  !
  ! in the interior of the domain, assuming a source 4piGrho
  !.............................................................................
  error =  0.0
  scratch = FourPiGrho
  !
  ! Add up contributions from Phi's
associate(n=>self%n,ds=>self%ds,l=>self%li,u=>self%ui)
  if (n(1)>1) then
    b = 1. / ds(1)**2
    do k=l(3),u(3)
    do j=l(2),u(2)
    do i=l(1),u(1)
      scratch(i,j,k) = scratch(i,j,k) - b * (phi(i+1,j,k)-2.*phi(i,j,k)+phi(i-1,j,k))
    enddo
    enddo
    enddo
  endif
  if (n(2)>1) then
    b = 1. / ds(2)**2
    do k=l(3),u(3)
    do j=l(2),u(2)
    do i=l(1),u(1)
      scratch(i,j,k) = scratch(i,j,k) - b * (phi(i,j+1,k)-2.*phi(i,j,k)+phi(i,j-1,k))
    enddo
    enddo
    enddo
  endif
  if (n(3)>1) then
    b = 1. / ds(3)**2
  do k=l(3),u(3)
  do j=l(2),u(2)
  do i=l(1),u(1)
      scratch(i,j,k) = scratch(i,j,k) - b * (phi(i,j,k+1)-2.*phi(i,j,k)+phi(i,j,k-1))
  enddo
  enddo
  enddo
  endif  
  error = maxval(abs(scratch(l(1):u(1),l(2):u(2),l(3):u(3))))
  error = error / max(maxval(abs(phi)),maxval(FourPiGrho))
end associate
END FUNCTION error_2nd
!===============================================================================
END MODULE laplace_mod
