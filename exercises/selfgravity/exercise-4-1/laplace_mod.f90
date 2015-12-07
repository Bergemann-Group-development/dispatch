!===============================================================================
!> Different iteration methods for solving the Poisson equation
!> =============================================================================
MODULE laplace_mod
  USE mpi_mod
  USE io_mod
  USE trace_mod
  USE patches_mod
  implicit none
  private
  type, public:: laplace_t
    integer :: li(3), ui(3), n(3), gn(3)
    real    :: ds(3)
    real    :: omega_sor
    procedure(jacobi_2nd), pointer :: solver => null()
  contains
    procedure:: init
    procedure:: jacobi_2nd
    procedure:: jacobi_4th
  end type
  logical :: first_time = .true.
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
  case('jacobi_2nd')
    self%solver => jacobi_2nd
  case('jacobi_4th')
    self%solver => jacobi_4th
  case default
    write(*,*) 'Poisson solver '//trim(solver_name)//' is unknown'
    write(*,*) 'Only valid solvers are: jacobi_2nd and jacobi_4th'
  end select
  first_time = .false.
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
  real    :: a,b,b1,b2,b3
  !.............................................................................
  ! Make one iterative update for
  !
  !   (2/dx^2 + 2/dy^2 + 2/dz^2) phi_new_ijk =
  !                                 (phi_i+1,j,k + phi_i-1,j,k)/dx^2
  !                               + (phi_i,j+1,k + phi_i,j-1,k)/dy^2
  !                               + (phi_i,j,k+1 + phi_i,j,k-1)/dz^2
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
  ! Add up contributions. First deal with bona fide 3D
  if (all(n > 1)) then
    b1 = a / ds(1)**2
    b2 = a / ds(2)**2
    b3 = a / ds(3)**2
    do k=l(3),u(3)
    do j=l(2),u(2)
    do i=l(1),u(1)
      phi_new(i,j,k) = phi_new(i,j,k) + b1 * (phi(i+1,j,k)+phi(i-1,j,k)) &
                     +                  b2 * (phi(i,j+1,k)+phi(i,j-1,k)) &
                     +                  b2 * (phi(i,j,k+1)+phi(i,j,k-1))
    enddo
    enddo
    enddo
  else ! If's are needed to cover Laplace for all possible 1D and 2D patch configurations
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
  endif
  ! Calculate relative error using old array as scratch
  !error = error_2nd(self, phi_new, phi, FourPiGrho)
  error = error_direct (self, phi_new, phi, FourPiGrho, a)
end associate
END SUBROUTINE jacobi_2nd
!===============================================================================
! 5pt 4th order Jacobi iteration, for filtering high-k
!
! On output phi is destroyed in the error calculation while
! phi_new contains the new guess for phi
!===============================================================================
SUBROUTINE jacobi_4th (self, phi, phi_new, FourPiGrho, error)
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
  real    :: a,b,b1,b2,b3
  !.............................................................................
  ! Make one iterative update for phi_new.
  !
  ! The SECOND ORDER equation is
  !
  !   (2/dx^2 + 2/dy^2 + 2/dz^2) phi_new_ijk =
  !                                 (phi_i+1,j,k + phi_i-1,j,k)/dx^2
  !                               + (phi_i,j+1,k + phi_i,j-1,k)/dy^2
  !                               + (phi_i,j,k+1 + phi_i,j,k-1)/dz^2
  !                               - FourPiGrho
  !
  ! WHAT IS THE FOURTH ORDER EQUIVALENT?
  !
  ! ****************************** HINT HINT HINT ******************************
  ! Check slides and Absalon page
  !
  ! Write down the form of the fourth order operator on paper BEFORE touching code
  !
  ! You are free to assume the you are solving it in the 3D case.
  !
  ! For an extra challenge: What should be done to handle 1D, 2D, and 3D, like in
  ! the 2nd order version? Split up terms in "face", "edge", and "corner" contributions.
  ! ****************************** HINT HINT HINT ******************************
  !
  ! Calculation has to be done in the interior of the domain, assuming a source 4piGrho, and that
  ! the boundaries are fixed and therefore copied
  !.............................................................................
associate(n=>self%n,ds=>self%ds,l=>self%li,u=>self%ui)
  !
  if (all(n > 1)) then
    ! ****************************** HINT HINT HINT ******************************
    ! Insert a, b1, b2, b3 are for fourth order case
    ! ****************************** HINT HINT HINT ******************************

    a = 1. ! WRITE CORRECT PREFACTOR

    ! Init iteration and copy boundary conditions
    phi_new = -a * FourPiGrho

    phi_new(     1:l(1)-1,:,:) = phi(     1:l(1)-1,:,:)
    phi_new(u(1)+1:n(1),  :,:) = phi(u(1)+1:n(1),  :,:)
    phi_new(:,     1:l(2)-1,:) = phi(:,     1:l(2)-1,:)
    phi_new(:,u(2)+1:n(2),  :) = phi(:,u(2)+1:n(2),  :)
    phi_new(:,     1:l(3)-1,:) = phi(:,     1:l(3)-1,:)
    phi_new(:,u(3)+1:n(3),  :) = phi(:,u(3)+1:n(3),  :)

    ! Add up contributions
    b1 = a / ds(1)**2 ! IS THIS RIGHT? ARE THERE MORE TERMS?
    b2 = a / ds(2)**2
    b3 = a / ds(3)**2
    do k=l(3),u(3)
    do j=l(2),u(2)
    do i=l(1),u(1)
      ! THIS IS SECOND ORDER. GENERALIZE TO FOURTH ORDER
      !phi_new(i,j,k) = phi_new(i,j,k) + b1 * (phi(i+1,j,k)+phi(i-1,j,k)) &
      !               +                  b2 * (phi(i,j+1,k)+phi(i,j-1,k)) &
      !               +                  b2 * (phi(i,j,k+1)+phi(i,j,k-1))
    enddo
    enddo
    enddo
  else
    call mpi%abort('Laplace_mod: jacobi_fourth; patch has to be 3D')
    if (n(1)>1) then
      !...
    endif
    if (n(2)>1) then
      !...
    endif
    if (n(3)>1) then
      !...
    endif
  endif
  ! Calculate relative error using old array as scratch
  !error = error_2nd(self, phi_new, phi, FourPiGrho)
  error = error_direct (self, phi_new, phi, FourPiGrho, a)
end associate
END SUBROUTINE jacobi_4th
!===============================================================================
! 3pt 2nd order Error calculation -- scratch is used as temp scratch space
!===============================================================================
FUNCTION error_direct (self, phi, scratch, FourPiGrho, a) result(error)
  implicit none
  class(laplace_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)) :: scratch
  real, dimension(self%gn(1),self%gn(2),self%gn(3)), intent(in)  :: &
                     Phi, FourPiGrho
  real    :: error, a
  !
  integer :: i,j,k
  !.............................................................................
  ! Calculate error as
  !
  !   error = maxval(abs(phi-phi_old)/FourPiGrho)/a
  !
  ! which is mathematically equivalent to maxval(abs(residue)/FourPiGrho)
  !.............................................................................
  !
associate(n=>self%n,ds=>self%ds,l=>self%li,u=>self%ui)
  do k=l(3),u(3)
  do j=l(2),u(2)
  do i=l(1),u(1)
    scratch(i,j,k) = abs(phi(i,j,k)-scratch(i,j,k))/FourPiGrho(i,j,k)
  enddo
  enddo
  enddo
  error = maxval(scratch(l(1):u(1),l(2):u(2),l(3):u(3)))/a
end associate
END FUNCTION error_direct
!===============================================================================
END MODULE laplace_mod
