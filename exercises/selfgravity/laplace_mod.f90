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
    procedure:: gauss_seidel_2nd
    procedure:: sor_2nd
    procedure:: error_2nd
  end type
  logical :: first_time = .true.
CONTAINS
!> =============================================================================
SUBROUTINE init (self, patch, solver_name, omega_sor)
  implicit none
  class(laplace_t)  :: self
  class(patch_t)    :: patch
  character(len=12) :: solver_name
  real              :: omega_sor
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
  case('sor')
    self%solver => sor_2nd
    ! Calculate optimal guess for the SOR prefactor depending
    ! dimensionality, ds, and grid points n:
    !
    !    delta = (1/sum(1/ds^2)) sum(1/ds^2 * cos(pi/n))
    !    omega_sor = 2 / (1 + (1 - delta^2)^0.5)
    !
    ! For an even dimensionality grid
    !
    !    delta = cos(pi / n)
    !
    ! which for n > 10 is well approximated by
    !
    !    delta ~ 1 - 1/2 (pi/n)**2
    !
    ! and then
    !
    !    omega_sor = 2 / (1 + pi/n)
    !
block                                                                           ! Block can be used for local scope variables
    integer :: i
    real(kind=8)            :: delta, pref
    real(kind=8), parameter :: pi=3.14159265358979323846
    character(len=16)       :: comment
    if (omega_sor==0.) then
      delta = 0.0; pref = 0.0
      do i=1,3
        if (self%n(i) > 1) then
          pref  = pref + 1. / self%ds(i)**2
          delta = delta + 1. / self%ds(i)**2 * cos(pi / self%n(i))
        endif
      enddo
      delta = delta / pref
      self%omega_sor = 2.0 / (1. + sqrt(1.0 - delta**2))
      comment =' (computed)'
    else
      self%omega_sor = omega_sor
      comment =' (from namelist)'
    endif
    if (self%omega_sor < 0. .or. self%omega_sor>2) &
      call mpi%abort('omega_sor in Laplace_mod has to be between 0 and 2 for stability')
    if (first_time) write(*,'(a,f6.4,a)') 'Laplace_mod: SOR prefactor set to: ', self%omega_sor, comment
end block
  case default
    write(*,*) 'Poisson solver '//trim(solver_name)//' is unknown'
    write(*,*) 'Only valid solvers are: jacobi, gauss-seidel, and sor'
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
  ! Add up contributions
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
  else
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
  integer :: i,j,k,sweep,ss
  real    :: a,b,b1,b2,b3
  real, dimension(:,:,:), pointer :: source
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
  do sweep=0,1
    if (sweep==0) then
      source => phi
    else
      source => phi_new
    endif
    ! Add up contributions
    if (all(n > 1)) then
      b1 = a / ds(1)**2
      b2 = a / ds(2)**2
      b3 = a / ds(3)**2
      do k=l(3),u(3)
      do j=l(2),u(2)
        ss = modulo(j+k+sweep,2)                                                ! n(1)>1: We can move sweep factor outside
      do i=l(1)+ss,u(1),2
        phi_new(i,j,k) = phi_new(i,j,k) + b1 * (source(i+1,j,k)+source(i-1,j,k)) & ! The new solution
                                       + b2 * (source(i,j+1,k)+source(i,j-1,k)) &
                                       + b3 * (source(i,j,k+1)+source(i,j,k-1))
      enddo
      enddo
      enddo
    else
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
    endif
  enddo
  ! Calculate relative error using old array as scratch
  !error = error_2nd (self, phi_new, phi, FourPiGrho)
  error = error_direct (self, phi_new, phi, FourPiGrho, a)
end associate
END SUBROUTINE gauss_seidel_2nd
!===============================================================================
! 3pt 2nd order Red-Black Gauss-Seidel iteration with SOR acceleration, for filtering high-k
!
! On output phi is destroyed in the error calculation while
! phi_new contains the new guess for phi
!===============================================================================
SUBROUTINE sor_2nd (self, phi, phi_new, FourPiGrho, error)
  implicit none
  class(laplace_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)), target :: phi
  real, dimension(self%gn(1),self%gn(2),self%gn(3)), intent(out), target :: &
                     phi_new
  real, dimension(self%gn(1),self%gn(2),self%gn(3)), intent(in)  :: &
                     FourPiGrho
  real    :: error
  !
  integer :: i,j,k,sweep,ss
  real    :: a,b,b1,b2,b3,sor
  real, dimension(:,:,:), pointer :: source
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
associate(n=>self%n,ds=>self%ds,l=>self%li,u=>self%ui,omega_sor=>self%omega_sor)
  ! Calculate 1 / [(2/dx + 2/dy + 2/dz)] prefactor according to dimensionality
  a = 0.0
  if (n(1)>1) a = a + 2. / ds(1)**2
  if (n(2)>1) a = a + 2. / ds(2)**2
  if (n(3)>1) a = a + 2. / ds(3)**2
  a = 1. / a
  ! Init iteration and copy boundary conditions
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
    if (all(n > 1)) then
      b1 = a / ds(1)**2
      b2 = a / ds(2)**2
      b3 = a / ds(3)**2
      do k=l(3),u(3)
      do j=l(2),u(2)
        ss = modulo(j+k+sweep,2)                                                ! n(1)>1: We can move sweep factor outside
      do i=l(1)+ss,u(1),2
        phi_new(i,j,k) = (1. - omega_sor) * phi(i,j,k) + &
                               omega_sor * ( &
                                 b1 * (source(i+1,j,k)+source(i-1,j,k)) & ! The new solution
                               + b2 * (source(i,j+1,k)+source(i,j-1,k)) &
                               + b3 * (source(i,j,k+1)+source(i,j,k-1)) &
                               - a * FourPiGrho(i,j,k) )
      enddo
      enddo
      enddo
    else
      ! Correction from last iteration and source term 4 pi G rho
      do k=l(3),u(3)
      do j=l(2),u(2)
      do i=l(1),u(1)
        if (modulo(i+j+k,2)==sweep) cycle                                     ! red-black sweep
        phi_new(i,j,k) = (1 - omega_sor)*phi(i,j,k) &
                            - a * omega_sor * FourPiGrho(i,j,k)
      enddo
      enddo
      enddo
      if (n(1)>1) then ! x-direction of Laplace operator
        b = omega_sor * a / ds(1)**2
        do k=l(3),u(3)
        do j=l(2),u(2)
        do i=l(1),u(1)
          if (modulo(i+j+k,2)==sweep) cycle                                     ! red-black sweep
          phi_new(i,j,k) = phi_new(i,j,k) + b * (source(i+1,j,k)+source(i-1,j,k))
        enddo
        enddo
        enddo
      endif
      if (n(2)>1) then ! y-direction of Laplace operator
        b = omega_sor * a / ds(2)**2
        do k=l(3),u(3)
        do j=l(2),u(2)
        do i=l(1),u(1)
          if (modulo(i+j+k,2)==sweep) cycle
          phi_new(i,j,k) = phi_new(i,j,k) + b * (source(i,j+1,k)+source(i,j-1,k))
        enddo
        enddo
        enddo
      endif
      if (n(3)>1) then ! z-direction of Laplace operator
        b = omega_sor * a / ds(3)**2
        do k=l(3),u(3)
        do j=l(2),u(2)
        do i=l(1),u(1)
          if (modulo(i+j+k,2)==sweep) cycle
          phi_new(i,j,k) = phi_new(i,j,k) + b * (source(i,j,k+1)+source(i,j,k-1))
        enddo
        enddo
        enddo
      endif
      do k=l(3),u(3)
      do j=l(2),u(2)
      do i=l(1),u(1)
        if (modulo(i+j+k,2)==sweep) cycle
        phi_new(i,j,k) = source(i,j,k) + omega_sor*(phi_new(i,j,k) - source(i,j,k))
      enddo
      enddo
      enddo
    endif
  enddo

  ! Calculate relative error using old array as scratch
  !error = error_2nd (self, phi_new, phi, FourPiGrho)
  error = error_direct (self, phi_new, phi, FourPiGrho, a)
end associate
END SUBROUTINE sor_2nd
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
  !         - (phi_i+1,j,k - 2 phi_i,j,k + phi_i-1,j,k)/dx^2
  !         - (phi_i,j+1,k - 2 phi_i,j,k + phi_i,j-1,k)/dy^2
  !         - (phi_i,j,k+1 - 2 phi_i,j,k + phi_i,j,k-1)/dz^2
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
