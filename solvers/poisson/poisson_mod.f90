!> =============================================================================
!> Methods for solving the Poisson equation; derived from routines developed by
!> Troels Haugbølle for the course Computational Astrophysics
!> =============================================================================
MODULE poisson_mod
USE mpi_mod
USE io_mod
USE trace_mod
USE patch_mod
USE mesh_mod
!USE boundary_mod
  implicit none
  private
  type, public:: poisson_t
    real:: omega_sor
    procedure(sor_2nd), pointer :: update => sor_2nd
  contains
    procedure:: init
    procedure:: init_sor
    procedure:: sor_2nd
  end type
  real , save   :: omega_sor=0.0
  logical, save :: first_time = .true.
CONTAINS

!> =============================================================================
!> Initialize poisson_params
!> =============================================================================
SUBROUTINE init (self, patch)
  implicit none
  class(poisson_t)  :: self
  class(patch_t)    :: patch
  character(len=12) :: solver='sor'
  integer           :: iostat
  namelist /poisson_params/ solver, omega_sor
  !.............................................................................
  !$omp critical (poisson_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, poisson_params, iostat=iostat)
    if (io%master) write (io%output, poisson_params)
  end if
  !$omp end critical (poisson_cr)
  select case(trim(solver))
  case('sor')
    self%update => sor_2nd
  case default
    if (io%master) then
      write(*,*) 'Poisson solver '//trim(solver)//' is unknown'
      write(*,*) 'Currently, the only valid solver is: sor'
    end if
  end select
END SUBROUTINE init

!> =============================================================================
!> Initialize SOR factor
!> =============================================================================
SUBROUTINE init_sor (self, patch)
  implicit none
  class(poisson_t)  :: self
  class(patch_t)    :: patch
  !.............................................................................
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
  integer           :: i
  real              :: delta, pref
  real, parameter   :: pi=3.1415926536
  character(len=16) :: comment
  !-----------------------------------------------------------------------------
  if (omega_sor==0.) then
    delta = 0.0; pref = 0.0
    do i=1,3
      if (patch%n(i) > 1) then
        pref  = pref + 1.
        delta = delta + cos(pi / patch%n(i))
      endif
    enddo
    delta = delta / pref
    self%omega_sor = 2.0 / (1. + sqrt(1.0 - delta**2))
    comment =' (computed)'
  else
    self%omega_sor = omega_sor
    comment =' (from namelist)'
  endif
  if (self%omega_sor < 0. .or. self%omega_sor>2) call mpi%abort( &
    'omega_sor in Laplace_mod has to be between 0 and 2 for stability')
  if (first_time .and. io%master) print '(a,f6.4,a)', &
    'poisson_mod: SOR prefactor set to: ', self%omega_sor, comment
END SUBROUTINE init_sor

!> =============================================================================
!> 3pt 2nd order Red-Black Gauss-Seidel iteration with SOR acceleration
!> =============================================================================
SUBROUTINE sor_2nd (self, patch, phi, FourPiGrho, d, error, detailed_timer)
  class(poisson_t) :: self
  class(patch_t)   :: patch
  real, dimension(:,:,:), pointer :: phi
  real, dimension(:,:,:), pointer, intent(in) :: FourPiGrho, d
  real:: error
  logical:: detailed_timer
  !
  integer :: i,j,k,sweep,ss
  real    :: a,b,b1,b2,b3,sor
  real, dimension(:,:,:), pointer :: source
  integer, save          :: itimer=0
associate(n=>patch%n,ds=>patch%ds,l=>patch%li,u=>patch%ui,gn=>patch%gn,omega_sor=>self%omega_sor)
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('poisson_t%sor_2nd_abs', itimer=itimer)
  !.............................................................................
  ! Make one iterative update for
  !
  !   (2/dx^2 + 2/dy^2 + 2/dz^2) phi_new_ijk =
  !                                 (phi_i+1,j,k + phi_i-1,j,k)/dx^2
  !                               + (phi_i,j+1,k + phi_i,j-1,k)/dy^2
  !                               + (phi_i,j,k+1 + phi_i,j,k-1)/dz^2
  !                               - FourPiGrho
  !.............................................................................
  ! Calculate 1 / [(2/dx + 2/dy + 2/dz)] prefactor according to dimensionality
  a = 0.0
  if (n(1)>1) a = a + 2. / ds(1)**2
  if (n(2)>1) a = a + 2. / ds(2)**2
  if (n(3)>1) a = a + 2. / ds(3)**2
  a = 1. / a
  error = 0.0
  do sweep=0,1
    ! Add up contributions
    b1 = a / ds(1)**2
    b2 = a / ds(2)**2
    b3 = a / ds(3)**2
    b = 2.0*(b1+b2+b3)
    do k=l(3),u(3)
      do j=l(2),u(2)
        ss = modulo(j+k+sweep,2)
        do i=l(1)+ss,u(1),2
          phi(i,j,k)     = (1. - omega_sor) * phi(i,j,k) + &
                                 omega_sor * ( &
                                   b1 * (phi(i+1,j,k)+phi(i-1,j,k)) &
                                 + b2 * (phi(i,j+1,k)+phi(i,j-1,k)) &
                                 + b3 * (phi(i,j,k+1)+phi(i,j,k-1)) &
                                 - a * FourPiGrho(i,j,k) )
        end do
        !-----------------------------------------------------------------------
        ! At least with gfortran, keeping these loops separate is faster
        !-----------------------------------------------------------------------
        do i=l(1)+ss,u(1),2
          error = max(error, abs( &
            ((b1*(phi(i+1,j,k)+phi(i-1,j,k)) + &
              b2*(phi(i,j+1,k)+phi(i,j-1,k)) + &
              b3*(phi(i,j,k+1)+phi(i,j,k-1)) - &
              b * phi(i,j,k)) - a*FourPiGrho(i,j,k)))/d(i,j,k))
        end do
      end do
    end do
  end do
  error = error/a
  if (detailed_timer) call trace%end (itimer)
end associate
END SUBROUTINE sor_2nd

END MODULE poisson_mod
