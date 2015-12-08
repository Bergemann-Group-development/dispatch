!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the Poisson problem in 3-D
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE ghd_mod
  USE io_mod
  USE omp_timer_mod
  USE trace_mod
  USE patches_mod
  USE scalar_mod
  USE laplace_mod
  USE hd_mod
  USE operators_mod
  implicit none
  private
  integer, parameter:: i_phi=6
  type, public, extends(hd_t):: ghd_t
    type(scalar_t):: phi
    real(kind=8):: G                                                                    ! Strength of gravity
    real(kind=8):: tolerance
    integer:: min_iter, max_iter
    type(laplace_t)  :: laplace
  contains
    procedure:: init
    procedure:: update
    procedure:: solve
  end type
  type(operator_t), pointer:: o
CONTAINS

!> =============================================================================
!> Setup parameters for the Poisson solver
!> =============================================================================
SUBROUTINE init (self)
  class(ghd_t):: self
  real(kind=8) :: tolerance=1e-4, G=1., omega_sor=0.0
  integer      :: min_iter=0, max_iter=100
  character(len=12):: solver='sor'
  logical:: first_time = .true.
  namelist /ghd_params/ G, tolerance, min_iter, max_iter, solver, omega_sor
  !-----------------------------------------------------------------------------
  ! Read standard parameters and copy to structure
  !-----------------------------------------------------------------------------
  call trace_begin ('ghd_mod:ghd_t%init', 2)
  call self%hd_t%init
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, ghd_params)
    if (io%master) write (*,ghd_params)
  end if
  self%size      = 1.
  self%position  = 0.
  self%G         = G
  self%tolerance = tolerance
  self%min_iter  = min_iter
  self%max_iter  = max_iter
  self%nv        = 6
  !-----------------------------------------------------------------------------
  ! Set up the right dimensions from current patch, etc for laplace iteration method
  !-----------------------------------------------------------------------------
  call self%laplace%init (self, solver, omega_sor)
  !-----------------------------------------------------------------------------
  ! Memory mapping of density and potentiall field
  !-----------------------------------------------------------------------------
  call self%phi%init (self%gn, self%size/self%n, 'd')
  if (self%phi%allocated) deallocate(self%phi%f)
  self%phi%f => self%mem(:,:,:,i_phi ,self%it,1)
  self%phi%f = 0.0
  !
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Compute iteratively the solution to the Poisson equation, stopping when
!> reaching the tolerance tol
!> =============================================================================
SUBROUTINE solve (self)
  class(ghd_t):: self
  !
  real(kind=8), allocatable, dimension(:,:,:)         :: FourPiGrho
  real(kind=8), allocatable, dimension(:,:,:), target :: phi_1, phi_2
  real(kind=8), pointer,     dimension(:,:,:)         :: phi_old, phi_new, phi_swap
  real(kind=8), parameter :: FourPi = 4.0 * 3.14159265358979323846
  real(kind=8)            :: error
  integer         :: it
  !.............................................................................
  call trace_begin ('ghd_mod:ghd_t%update', 2)
associate(gn=>self%gn)
  allocate(FourPiGrho(gn(1),gn(2),gn(3)), &
           phi_1(gn(1),gn(2),gn(3)),      &
           phi_2(gn(1),gn(2),gn(3)))            ! allocate scratch memory
end associate
  !-----------------------------------------------------------------------------
  ! Set up source of gravity
  !-----------------------------------------------------------------------------
  FourPiGrho = FourPi * self%G * self%d%f
  phi_1 = self%phi%f
  phi_old => phi_1                              ! initial pointer arrangement
  phi_new => phi_2
  error = 2.*self%tolerance
  it    = 0
  !-----------------------------------------------------------------------------
  ! Do at least min_iter iterations, stop if error<=tolerance or max iterations
  ! reached.  Swap phi_old and phi_new pointers before the next iteration.
  !-----------------------------------------------------------------------------
  do while((it<self%min_iter) .or. &
           (error >= self%tolerance .and. it<self%max_iter))
    call self%laplace%solver (phi_old, phi_new, FourPiGrho, error)
    if (io%verbose>2) write(*,'(a,i4,2(a,e10.3))'), 'Poisson: Iter=', it, &
      ' Tolerance=', self%tolerance, ' Error=', error
    phi_swap => phi_new; phi_new=>phi_old; phi_old=>phi_swap
    it = it + 1
  enddo
  !-----------------------------------------------------------------------------
  ! Check if we reached solution in time
  !-----------------------------------------------------------------------------
  if (it == self%max_iter .and. io%verbose>0) &
    print '(a,i5,2(a,e10.3))', &
     'Poisson: Did not reach convergence in patch. Max iter=', it, &
     ' Tolerance=', self%tolerance, ' Error=', error
  if (it < self%max_iter .and. io%verbose>2) &
    print '(2(a,i4,2x),2(a,e10.3,2x))', &
     'Poisson: Reached convergence in patch. Iter=', it, &
     'Max iter=', self%max_iter, 'Tolerance=', self%tolerance, 'Error=', error
  !-----------------------------------------------------------------------------
  ! Because we swap pointers each iteration, newest solution is in phi_old
  !-----------------------------------------------------------------------------
  self%phi%f = phi_old
  call trace_end
END SUBROUTINE solve

!> =============================================================================
!> Add the gravitional force to the equation of motion
!> =============================================================================
SUBROUTINE update (self)
  class(ghd_t):: self
  !.............................................................................
  call trace_begin ('ghd_mod:ghd_t%update', 2)
  call self%hd_t%update
  self%dpdt%x = self%dpdt%x - o%ddxdn(self%phi%f)/exp(o%xdn(log(self%d%f)))
  self%dpdt%y = self%dpdt%y - o%ddydn(self%phi%f)/exp(o%ydn(log(self%d%f)))
  self%dpdt%z = self%dpdt%z - o%ddzdn(self%phi%f)/exp(o%zdn(log(self%d%f)))
  call trace_end
END SUBROUTINE update
!> =============================================================================
END MODULE ghd_mod
