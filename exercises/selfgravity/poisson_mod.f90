!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the Poisson problem in 3-D
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE pde_mod
  USE io_mod
  USE omp_timer_mod
  USE trace_mod
  USE patches_mod
  USE scalar_mod
  USE laplace_mod
  implicit none
  private
  integer, parameter:: i_d=1, i_phi=2
  type, public, extends(patch_t):: pde_t
    type(scalar_t):: d, phi
    real:: gamma=5./3.
    real:: w_gauss
    real:: p_gauss
    real:: amplitude
    real:: G                                                                    ! Strength of gravity
    real:: tolerance
    integer:: min_iter, max_iter
    type(laplace_t)  :: laplace
  contains
    procedure:: init
    procedure:: solve
  end type
  real    :: tolerance=1e-4, G=1., omega_sor=0.0
  integer :: min_iter=0, max_iter=100
  character(len=12):: solver='jacobi'
CONTAINS

!> =============================================================================
!> Setup a patch with density and gravitational potential for the Poisson solver
!> =============================================================================
SUBROUTINE init (self)
  class(pde_t):: self
  integer:: i, j, k
  real:: x, y, z
  logical:: first_time = .true.
  real:: w_gauss=5.0, p_gauss=1.0, amplitude=1.0
  namelist /poisson_params/ w_gauss, p_gauss, amplitude, G, tolerance, &
                            min_iter, max_iter, solver, omega_sor
  !-----------------------------------------------------------------------------
  ! Read standard parameters and copy to structure
  !-----------------------------------------------------------------------------
  experimentname = 'poisson'
  call io%init (experimentname)
  call trace_begin ('hd_mod:pde_t%init', 2)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, poisson_params)
    if (io%master) write (*,poisson_params)
  end if
  self%size      = 1.
  self%position  = 0.
  self%w_gauss   = w_gauss
  self%p_gauss   = p_gauss
  self%amplitude = amplitude
  self%G         = G
  self%tolerance = tolerance
  self%min_iter  = min_iter
  self%max_iter  = max_iter
  self%nv        = 2
  call self%patch_t%init
  !-----------------------------------------------------------------------------
  ! Set up the right dimensions from current patch, etc for laplace iteration method
  !-----------------------------------------------------------------------------
  call self%laplace%init(self,solver,omega_sor)
  !-----------------------------------------------------------------------------
  ! Memory mapping og density and potentiall field
  !-----------------------------------------------------------------------------
  call    self%d%init (self%gn, self%size/self%n, 'd')
  call  self%phi%init (self%gn, self%size/self%n, 'phi')
  if (self%d%allocated)   deallocate(self%d%f)
  if (self%phi%allocated) deallocate(self%phi%f)
  self%d%f   => self%mem(:,:,:,i_d   ,self%it,1)
  self%phi%f => self%mem(:,:,:,i_phi ,self%it,1)
  !-----------------------------------------------------------------------------
  ! Physical initial conditions: gaussian top hat
  !-----------------------------------------------------------------------------
  if (self%time==0.0) then
    do k=1,self%gn(3)
    do j=1,self%gn(2)
    do i=1,self%gn(1)
      x = (i-(self%gn(1)+1)*0.5)*self%ds(1) + self%position(1)
      y = (j-(self%gn(2)+1)*0.5)*self%ds(2) + self%position(2)
      z = (k-(self%gn(3)+1)*0.5)*self%ds(3) + self%position(3)
      self%d%f(i,j,k) = &
        exp(self%amplitude*exp(-(sqrt(x**2+y**2+z**2)/self%w_gauss)**self%p_gauss))
    end do
    end do
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Initial and boundary conditions for phi
  !-----------------------------------------------------------------------------
  self%phi%f = 0.0
  !
  call trace_end
END SUBROUTINE init
!> =============================================================================
!> Compute iteratively the solution to the Poisson equation stoppin when reaching the tolerance tol
!> =============================================================================
SUBROUTINE solve (self)
  class(pde_t):: self
  !
  real, allocatable, dimension(:,:,:)         :: FourPiGrho
  real, allocatable, dimension(:,:,:), target :: phi_scratch
  real, pointer,     dimension(:,:,:)         :: phi_old, phi_new, phi_swap
  real, parameter :: FourPi = 4.0 * 3.14159265358979323846
  real            :: error
  integer         :: it
  !.............................................................................
  call trace_begin ('poisson_mod:pde_t%update', 2)

associate(gn=>self%gn)
  Allocate(FourPiGrho(gn(1),gn(2),gn(3)),phi_scratch(gn(1),gn(2),gn(3)))        ! Allocate scratch memory
end associate
  FourPiGrho = FourPi * self%G * self%d%f                                       ! Set up source of gravity

  phi_old => self%phi%f                                                         ! Do initial pointer arrangement
  phi_new => phi_scratch

  error = 2.*self%tolerance
  it    = 0
  call tic
  do while((it<self%min_iter) .or. &                                            ! Do at least min_iter iterations 
           (error >= self%tolerance .and. it<self%max_iter))                    ! Stop if error<=tolerance or max iterations reached

    call self%laplace%solver (phi_old, phi_new, FourPiGrho, error)              ! do one iteration and get error
    if (io%verbose>2) write(*,'(a,i4,2(a,e10.3))'), 'Poisson: Iter=', it, ' Tolerance=', self%tolerance, ' Error=', error
    phi_swap => phi_new; phi_new=>phi_old; phi_old=>phi_swap                    ! swap phi_old and phi_new pointers
    it = it + 1

  enddo
  call toc ('poisson', product(self%n))

  ! Check if we reached solution in time
  if (it == self%max_iter .and. io%verbose>0) &
    write(*,'(a,i4,2(a,e10.3))'), 'Poisson: Did not reach convergence in patch. Max iter=', it, &
                                  ' Tolerance=', self%tolerance, ' Error=', error
  if (it < self%max_iter .and. io%verbose>1) &
    write(*,'(2(a,i4),2(a,e10.3))'), 'Poisson: Reached convergence in patch. Iter=', it, &
                                  ' Max iter=', self%max_iter, ' Tolerance=', self%tolerance, ' Error=', error

  ! Because we swap pointers at the end of each iteration, newest solution is in phi_old
  ! Check if we have to copy over solution
  if (associated(phi_old,phi_scratch)) self%phi%f = phi_scratch
  
  call trace_end
END SUBROUTINE solve
!> =============================================================================
END MODULE pde_mod
