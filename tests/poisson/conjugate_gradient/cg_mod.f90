!> =============================================================================
!> Conjugate gradient Poisson solver module
!> =============================================================================
MODULE cg_mod
  implicit none
  private
  type, public:: cg_t
    real(8)   :: ds(3)
    real      :: tolerance, floor, width, courant
    integer   :: max_iter, dims(3), move, nt
  contains
    procedure :: init
    procedure :: laplace
    procedure :: solve
    procedure :: test
  end type
  integer, parameter:: verbose=1
  type(cg_t), public:: cg
CONTAINS

!> =============================================================================
!> Conjugate gradient parameters
!> =============================================================================
SUBROUTINE init (self)
  class(cg_t)   :: self
  !.............................................................................
  real, save    :: tolerance=1e-6, floor=1e-6, width=0.2, courant=0.2
  integer, save :: max_iter=1000, dims(3)=64, move=8, nt=5
  !.............................................................................
  namelist /cg_params/ tolerance, floor, max_iter, width, dims
  !-----------------------------------------------------------------------------
  !read (*,cg_params)
  self%tolerance = tolerance
  self%floor     = floor
  self%max_iter  = max_iter
  self%courant   = courant
  self%width     = width
  self%dims      = dims
  self%move      = move
  self%nt        = nt
  self%ds        = 1d0/dims
  write (*,cg_params)
END SUBROUTINE init

!> =============================================================================
SUBROUTINE zero_boundaries (f)
  real:: f(:,:,:)
  !.............................................................................
  f(1,:,:) = 0.0
  f(:,1,:) = 0.0
  f(:,:,1) = 0.0
  f(size(f,1),:,:) = 0.0
  f(:,size(f,2),:) = 0.0
  f(:,:,size(f,3)) = 0.0
END SUBROUTINE zero_boundaries

!> =============================================================================
SUBROUTINE zero_boundaries8 (f)
  real(8):: f(:,:,:)
  !.............................................................................
  f(1,:,:) = 0.0
  f(:,1,:) = 0.0
  f(:,:,1) = 0.0
  f(size(f,1),:,:) = 0.0
  f(:,size(f,2),:) = 0.0
  f(:,:,size(f,3)) = 0.0
END SUBROUTINE zero_boundaries8

!> =============================================================================
!> Laplace operator
!> =============================================================================
SUBROUTINE laplace (self, phi, out)
  class(cg_t)               :: self
  real(8), dimension(:,:,:) :: phi, out
  !.............................................................................
  integer                   :: ix, iy, iz
  real(8)                   :: c(3)
  !-----------------------------------------------------------------------------
  associate (m=>self%dims)
  c = 1d0/self%ds**2
  do iz=2,m(3)-1
  do iy=2,m(2)-1
  do ix=2,m(1)-1
    out(ix,iy,iz) = c(1)*(phi(ix+1,iy,iz)+phi(ix-1,iy,iz)-2d0*phi(ix,iy,iz)) &
                  + c(2)*(phi(ix,iy+1,iz)+phi(ix,iy-1,iz)-2d0*phi(ix,iy,iz)) &
                  + c(3)*(phi(ix,iy,iz+1)+phi(ix,iy,iz-1)-2d0*phi(ix,iy,iz))
  end do
  end do
  end do
  call zero_boundaries8 (out)
  end associate
END SUBROUTINE laplace

!> =============================================================================
!> Conjugate gradient Poisson solver
!> =============================================================================
SUBROUTINE solve (self, phi1, source, niter, error)
  class(cg_t)            :: self
  real, dimension(:,:,:) :: phi1, source
  integer                :: niter
  real                   :: error
  !.............................................................................
  real(8)                :: pap, resnorm, resnorm2, phi_tol
  real(8)                :: alpha, beta
  integer                :: iter, m(3)
  real(8), dimension(:,:,:), allocatable:: phi, res, pes, Ap
  !-----------------------------------------------------------------------------
  m = shape(phi1)
  allocate (res(m(1),m(2),m(3)), pes(m(1),m(2),m(3)))
  allocate ( Ap(m(1),m(2),m(3)), phi(m(1),m(2),m(3)))
  !-----------------------------------------------------------------------------
  ! Calculate the residual = source - Laplace(phi)
  !-----------------------------------------------------------------------------
  phi = phi1
  call self%laplace (phi, res)
  res = source - res
  resnorm = sum(res**2)
  !-----------------------------------------------------------------------------
  ! Main loop
  !-----------------------------------------------------------------------------
  pes = res
  do iter = 1,self%max_iter
    call self%laplace (pes, Ap)
    pap = sum(pes*Ap)
    alpha = resnorm/pap
    !---------------------------------------------------------------------------
    ! Update phi and res
    !---------------------------------------------------------------------------
    phi = phi + alpha*pes
    call zero_boundaries8 (phi)
    res = res - alpha*Ap
    resnorm2 = sum(res**2)
    error = maxval(abs(res/(source+self%floor)))
    if (error < self%tolerance) exit
    beta = resnorm2/resnorm
    pes = res + beta*pes
    resnorm = resnorm2
    if (verbose>1) print *,iter,error,resnorm,alpha,beta
  end do
  niter = iter
  if (iter==self%max_iter+1) print *, &
    'cg_mod%solve: no convergence, max_iter, error =', self%max_iter, error
  phi1 = phi
  deallocate (res, pes, Ap, phi)
END SUBROUTINE solve

!> =============================================================================
!> Setup a test problem for the Poisson equation
!> =============================================================================
SUBROUTINE test (self)
  class(cg_t)                       :: self
  real, dimension(:,:,:,:), pointer :: phii
  real, dimension(:,:,:), pointer   :: source, phi
  real                              :: x, y, z, ds(3), c, r2, dx, error
  integer                           :: ix, iy, iz, i0, i, it, new, niter
  !-----------------------------------------------------------------------------
  ds = 1d0/self%dims
  c = 0.5d0/self%width**2
  associate (m=>self%dims, nt=>self%nt)
  allocate (phii(m(1),m(2),m(3),nt))
  allocate (source(m(1),m(2),m(3)))
  phii = 0.0
  it = 1
  i0 = self%move/self%courant
  do i=-i0,i0
    phi => phii(:,:,:,it)
    dx = i*self%courant
    do iz=1,m(3)
      z = (iz-(m(3)+1d0)/2d0)*ds(3)
      do iy=1,m(2)
        y = (iy-(m(2)+1d0)/2d0)*ds(2)
        do ix=1,m(1)
          x = (ix-dx-(m(1)+1d0)/2d0)*ds(1)
          r2 = x**2 + y**2 + z**2
          source(ix,iy,iz) = exp(-r2*c)
        end do
      end do
    end do
    call zero_boundaries (source)
    call self%solve (phi, source, niter, error)
    print *, 'position, niter, error =', dx, niter, error
    new = 1 + mod(it,nt)
    phii(:,:,:,new) = phii(:,:,:,it)             ! 0th order extrapolation
    it = new
  end do
  end associate
END SUBROUTINE test

END MODULE
