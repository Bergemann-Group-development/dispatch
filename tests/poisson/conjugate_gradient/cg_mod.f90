!===============================================================================
!> Conjugate gradient Poisson solver module
!===============================================================================
MODULE cg_mod
  USE lagrange_mod
  USE omp_timer_mod
  USE timer_mod
  implicit none
  private
  integer, parameter:: nt=7
  type, public:: cg_t
    real(8)   :: ds(3), rjac2
    real      :: tolerance, floor, width, courant, omega_sor
    integer   :: max_iter, dims(3), move, order
    integer   :: iit(nt)
    procedure(aslaka), pointer:: solver => null()
    logical   :: precond, red_black
    character(len=16):: method='solve'
  contains
    procedure :: init
    procedure :: laplace
    procedure :: precondition
    procedure :: test
  end type
  logical, save:: dump=.false.
  logical, save:: chebyshev=.false.
  integer, save:: verbose=1
  type(cg_t), public:: cg
CONTAINS

!===============================================================================
!> Conjugate gradient parameters
!===============================================================================
SUBROUTINE init (self)
  class(cg_t)     :: self
  !.............................................................................
  real            :: delta, pref
  real, parameter :: pi = 4.*atan(1.)
  integer         :: i
  real(8)         :: dsi(3), q(3)
  !.............................................................................
  real, save      :: tolerance=1e-6, floor=1e-6, width=0.2, courant=0.2
  real, save      :: omega_sor=0.0
  integer, save   :: max_iter=1000, dims(3)=64, move=8, order=2
  character(len=16), save:: method='solve'
  logical, save   :: precondition=.false., red_black=.true.
  namelist /cg_params/ tolerance, floor, max_iter, width, dims, move, order, &
                       omega_sor, method, precondition, red_black, chebyshev, &
                       verbose, dump
  !-----------------------------------------------------------------------------
  read (*,cg_params)
  self%red_black = precondition     ! for printout
  self%ds        = 1d0/dims
  select case (trim(method))
  case ('jacobi')
    self%solver => jacobi
  case ('ramses')
    self%solver => ramses
  case ('gv1')
    self%solver => gv1
  case ('aslakp')
    self%solver => aslakp
  case ('aslaka')
    self%solver => aslaka
  case default
    self%solver => sor_2nd
    if (omega_sor==0.0) then
      delta = 0.
      pref  = 0.
      do i=1,3
        delta = delta + cos(pi/dims(i))
        pref  = pref + 1.
      end do
      delta = delta/pref
      omega_sor = 2./(1.+sqrt(1.-delta**2))
    end if
    if (chebyshev) then
      dsi = 1d0/self%ds
      q = cos(pi/dims)*dsi**2
      self%rjac2 = (sum(q)/sum(dsi**2))**2
      print *,'dsi, q, rjac2', dsi(1), q(1), self%rjac2
    end if
    self%red_black = red_black
  end select
  write (*,cg_params)
  self%tolerance = tolerance
  self%floor     = floor
  self%max_iter  = max_iter
  self%courant   = courant
  self%width     = width
  self%dims      = dims
  self%move      = move
  self%order     = order
  self%omega_sor = omega_sor
  self%precond   = precondition
  self%method    = method
  !print *, lagrange%weights (3.,[0.,1.,2.])
  !print *, lagrange%weights (4.,[0.,1.,2.,3.])
END SUBROUTINE init

!===============================================================================
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

!===============================================================================
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

!===============================================================================
!> Laplace operator
!===============================================================================
SUBROUTINE laplace (self, phi, out, source)
  class(cg_t)               :: self
  real(8), dimension(:,:,:) :: phi, out
  real, dimension(:,:,:)    :: source
  optional                  :: source
  !.............................................................................
  integer                   :: ix, iy, iz
  real(8)                   :: c(3)
  !-----------------------------------------------------------------------------
  associate (m=>self%dims)
  c = 1d0/self%ds**2
  if (present(source)) then
    do iz=2,m(3)-1
    do iy=2,m(2)-1
    do ix=2,m(1)-1
      out(ix,iy,iz) = c(1)*(phi(ix+1,iy,iz)+phi(ix-1,iy,iz)-2d0*phi(ix,iy,iz)) &
                    + c(2)*(phi(ix,iy+1,iz)+phi(ix,iy-1,iz)-2d0*phi(ix,iy,iz)) &
                    + c(3)*(phi(ix,iy,iz+1)+phi(ix,iy,iz-1)-2d0*phi(ix,iy,iz)) &
                    - source(ix,iy,iz)
    end do
    end do
    end do
  else
    do iz=2,m(3)-1
    do iy=2,m(2)-1
    do ix=2,m(1)-1
      out(ix,iy,iz) = c(1)*(phi(ix+1,iy,iz)+phi(ix-1,iy,iz)-2d0*phi(ix,iy,iz)) &
                    + c(2)*(phi(ix,iy+1,iz)+phi(ix,iy-1,iz)-2d0*phi(ix,iy,iz)) &
                    + c(3)*(phi(ix,iy,iz+1)+phi(ix,iy,iz-1)-2d0*phi(ix,iy,iz))
    end do
    end do
    end do
  end if
  call zero_boundaries8 (out)
  end associate
END SUBROUTINE laplace

!===============================================================================
!> Laplace operator
!===============================================================================
SUBROUTINE precondition (self, res, out)
  class(cg_t)               :: self
  real(8), dimension(:,:,:) :: res, out
  !.............................................................................
  integer                   :: ix, iy, iz
  real(8)                   :: a, b, c
  !-----------------------------------------------------------------------------
  if (self%precond) then
    associate (m=>self%dims)
    a = 1d0 + 1d0/12d0
    b = 1d0/6d0
    if (self%method(1:5)=='aslak') then
      c = 6*b + a
      a = a/c
      b = b/c
    end if
    do iz=2,m(3)-1
    do iy=2,m(2)-1
    do ix=2,m(1)-1
      out(ix,iy,iz) = b*(res(ix+1,iy,iz)+res(ix-1,iy,iz)  &
                    +    res(ix,iy+1,iz)+res(ix,iy-1,iz)  &
                    +    res(ix,iy,iz+1)+res(ix,iy,iz-1)) &
                    + a*res(ix,iy,iz)
    end do
    end do
    end do
    call zero_boundaries8 (out)
    end associate
  else
    out = res
  end if
END SUBROUTINE precondition

!===============================================================================
!> Conjugate gradient Poisson solver, Ghysels & Vanroose 2014, Algorithm 1
!===============================================================================
SUBROUTINE gv1 (self, phi1, source, niter, error)
  class(cg_t)            :: self
  real, dimension(:,:,:) :: phi1, source
  intent(in)             :: source
  integer                :: niter
  real                   :: error
  !.............................................................................
  real(8)                :: pap, ru, ru2
  real(8)                :: alpha, beta
  integer                :: iter, m(3)
  real(8), dimension(:,:,:), allocatable:: phi, r, p, s, u
  !-----------------------------------------------------------------------------
  m = shape(phi1)
  allocate (r(m(1),m(2),m(3)), p(m(1),m(2),m(3)))
  allocate (s(m(1),m(2),m(3)), u(m(1),m(2),m(3)))
  allocate (phi(m(1),m(2),m(3)))
  !-----------------------------------------------------------------------------
  ! Calculate the residual = source - Laplace(phi)
  !-----------------------------------------------------------------------------
  phi = phi1
  call self%laplace (phi, r)
  r = source - r
  call self%precondition (r, u)
  ru = sum(r*u)
  !-----------------------------------------------------------------------------
  ! Main loop
  !-----------------------------------------------------------------------------
  p = u
  do iter = 1,self%max_iter
    call self%laplace (p, s)
    alpha = ru/sum(s*p)
    phi = phi + alpha*p
    call zero_boundaries8 (phi)
    r = r - alpha*s
    call self%precondition (r, u)
    ru2 = sum(r*u)
    error = maxval(abs(r/(source+self%floor)))
    if (error < self%tolerance) exit
    beta = ru2/ru
    p = u + beta*p
    ru = ru2
    if (verbose>2) print *,iter,error,alpha,beta
  end do
  niter = iter
  if (iter==self%max_iter+1) print *, &
    'cg_mod%solve: no convergence, max_iter, error =', self%max_iter, error
  phi1 = phi
  if (verbose>1) then
    call self%laplace (phi, r, source)
    print *, 'min/maxval residual:', minval(r), maxval(r)
  end if
  if (dump) write (11) shape(phi)
  if (dump) write (11) phi, r
END SUBROUTINE gv1

!===============================================================================
!> Conjugate gradient Poisson solver
!===============================================================================
SUBROUTINE aslaka (self, phi1, source, niter, error)
  class(cg_t)            :: self
  real, dimension(:,:,:) :: phi1, source
  intent(in)             :: source
  integer                :: niter
  real                   :: error
  !.............................................................................
  real(8)                :: pap, resnorm, resnorm2
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
    phi = phi + alpha*pes
    call zero_boundaries8 (phi)
    res = res - alpha*Ap
    resnorm2 = sum(res**2)
    error = maxval(abs(res/(source+self%floor)))
    if (error < self%tolerance) exit
    beta = resnorm2/resnorm
    pes = res + beta*pes
    resnorm = resnorm2
    if (verbose>2) print *,iter,error,resnorm,alpha,beta
  end do
  niter = iter
  if (iter==self%max_iter+1) print *, &
    'cg_mod%solve: no convergence, max_iter, error =', self%max_iter, error
  phi1 = phi
  call self%laplace (phi, res, source)
  if (verbose>1) print *, 'min/maxval residual:', minval(res), maxval(res)
  if (dump) write (11) shape(phi)
  if (dump) write (11) phi, res
  deallocate (res, pes, Ap, phi)
END SUBROUTINE aslaka

!===============================================================================
!> Conjugate gradient Poisson solver
!===============================================================================
SUBROUTINE aslakp (self, phi1, source, niter, error)
  class(cg_t)            :: self
  real, dimension(:,:,:) :: phi1, source
  intent(in)             :: source
  integer                :: niter
  real                   :: error
  !.............................................................................
  real(8)                :: pap, resnorm, resnorm2
  real(8)                :: alpha, beta
  integer                :: iter, m(3)
  real(8), dimension(:,:,:), pointer:: phi, res, pes, Ap
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
    phi = phi + alpha*pes
    call zero_boundaries8 (phi)
    res = res - alpha*Ap
    resnorm2 = sum(res**2)
    error = maxval(abs(res/(source+self%floor)))
    if (error < self%tolerance) exit
    beta = resnorm2/resnorm
    pes = res + beta*pes
    resnorm = resnorm2
    if (verbose>2) print *,iter,error,resnorm,alpha,beta
  end do
  niter = iter
  if (iter==self%max_iter+1) print *, &
    'cg_mod%solve: no convergence, max_iter, error =', self%max_iter, error
  call self%laplace (phi, res, source)
  if (verbose>1) print *, 'min/maxval residual:', minval(res), maxval(res)
  phi1 = phi
  if (dump) write (12) shape(phi)
  if (dump) write (12) phi, res
  deallocate (res, pes, Ap, phi)
END SUBROUTINE aslakp

!===============================================================================
!> Conjugate gradient Poisson solver
!===============================================================================
SUBROUTINE ramses (self, phi1, source, niter, error1)
  class(cg_t)            :: self
  real, dimension(:,:,:) :: phi1, source
  intent(in)             :: source
  integer                :: niter
  real                   :: error1, error2
  !.............................................................................
  real(8)                :: gamma, gamma_old, delta
  real(8)                :: alpha, alpha_old, beta
  integer                :: iter, mm(3)
  real(8), dimension(:,:,:), allocatable:: phi, p, q, r, s, w, n, z, res
  real(8), dimension(:,:,:), pointer:: u, m
  !-----------------------------------------------------------------------------
  ! Initially :
  !   r = b - A x
  !   u = M^-1 r
  !   w = A u
  !-----------------------------------------------------------------------------
  mm = shape(phi1)
  allocate (p(mm(1),mm(2),mm(3)), q(mm(1),mm(2),mm(3)))
  allocate (r(mm(1),mm(2),mm(3)), s(mm(1),mm(2),mm(3)))
  allocate (u(mm(1),mm(2),mm(3)), w(mm(1),mm(2),mm(3)))
  allocate (m(mm(1),mm(2),mm(3)), n(mm(1),mm(2),mm(3)))
  allocate (z(mm(1),mm(2),mm(3)), phi(mm(1),mm(2),mm(3)))
  allocate (res(mm(1),mm(2),mm(3)))
  phi = phi1
  call self%laplace (phi, r, source)
  r = -r
  call self%precondition (r, u)
  call self%laplace (u, w)
  call zero_boundaries8 (w)
  !-----------------------------------------------------------------------------
  ! Comments copied verbatim from ramses/poisson/phi_fine_cg.f90:
  ! while (iter<cg_itermax .and. error > error_ini * epsilon)
  !   gamma = (r.u)
  !   delta = (w.u)
  !   if (iter=1) error_ini = error
  !   if (iter>1) then
  !     alpha = gamma / (delta - beta * gamma / alpha_old)
  !     beta  = gamma / gamma_old
  !   else
  !     alpha = gamma / delta
  !     beta = 0
  !   endif
  !   alpha_old = alpha
  !   gamma_old = gamma
  !   m = M^-1 w
  !   n = A m
  !   z = n + beta z
  !   q = m + beta q
  !   s = w + beta s
  !   p = u + beta p
  !   phi = phi + alpha p
  !   r   = r   - alpha s
  !   u   = u   - alpha q
  !   w   = w   - alpha z
  !-----------------------------------------------------------------------------
  do iter=1,self%max_iter
    gamma = sum(r*u)
    delta = sum(w*u)
    if (iter>1) then
      beta  = gamma / gamma_old
      alpha = gamma / (delta - beta * gamma / alpha_old)
    else
      beta = 0
      alpha = gamma / delta
    end if
    alpha_old = alpha
    gamma_old = gamma
    call self%precondition (w, m)
    call zero_boundaries8 (m)
    call self%laplace (m, n)
    call zero_boundaries8 (n)
    z = n + beta*z
    q = m + beta*q
    s = w + beta*s
    p = u + beta*p
    phi = phi + alpha*p
    r   = r   - alpha*s
    u   = u   - alpha*q
    w   = w   - alpha*z
    error1 = maxval(abs(r/(source+self%floor)))
    !call self%laplace (phi, res, source)
    !error2 = maxval(abs(res/(source+self%floor)))
    if (verbose>2) print *, iter, error1, error1
    if (error1 < self%tolerance) exit
  end do
  niter = iter
  if (iter==self%max_iter+1) print *, &
    'cg_mod%ramses: no convergence, max_iter, error =', self%max_iter, error1
  phi1 = phi
  call self%laplace (phi, r, source)
  if (verbose>1) print *, 'min/maxval residual:', minval(r), maxval(r)
  if (dump) write (13) shape(phi)
  if (dump) write (13) phi, r
  deallocate (p,q,r,s,u,w,m,n,z,phi)
END SUBROUTINE ramses

!===============================================================================
!> Iterative Poisson solver
!===============================================================================
SUBROUTINE jacobi (self, phi1, source, niter, error)
  class(cg_t)            :: self
  real, dimension(:,:,:) :: phi1, source
  intent(in)             :: source
  integer                :: niter
  real                   :: error, alpha
  !.............................................................................
  integer                :: iter
  real(8), dimension(:,:,:), allocatable:: phi, res
  !-----------------------------------------------------------------------------
  associate (m=>self%dims)
  allocate (phi(m(1),m(2),m(3)), res(m(1),m(2),m(3)))
  alpha = self%ds(1)**2/12.
  phi = phi1
  do iter=1,self%max_iter
    call self%laplace (phi, res, source)
    error = maxval(abs(res/(source+self%floor)))
    if (error < self%tolerance) exit
    phi = phi + alpha*res
  end do
  niter = iter
  phi1 = phi
  deallocate (phi, res)
  end associate
END SUBROUTINE jacobi

!===============================================================================
!> 3pt 2nd order Red-Black Gauss-Seidel iteration with SOR acceleration
!===============================================================================
SUBROUTINE sor_2nd (self, phi1, source, niter, error)
  class(cg_t)            :: self
  real, dimension(:,:,:) :: phi1, source
  intent(in)             :: source
  integer                :: niter
  real                   :: error
  !.............................................................................
  integer :: i, j, k, sweep, ss, iter, i1
  real    :: a, b, b1, b2, b3, sor
  real(8), dimension(:,:,:), allocatable:: phi, res
  associate (n=>self%dims, ds=>self%ds, l=>[2,2,2], u=>self%dims-1)
  !-----------------------------------------------------------------------------
  ! Make one iterative update for
  !
  !   (2/dx^2 + 2/dy^2 + 2/dz^2) phi_new_ijk =
  !                                 (phi_i+1,j,k + phi_i-1,j,k)/dx^2
  !                               + (phi_i,j+1,k + phi_i,j-1,k)/dy^2
  !                               + (phi_i,j,k+1 + phi_i,j,k-1)/dz^2
  !                               - source
  !.............................................................................
  ! Calculate 1 / [(2/dx + 2/dy + 2/dz)] prefactor according to dimensionality
  allocate (phi(n(1),n(2),n(3)), res(n(1),n(2),n(3)))
  phi = phi1
  a = 0.0
  if (n(1)>1) a = a + 2. / ds(1)**2
  if (n(2)>1) a = a + 2. / ds(2)**2
  if (n(3)>1) a = a + 2. / ds(3)**2
  a = 1. / a
  !print *,a,n,ds
  if (chebyshev) then
    self%omega_sor = 1d0/(1d0-0.5d0*self%rjac2)
  end if
  do iter=1,self%max_iter
    error = 0.0
    do sweep=0,merge(1,0,self%red_black)
      ! Add up contributions
      b1 = a / ds(1)**2
      b2 = a / ds(2)**2
      b3 = a / ds(3)**2
      b = 2.0*(b1+b2+b3)
      do k=l(3),u(3)
        do j=l(2),u(2)
          ss = modulo(j+k+sweep,2)
          ss = merge(ss,0,self%red_black)
          i1 = merge(2,1,self%red_black)
          do i=l(1)+ss,u(1),i1
            phi(i,j,k)     = (1. - self%omega_sor) * phi(i,j,k) + &
                                   self%omega_sor * ( &
                                     b1 * (phi(i+1,j,k)+phi(i-1,j,k)) &
                                   + b2 * (phi(i,j+1,k)+phi(i,j-1,k)) &
                                   + b3 * (phi(i,j,k+1)+phi(i,j,k-1)) &
                                   - a * source(i,j,k) )
          end do
          !-----------------------------------------------------------------------
          ! At least with gfortran, keeping these loops separate is faster
          !-----------------------------------------------------------------------
          do i=l(1)+ss,u(1),2
            error = max(error, abs( &
              ((b1*(phi(i+1,j,k)+phi(i-1,j,k)) + &
                b2*(phi(i,j+1,k)+phi(i,j-1,k)) + &
                b3*(phi(i,j,k+1)+phi(i,j,k-1)) - &
                b * phi(i,j,k)) - a*source(i,j,k)))/(a*(source(i,j,k)+self%floor)))
            block
            real a1,a2,a3
            a1 = b1*(phi(i+1,j,k)+phi(i-1,j,k)) + &
                 b2*(phi(i,j+1,k)+phi(i,j-1,k)) + &
                 b3*(phi(i,j,k+1)+phi(i,j,k-1))
            a2 = b*phi(i,j,k)
            a3 = a*source(i,j,k)
            end block
          end do
        end do
      end do
      if (chebyshev) then
        self%omega_sor = 1d0/(1d0-self%rjac2*0.25d0*self%omega_sor)
      end if
    end do
    if (error < self%tolerance) exit
    if (verbose>2) print *, iter, error, self%omega_sor, self%rjac2
  end do
  call self%laplace (phi, res, source)
  if (verbose>1) print *, 'min/maxval residual:', minval(res), maxval(res)
  if (dump) write (14) shape(phi)
  if (dump) write (14) phi, res
  niter = iter
  phi1 = phi
  end associate
END SUBROUTINE sor_2nd

!===============================================================================
!> Setup a test problem for the Poisson equation
!===============================================================================
SUBROUTINE test (self)
  class(cg_t)                       :: self
  !.............................................................................
  real, dimension(:,:,:,:), pointer :: phii
  real, dimension(:,:,:), pointer   :: source, phi
  real    :: x, y, z, ds(3), c, r2, dx, error
  integer :: ix, iy, iz, i0, i, it, jt, new, niter, order, sum_iter
  real(8) :: wc
  character(len=16):: option
  !-----------------------------------------------------------------------------
  ds = 1d0/self%dims
  c = 0.5d0/self%width**2
  associate (m=>self%dims)
  allocate (phii(m(1),m(2),m(3),nt))
  allocate (source(m(1),m(2),m(3)))
  phii = 0.0
  it = 1
  jt = 1
  sum_iter = 0
  self%iit = 1
  i0 = self%move/self%courant
  wc = wallclock()
  call timer%tic
  do i=-i0,i0
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
    phi => phii(:,:,:,it)
    call self%solver (phi, source, niter, error)
    sum_iter = sum_iter+niter
    new = 1 + mod(it,nt)
    self%iit(1:nt-1) = self%iit(2:nt)
    self%iit(nt) = new
    order = min(self%order,jt-1)
    if (order==-1) then
      phii(:,:,:,new) = 0.0
    else if (order==0) then
      phii(:,:,:,new) = phii(:,:,:,it)
    else if (order==1) then
      phii(:,:,:,self%iit(nt)) &
      = phii(:,:,:,self%iit(nt-1))*2.0 &
      - phii(:,:,:,self%iit(nt-2))
    else if (order==2) then
      phii(:,:,:,self%iit(nt)) &
      = phii(:,:,:,self%iit(nt-1))*3.0 &
      - phii(:,:,:,self%iit(nt-2))*3.0 &
      + phii(:,:,:,self%iit(nt-3))
    else
      phii(:,:,:,self%iit(nt)) &
      = phii(:,:,:,self%iit(nt-1))*4.0 &
      - phii(:,:,:,self%iit(nt-2))*6.0 &
      + phii(:,:,:,self%iit(nt-3))*4.0 &
      - phii(:,:,:,self%iit(nt-4))
    end if
    print '(a,f8.2,i5,1p,e10.2,2i3,2x,7i3)', &
      'position, niter, error =', dx, niter, error, it, new, self%iit
    it = new
    jt = jt+1
  end do
  print *, 'min/maxval(phi) =', minval(phi), maxval(phi)
  print *, 'total number of iterations =', sum_iter
  call timer%toc('Poisson',jt*product(self%dims))
  if (trim(self%method)=='sor') then
    option = '   red_black:'
  else
    option = 'precondition:'
  end if
  print '(a,a6,4x,a13,l3,4x,a,i3,4x,a,f7.2)', &
    'method: ', self%method, &
    option, self%red_black, &
    'order:', order, &
    'TIME:', wallclock()-wc
  end associate
END SUBROUTINE test

END MODULE
