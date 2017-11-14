!> =============================================================================
!> Methods for solving the Poisson equation; derived from routines developed by
!> Troels Haugbølle for the course Computational Astrophysics
!> =============================================================================
MODULE poisson_mod
  USE mpi_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE patch_mod
  USE mesh_mod
  USE bits_mod
  implicit none
  private
  type, public:: poisson_t
    character(len=32):: kind
    real:: omega_sor
    procedure(sor_2nd), pointer :: update => null()
    ! for use with ZEUS SOR solver
    real(8):: rjac2 ! spectral radius
    real(kind=KindScalarVar), dimension(:), allocatable:: dx1m2i, dx2m2i, dx1p2i, dx2p2i
    real(kind=KindScalarVar), dimension(:,:,:), allocatable:: eta
  contains
    procedure:: init
    procedure:: cg
    procedure:: sor_iter
    procedure, nopass:: laplace
    procedure:: check => check_poisson
  end type
  real(8), parameter:: pi = asin(1.0) * 2.0
  real,    save:: floor=1e-2, tolerance=1e-4
  integer, save:: verbose=0
  integer, save:: max_iter=500
CONTAINS

!> =============================================================================
!> Initialize poisson_params
!> =============================================================================
SUBROUTINE init (self, patch, selfgravity)
  implicit none
  class(poisson_t)  :: self
  class(patch_t)    :: patch
  logical           :: selfgravity
  !.............................................................................
  integer           :: iostat
  !.............................................................................
  character(len=12), save :: solver='cg'
  real,              save :: omega_sor=0.0
  logical,           save :: first_time = .true.
  namelist /poisson_params/ solver, omega_sor, tolerance, floor, max_iter, verbose
  !-----------------------------------------------------------------------------
  !$omp critical (poisson_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, poisson_params, iostat=iostat)
    if (io%master) write (io%output, poisson_params)
  end if
  !$omp end critical (poisson_cr)
  select case(trim(solver))
  case('sor_2nd')
    self%omega_sor = omega_sor
    call init_sor_2nd(self, patch)
    self%update => sor_2nd
  case('sor_iter')
    self%omega_sor = omega_sor
    call init_sor_2nd(self, patch)
  case('zeus_sor')
    call init_zeus_sor(self, patch)
    self%update => zeus_sor
  case('cg')
  case default
    if (io%master) then
      write(*,*) 'Poisson solver '//trim(solver)//' is unknown'
      write(*,*) 'Valid solvers are: cg, sor_2nd, zeus_sor'
    end if
    selfgravity = .false.
  end select
  self%kind = solver
END SUBROUTINE init

!> =============================================================================
!> Initialize SOR factor
!> =============================================================================
SUBROUTINE init_sor_2nd (self, patch)
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
  logical, save     :: first_time = .true.
  !-----------------------------------------------------------------------------
  ! Note that, when not using red-black ordering, cos(pi/n) -> cos(2*pi/n)
  !-----------------------------------------------------------------------------
  if (self%omega_sor==0.) then
    delta = 0.0; pref = 0.0
    do i=1,3
      if (patch%n(i) > 1) then
        pref  = pref + 1.
        delta = delta + cos(pi/patch%n(i))
      endif
    enddo
    delta = delta / pref
    self%omega_sor = 2.0 / (1. + sqrt(1.0 - delta**2))
    comment =' (computed)'
  else
    comment =' (from namelist)'
  endif
  if (self%omega_sor < 0. .or. self%omega_sor>2) call mpi%abort( &
    'omega_sor in Laplace_mod has to be between 0 and 2 for stability')
  !$omp critical (poisson_cr2)
  if (first_time .and. io%master) then
    print '(a,f6.4,a)', 'poisson_mod: SOR prefactor set to: ', self%omega_sor, comment
    first_time = .false.
  endif
  !$omp end critical (poisson_cr2)
END SUBROUTINE init_sor_2nd

!> =============================================================================
!> 3pt 2nd order Red-Black Gauss-Seidel iteration with SOR acceleration
!>
!> Make one iterative update for
!>
!>   (2/dx^2 + 2/dy^2 + 2/dz^2) phi_new_ijk =
!>                                 (phi_i+1,j,k + phi_i-1,j,k)/dx^2
!>                               + (phi_i,j+1,k + phi_i,j-1,k)/dy^2
!>                               + (phi_i,j,k+1 + phi_i,j,k-1)/dz^2
!>                               - FourPiGrho
!>
!> Using double precision for "res" has no noticeable effect, but using double
!> for the b[123] coefficients does have a positive effect.
!> =============================================================================
SUBROUTINE sor_2nd (self, patch, phi, resid, FourPiGrho, d, error, detailed_timer)
  class(poisson_t) :: self
  class(patch_t)   :: patch
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: phi, resid
  real(kind=KindScalarVar), dimension(:,:,:), pointer, intent(in):: FourPiGrho, d
  logical:: detailed_timer
  integer :: i, j, k, sweep, ss, ione, jone, kone
  real    :: error, res
  real(8) :: a, b, b1, b2, b3
  integer, save          :: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('poisson_t%sor_2nd_abs', itimer=itimer)
  associate(n=>patch%n,ds=>patch%ds,l=>patch%li,u=>patch%ui,gn=>patch%gn,omega_sor=>self%omega_sor)

  ione = 1
  jone = 1
  kone = 1
  if (n(1) <= 1) ione = 0
  if (n(2) <= 1) jone = 0
  if (n(3) <= 1) kone = 0

  ! Calculate 1 / [(2/dx + 2/dy + 2/dz)] prefactor according to dimensionality
  a = 0.0
  if (n(1) > 1) a = a + 2. / ds(1)**2
  if (n(2) > 1) a = a + 2. / ds(2)**2
  if (n(3) > 1) a = a + 2. / ds(3)**2
  a = 1. / a

  error = 0.0
  !do sweep=0,1
    ! Add up contributions
    b1 = a / ds(1)**2 * ione
    b2 = a / ds(2)**2 * jone
    b3 = a / ds(3)**2 * kone
    do k=l(3),u(3)
      do j=l(2),u(2)
        !ss = modulo(j + k + sweep,2)
        !do i=l(1)+ss,u(1),2
        do i=l(1),u(1)
          !---------------------------------------------------------------------
          ! This form of the residual has no factor in front of phi, which makes
          ! it clear how to apply the over-relaxation factor
          !---------------------------------------------------------------------
          res = b1 * (phi(i+ione,j,k) + phi(i-ione,j,k)) &
              + b2 * (phi(i,j+jone,k) + phi(i,j-jone,k)) &
              + b3 * (phi(i,j,k+kone) + phi(i,j,k-kone)) &
              - phi(i,j,k)- a * FourPiGrho(i,j,k)
          phi(i,j,k) = phi(i,j,k) + self%omega_sor * res
          resid(i,j,k) = res/d(i,j,k)
          error = max(abs(resid(i,j,k)), error)
        end do
      end do
    end do
  !end do
  !-----------------------------------------------------------------------------
  ! We are supposed to return the residual as (Laplace(phi)-source)/d, so we
  ! need to divide with the factor that source (FourPiGrho) was multiplied by
  !-----------------------------------------------------------------------------
  error = error / a
  end associate

  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE sor_2nd

SUBROUTINE sor_iter (self, patch, phi1, source, error0, error, niter, detailed_timer)
  class(poisson_t)       :: self
  class(patch_t)         :: patch
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: phi1
  real(kind=KindScalarVar), dimension(:,:,:), pointer, intent(in):: source
  real                   :: error, error0
  integer                :: niter
  logical                :: detailed_timer
  !.............................................................................
  integer                :: i, j, k, sweep, ss, ione, jone, kone, m(3), iter
  real(8)                :: a, b, b1, b2, b3, res
  integer, save          :: itimer=0
  real(8), allocatable   :: phi(:,:,:)
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('poisson_t%sor_iter', itimer=itimer)
  associate(n=>patch%n,ds=>patch%ds,l=>patch%li,u=>patch%ui,gn=>patch%gn)
  allocate (phi(gn(1),gn(2),gn(3)))
  phi = phi1
  ione = 1
  jone = 1
  kone = 1
  if (n(1) <= 1) ione = 0
  if (n(2) <= 1) jone = 0
  if (n(3) <= 1) kone = 0
  ! Calculate 1 / [(2/dx + 2/dy + 2/dz)] prefactor according to dimensionality
  a = 0.0
  if (n(1) > 1) a = a + 2. / ds(1)**2
  if (n(2) > 1) a = a + 2. / ds(2)**2
  if (n(3) > 1) a = a + 2. / ds(3)**2
  a = 1. / a
  do iter=1,max_iter
    error = 0.0
    !do sweep=0,1
      ! Add up contributions
      b1 = a / ds(1)**2 * ione
      b2 = a / ds(2)**2 * jone
      b3 = a / ds(3)**2 * kone
      do k=l(3),u(3)
        do j=l(2),u(2)
          !ss = modulo(j + k + sweep,2)
          !do i=l(1)+ss,u(1),2
          do i=l(1),u(1)
            !---------------------------------------------------------------------
            ! This form of the residual has no factor in front of phi, which makes
            ! it clear how to apply the over-relaxation factor
            !---------------------------------------------------------------------
            res = b1 * (phi(i+ione,j,k) + phi(i-ione,j,k)) &
                + b2 * (phi(i,j+jone,k) + phi(i,j-jone,k)) &
                + b3 * (phi(i,j,k+kone) + phi(i,j,k-kone)) &
                - phi(i,j,k)- a * source(i,j,k)
            phi(i,j,k) = phi(i,j,k) + self%omega_sor * res
            error = max(abs(res), error)
          end do
        end do
      end do
    !end do
    !-----------------------------------------------------------------------------
    ! We are supposed to return the residual as (Laplace(phi)-source)/d, so we
    ! need to divide with the factor that source (FourPiGrho) was multiplied by
    !-----------------------------------------------------------------------------
    error = error / a
    if (iter==1) error0 = error
    if (verbose > 1) print *, 'poisson_t%sor_iter: iter, error =', iter, error
    if (error < tolerance) exit
  end do
  niter = iter
  phi1 = phi
  deallocate (phi)
  end associate
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE sor_iter

!> =====================================================================
!> Initialise variables for the ZEUS SOR solver which are independent of
!> the MHD update and iterative calls to the SOR solver.
!> See accompanying `zeus_sor` procedure.
!> =====================================================================
SUBROUTINE init_zeus_sor (self, patch)
  class(poisson_t):: self
  class(patch_t):: patch
  class(mesh_t), pointer:: m1, m2, m3
  real(8):: q1, q2, q3, etaik, etajk, etakk, dx1a2i, dx2a2i, dx3a2i, dsi(3)
  integer:: i, j, k
  !---------------------------------------------------------------------
  call trace%begin('poisson_t%init_zeus_sor')

  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)

  dsi = 1.0 / patch%ds

  if (allocated(self%dx1m2i)) then
    deallocate (self%dx1m2i, self%dx1p2i, self%dx2m2i, self%dx2p2i, self%eta)
  end if
  allocate (self%dx1m2i(m1%gn), self%dx1p2i(m1%gn), self%dx2m2i(m2%gn), self%dx2p2i(m2%gn))
  allocate (self%eta(m1%gn,m2%gn,m3%gn))

  ! Compute grid-related quantities that don't change from one cycle to the next.
  ! These are required for curvilinear coords and go 1.0/ds**2 in Cartesian coords.
  ! Uniform cell sizes are assumed.
  self%dx1m2i(:) = 0.0
  self%dx1p2i(:) = 0.0
  self%dx2m2i(:) = 0.0
  self%dx2p2i(:) = 0.0
  do i=m1%li,m1%ui
    self%dx1m2i(i) = m1%h2f(i  ) * m1%h31f(i  ) / m1%h2c(i  ) / m1%h31c(i  ) * dsi(1)**2
    self%dx1p2i(i) = m1%h2f(i+1) * m1%h31f(i+1) / m1%h2c(i  ) / m1%h31c(i  ) * dsi(1)**2
  end do
  do j=m2%li,m2%ui
    self%dx2m2i(j) = m2%h32f(j  ) / m2%h32c(j) * dsi(2)**2
    self%dx2p2i(j) = m2%h32f(j+1) / m2%h32c(j) * dsi(2)**2
  end do

  !  Evaluate quantities needed by the SOR iteration loop, but are independent
  !  of phi.  These need only be calculated once for each patch. The Chebyshev
  !  acceleration factor, `omega_sor`, is also initialised.
  !
  !  Note from D. Clarke: omega_sor is started off at omega(n+1/2) in the series given
  !  by Press et al (page 859), and not omega(0) = 1.  This seems to allow
  !  for more rapid convergence.
  etakk  = 0.0
  etajk  = 0.0
  etaik  = 0.0
  q1     = 0.0
  q2     = 0.0
  q3     = 0.0
  dx1a2i = 0.0
  dx2a2i = 0.0
  dx3a2i = 0.0

  if (patch%n(3) > 1) then
    dx3a2i = dsi(3)**2
    etakk  =-2.0 * dx3a2i
    q3     = cos(pi / patch%n(3)) * dx3a2i
  end if
  if (patch%n(2) > 1) then
    dx2a2i = dsi(2)**2
    q2 = cos(pi / patch%n(2)) * dx2a2i
  end if
  if (patch%n(1) > 1) then
    dx1a2i = dsi(1)**2
    q1 = cos(pi / patch%n(1)) * dx1a2i
  end if
  self%rjac2 = ( ( q1 + q2 + q3 ) / ( dx1a2i + dx2a2i + dx3a2i ) )**2 ! Jacobi spectral radius squared
  do k=m3%li,m3%ui
    do j=m2%li,m2%ui
      if (patch%n(2) > 1) etajk =-self%dx2m2i(j) - self%dx2p2i(j)
      do i=m1%li,m1%ui
        if (patch%n(1) > 1) etaik =-self%dx1m2i(i) - self%dx1p2i(i)
        self%eta(i,j,k) = etaik &
                        + etajk / m1%h2c (i)**2 &
                        + etakk / (m1%h31c(i) * m2%h32c(j))**2
      end do
    end do
  end do

  call trace%end
END SUBROUTINE init_zeus_sor

!> =====================================================================
!> ZEUS-3D's version of SOR for solving Poisson's eqn.
!> Extended for curvilinear coords.
!> Originally written by Chris Howard in 2003 and modified by Alexander
!> Men'shchikov (2003), David Clarke (2004) and Nick MacDonald (2005).
!> Ported to DISPATCH by JPR (2017).
!      The SOR algorithm from Numerical Recipes (Press et al.), has been
  !  generalised for curvilinear geometries (XYZ, ZRP, RTP) and with any
  !  symmetry (1-D, 2-D, or full 3-D).
!> =====================================================================
SUBROUTINE zeus_sor (self, patch, phi, resid, FourPiGrho, d, anorm, detailed_timer)
  class(poisson_t):: self
  class(patch_t):: patch
  class(mesh_t), pointer:: m1, m2, m3
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: phi, resid
  real(kind=KindScalarVar), dimension(:,:,:), pointer, intent(in):: FourPiGrho, d
  real   :: anorm
  logical:: detailed_timer
  integer, save:: itimer = 0
  integer:: ione, jone, kone, i, j, k, sweep, sw1, sw2, sw3 &
          , inci, incj, inck, swi, swj, swk
  real(8):: da, dx3a2i
  !-----------------------------------------------------------------------
  if (detailed_timer) call trace%begin('poisson_mod%zeus_sor', itimer=itimer)
  associate (l => patch%li, u => patch%ui)

  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)

  ione = 1
  jone = 1
  kone = 1
  if (patch%n(1) <= 1) ione = 0
  if (patch%n(2) <= 1) jone = 0
  if (patch%n(3) <= 1) kone = 0
  ! loop strides
  inci = 1 + ione
  inck = 1 + (1-ione) * (1-jone)
  incj = 3 - inci*inck

  !  The checkerboarding described in Press et al. is performed on the
  !  innermost do-loop only, which is normally the i-loop.  However, if
  !  the problem has inherant i-symmetry, it is necessary to perform the
  !  checkerboarding on the j-loop, and if there is both i-symmetry and
  !  j-symmetry, on the k-loop.  The variables inc* and sw*, *=i,j,k, are
  !  introduced to perform the checkerboarding on the appropriate loop.
  anorm = 0.0d0
  dx3a2i = 1.0 / patch%ds(3)**2
  !do sweep=1,2
  do sweep=1,1
    sw1 = sweep - 1
    sw2 = sw1
    swk = (inck-1) * sw1
    !do k=l(3)+swk,u(3),inck
    do k=l(3),u(3)
      sw3 = sw2
      swj = (incj-1) * sw2
      !do j=l(2)+swj,u(2),incj
      do j=l(2),u(2)
        swi = (inci-1) * sw3
        !do i=l(1)+swi,u(1),inci
        do i=l(1),u(1)
          da = self%eta(i,j,k) * phi(i,j,k) - FourPiGrho(i,j,k) &
             + ( phi(i+ione,j,k) * self%dx1p2i(i) + phi(i-ione,j,k) * self%dx1m2i(i) ) &
               * ione &
             + ( phi(i,j+jone,k) * self%dx2p2i(j) + phi(i,j-jone,k) * self%dx2m2i(j) ) &
               / m1%h2c(i)**2 * jone &
             + ( phi(i,j,k+kone) + phi(i,j,k-kone) ) * dx3a2i &
               / m1%h31c(i)**2 / m2%h32c(j)**2 * kone
          resid(i,j,k) = da
          phi(i,j,k)   = phi(i,j,k) - self%omega_sor / self%eta(i,j,k) * da
          ! this presents the error as drho / rho, where drho (= resid) can be interpreted
          ! as the error in the mass.
          anorm = max(anorm, abs(da) / d(i,j,k))
        end do
        sw3 = 1 - sw3
      end do
      sw2 = 1 - sw2
    end do
    !  Optimise omega using "Chebyshev acceleration".
    self%omega_sor = 1.0 / ( 1.0 - self%rjac2 * 0.25 * self%omega_sor )
  end do

  end associate
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE zeus_sor

!> =============================================================================
!> Conjugate gradient Poisson solver (Ghysels & Vanroose 2014, Algorithm 1)
!> =============================================================================
SUBROUTINE cg (self, patch, phi1, source, error0, error, niter, detailed_timer)
  class(poisson_t)       :: self
  class(patch_t)         :: patch
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: phi1
  real(kind=KindScalarVar), dimension(:,:,:), pointer, intent(in):: source
  real                   :: error, error0
  integer                :: niter
  logical                :: detailed_timer
  !.............................................................................
  real(8)                :: pap, ru, ru2, alpha, beta, ds(3)
  integer                :: iter, m(3)
  integer, save          :: itimer=0
  real(8), dimension(:,:,:), allocatable:: phi, r, p, s, u
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin('poisson_mod%cg', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Allocate double precision scratch arrays
  !-----------------------------------------------------------------------------
  m = shape(phi1)
  allocate (r(m(1),m(2),m(3)), p(m(1),m(2),m(3)))
  allocate (s(m(1),m(2),m(3)), u(m(1),m(2),m(3)))
  allocate (phi(m(1),m(2),m(3)))
  !-----------------------------------------------------------------------------
  ! Calculate the initial residual = source - Laplace(phi)
  !-----------------------------------------------------------------------------
  ds = patch%ds
  phi = phi1
  call periodic_f (phi)
  call laplace (ds, phi, r)
  r = source - r
  error0 = patch%fmaxval(abs(r/(source+floor)))
  if (verbose > 1) print "(i6,2x,a,': iter = ',i4.1,3x,'error =',1p,g12.4)", &
    patch%id, trim(self%kind), 0, error0
  call precondition (r, u)
  ru = sum(r*u)
  !-----------------------------------------------------------------------------
  ! Iteratively improve the solution
  !-----------------------------------------------------------------------------
  p = u
  do iter = 1,max_iter
    call laplace (ds, p, s)
    alpha = ru/sum(s*p)
    phi = phi + alpha*p
    r = r - alpha*s
    error = patch%fmaxval(abs(r/(source+floor)))
    call precondition (r, u)
    ru2 = sum(r*u)
    if (error < tolerance) exit
    beta = ru2/ru
    p = u + beta*p
    ru = ru2
    if (verbose > 1) print "(i6,2x,a,': iter = ',i4.1,3x,'error =',1p,g12.4)", &
      patch%id, trim(self%kind), iter, error
  end do
  niter = iter
  if (iter==max_iter+1) then
    print *, 'poisson_t%cg: no convergence, niter, error =', max_iter, error
  else if (verbose>0 .and. io%master) then
    print *, 'poisson_t%cg: niter, error =', iter, error
  end if
  phi1 = phi
  if (detailed_timer) call trace%end (itimer)
contains
  !> =============================================================================
  !> Make a function period in the guard zones
  !> =============================================================================
  subroutine periodic_f1 (f)
    real(8):: f(:,:,:)
    class(mesh_t), pointer:: m(:)
    !.............................................................................
    if (patch%is_clear (bits%root_grid)) return
    m => patch%mesh
    f(m(1)%lb:m(1)%lo,:,:) = f(m(1)%lb+m(1)%n:m(1)%lo+m(1)%n,:,:)
    f(m(1)%uo:m(1)%ub,:,:) = f(m(1)%uo-m(1)%n:m(1)%ub-m(1)%n,:,:)
    f(:,m(2)%lb:m(2)%lo,:) = f(:,m(2)%lb+m(2)%n:m(2)%lo+m(2)%n,:)
    f(:,m(2)%uo:m(2)%ub,:) = f(:,m(2)%uo-m(2)%n:m(2)%ub-m(2)%n,:)
    f(:,:,m(3)%lb:m(3)%lo) = f(:,:,m(3)%lb+m(3)%n:m(3)%lo+m(3)%n)
    f(:,:,m(3)%uo:m(3)%ub) = f(:,:,m(3)%uo-m(3)%n:m(3)%ub-m(3)%n)
  end subroutine periodic_f1
  !> =============================================================================
  !> Make a function period in the guard zones
  !> =============================================================================
  subroutine periodic_f (f)
    real(8):: f(:,:,:)
    class(mesh_t), pointer:: m(:)
    integer:: ix, iy, iz
    integer, dimension(3):: lb, ub, li, ui, n, i, j
    !.............................................................................
    if (patch%is_clear (bits%root_grid)) return
    m => patch%mesh
    lb = m%lb; ub=m%ub
    li = m%li; ui=m%ui
    n = m%n
    do iz=lb(3),ub(3)
    do iy=lb(2),ub(2)
    do ix=lb(1),ub(1)
      i = [ix,iy,iz]
      if (all(i >= li .and. i <= ui)) cycle
      where (n > 0)
        j = modulo(i-li,n)+li
      else where
        j = 1
      end where
      f(ix,iy,iz) = f(j(1),j(2),j(3))
      f(ix,iy,iz) = f(j(1),j(2),j(3))
    end do
    end do
    end do
  end subroutine periodic_f
  !> =============================================================================
  !> Set scratch array outer boundary to zero
  !> =============================================================================
  subroutine zero_boundaries (f)
    real(8):: f(:,:,:)
    class(mesh_t), pointer:: m(:)
    !.............................................................................
    m => patch%mesh
    f(m(1)%lb:m(1)%lo,:,:) = 0d0
    f(m(1)%uo:m(1)%ub,:,:) = 0d0
    f(:,m(2)%lb:m(2)%lo,:) = 0d0
    f(:,m(2)%uo:m(2)%ub,:) = 0d0
    f(:,:,m(3)%lb:m(3)%lo) = 0d0
    f(:,:,m(3)%uo:m(3)%ub) = 0d0
  end subroutine zero_boundaries
  !> =============================================================================
  !> 2nd order Laplace operator
  !> =============================================================================
  subroutine laplace (ds, phi, out)
    real(8)                   :: ds(3)
    real(8), dimension(:,:,:) :: phi, out
    !.............................................................................
    integer                   :: ix, iy, iz
    real(8)                   :: c(3)
    !-----------------------------------------------------------------------------
    c = 1d0/ds**2
    do iz=patch%li(3),patch%ui(3)
    do iy=patch%li(2),patch%ui(2)
    do ix=patch%li(1),patch%ui(1)
      out(ix,iy,iz) = c(1)*(phi(ix+1,iy,iz)+phi(ix-1,iy,iz)-2d0*phi(ix,iy,iz)) &
                    + c(2)*(phi(ix,iy+1,iz)+phi(ix,iy-1,iz)-2d0*phi(ix,iy,iz)) &
                    + c(3)*(phi(ix,iy,iz+1)+phi(ix,iy,iz-1)-2d0*phi(ix,iy,iz))
    end do
    end do
    end do
    call periodic_f (out)
  end subroutine laplace
  !> =============================================================================
  !> Preconditioner for conjugate gradient method
  !> =============================================================================
  subroutine precondition (res, out)
    real(8), dimension(:,:,:) :: res, out
    !.............................................................................
    integer                   :: ix, iy, iz, m(3)
    real(8)                   :: a, b, c
    !-----------------------------------------------------------------------------
    m = shape(res)
    a = 1d0 + 1d0/12d0
    b = 1d0/6d0
    do iz=patch%li(3),patch%ui(3)
    do iy=patch%li(2),patch%ui(2)
    do ix=patch%li(1),patch%ui(1)
      out(ix,iy,iz) = b*(res(ix+1,iy,iz)+res(ix-1,iy,iz)  &
                    +    res(ix,iy+1,iz)+res(ix,iy-1,iz)  &
                    +    res(ix,iy,iz+1)+res(ix,iy,iz-1)) &
                    + a*res(ix,iy,iz)
    end do
    end do
    end do
    call periodic_f (out)
  end subroutine precondition
END SUBROUTINE cg

!> =============================================================================
!> 2nd order Laplace operator
!> =============================================================================
SUBROUTINE laplace (ds, phi, source, res)
  real(8)                   :: ds(3)
  real, dimension(:,:,:)    :: phi, source, res
  !.............................................................................
  integer                   :: ix, iy, iz, m(3)
  real(8)                   :: c(3)
  !-----------------------------------------------------------------------------
  c = 1d0/ds**2
  m = shape(phi)
  do iz=2,m(3)-1
  do iy=2,m(2)-1
  do ix=2,m(1)-1
    res(ix,iy,iz) = c(1)*(phi(ix+1,iy,iz)+phi(ix-1,iy,iz)-2d0*phi(ix,iy,iz)) &
                  + c(2)*(phi(ix,iy+1,iz)+phi(ix,iy-1,iz)-2d0*phi(ix,iy,iz)) &
                  + c(3)*(phi(ix,iy,iz+1)+phi(ix,iy,iz-1)-2d0*phi(ix,iy,iz)) &
                  - source(ix,iy,iz)
  end do
  end do
  end do
END SUBROUTINE laplace

!> =============================================================================
!> Given a grav. pot. and a density, calculate Poisson's equation and return
!> the residual.
!> =============================================================================
SUBROUTINE check_poisson (self, patch, FourPiG, resid)
  class(poisson_t):: self
  class(patch_t):: patch
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: phi
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d
  real(kind=KindScalarVar), dimension(patch%gn(1),patch%gn(2),patch%gn(3)):: resid
  real:: FourPiG
  integer:: i, j, k, ione, jone, kone
  real, dimension(3):: dsi2
  !.............................................................................
  call trace%begin('experiment_t%check_poisson')
  d => patch%mem(:,:,:,patch%idx%d,patch%it,1)
  phi => patch%mem(:,:,:,patch%idx%phi,patch%it,1)
  ione = 1
  jone = 1
  kone = 1
  if (patch%n(1) <= 1) ione = 0
  if (patch%n(2) <= 1) jone = 0
  if (patch%n(3) <= 1) kone = 0
  dsi2 = 1.0 / patch%ds**2
  resid = 0.0
  associate (l => patch%li, u => patch%ui)
  do k=l(3),u(3)
    do j=l(2),u(2)
      do i=l(1),u(1)
        resid(i,j,k) = (phi(i+ione,j,k) - 2.0 * phi(i,j,k) + phi(i-ione,j,k)) * dsi2(1) * ione &
                     + (phi(i,j+jone,k) - 2.0 * phi(i,j,k) + phi(i,j-jone,k)) * dsi2(2) * jone &
                     + (phi(i,j,k+kone) - 2.0 * phi(i,j,k) + phi(i,j,k-kone)) * dsi2(3) * kone &
                     - FourPiG * d(i,j,k)
        resid(i,j,k) = abs(resid(i,j,k))
      end do
    end do
  end do
  end associate
  print*,"Max. resid.:",maxval(resid)
  call trace%end
END SUBROUTINE check_poisson

END MODULE poisson_mod
