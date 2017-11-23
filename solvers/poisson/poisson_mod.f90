!> =============================================================================
!> Methods for solving the Poisson equation; partially derived from routines
!> developed by Troels Haugbølle for the course Computational Astrophysics.
!> =============================================================================
MODULE poisson_mod
  USE mpi_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE patch_mod
  USE mesh_mod
  USE bits_mod
  USE scaling_mod
  USE omp_timer_mod
  implicit none
  private
  type, public:: poisson_t
    real(8):: omega_sor ! Overrelaxation parameter
    real(8):: rjac2     ! Jacobi spectral radius squared
    procedure(cg), pointer:: update => null()
  contains
    procedure:: init
    procedure:: cg
    procedure:: sor
    procedure, nopass:: laplace
  end type
  real,    save:: floor=1e-2, tolerance=1e-4, fourPiG=1.0
  integer, save:: verbose=0
  integer, save:: max_iter=500
  logical, save:: precondition=.true.
  real(8), parameter:: pi = asin(1.0) * 2.0
CONTAINS

!> =============================================================================
!> Initialize poisson_params
!> =============================================================================
SUBROUTINE init (self)
  implicit none
  class(poisson_t)  :: self
  !.............................................................................
  integer           :: iostat
  !.............................................................................
  character(len=12), save :: solver='cg'
  logical,           save :: first_time = .true.
  namelist /poisson_params/ tolerance, floor, max_iter, fourPiG, precondition, verbose, solver
  !-----------------------------------------------------------------------------
  !$omp critical (poisson_cr)
  if (first_time) then
    first_time = .false.
    fourpiG = scaling%grav
    rewind (io%input)
    read (io%input, poisson_params, iostat=iostat)
    if (io%master) write (io%output, poisson_params)
  end if
  !$omp end critical (poisson_cr)
  select case(trim(solver))
  case ('cg')
    self%update => cg
  case ('sor')
    self%update => sor
  case default
    if (io%master) then
      write(*,*) 'Poisson solver '//trim(solver)//' is unknown'
      write(*,*) 'Valid solvers are: cg, sor'
      call mpi%abort()
    end if
  end select
END SUBROUTINE init

!> =============================================================================
!> Conjugate gradient Poisson solver (Ghysels & Vanroose 2014, Algorithm 1)
!>   http://www.sciencedirect.com/science/article/pii/S0167819113000719
!> With the "incomplete Poisson" preconditioner from:
!> http://www.vis.uni-stuttgart.de/~weiskopf/publications/pdp10.pdf 
!> or http://ieeexplore.ieee.org/document/5452414/
!> See also the CG method in RAMSES: `ramses/poisson/phi_fine_cg.f90`.
!> =============================================================================
SUBROUTINE cg (self, patch, phi1, d, d0, detailed_timer)
  class(poisson_t)       :: self
  class(patch_t)         :: patch
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: phi1
  real(kind=KindScalarVar), dimension(:,:,:), pointer, intent(in):: d
  real(8)                :: d0
  logical                :: detailed_timer
  !.............................................................................
  integer                :: niter
  real                   :: error, error0, usppt
  real(8)                :: pap, ru, ru2, alpha, beta, ds(3), wc
  integer                :: iter, m(3), ione, jone, kone
  integer, save          :: itimer=0
  real(8), dimension(:,:,:), allocatable:: source, phi, r, p, s, u
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin('poisson_mod%cg', itimer=itimer)
  ione = 1
  jone = 1
  kone = 1
  if (size(phi1,1) <= 1) ione = 0
  if (size(phi1,2) <= 1) jone = 0
  if (size(phi1,3) <= 1) kone = 0
  !-----------------------------------------------------------------------------
  ! Allocate double precision scratch arrays
  !-----------------------------------------------------------------------------
  m = shape(phi1)
  allocate (r(m(1),m(2),m(3)), p(m(1),m(2),m(3)))
  allocate (s(m(1),m(2),m(3)), u(m(1),m(2),m(3)))
  allocate (source(m(1),m(2),m(3)), phi(m(1),m(2),m(3)))
  !-----------------------------------------------------------------------------
  ! Calculate the initial residual = source - Laplace(phi)
  !-----------------------------------------------------------------------------
  wc = wallclock()
  ds = patch%ds
  phi = phi1
  source = fourPiG*(d-d0)
  call zero_boundaries (source)
  p = 0d0
  r = 0d0
  s = 0d0
  u = 0d0
  call laplace (ds, phi, r)
  r = source - r
  error0 = patch%fmaxval(abs(r/(d+floor)))
  if (verbose > 1) print "(i6,2x,'cg: iter = ',i4.1,3x,'error =',1p,g12.4)", &
    patch%id, 0, error0
  call precond (r, u)
  ru = patch%fsum(r*u)
  !-----------------------------------------------------------------------------
  ! Iteratively improve the solution
  !-----------------------------------------------------------------------------
  p = u
  do iter = 1,max_iter
    call laplace (ds, p, s)
    alpha = ru/patch%fsum(s*p)
    phi = phi + alpha*p
    r   = r   - alpha*s
    error = patch%fmaxval(abs(r/(d+floor)))/fourPiG
    call precond (r, u)
    ru2 = patch%fsum(r*u)
    if (error < tolerance) exit
    beta = ru2/ru
    p = u + beta*p
    ru = ru2
    if (verbose > 1) print "(i6,2x,'cg: iter = ',i4.1,3x,'error =',1p,g12.4)", &
      patch%id, iter, error
    !!!write (11) phi
  end do
  if (iter==max_iter+1) then
    print *, 'poisson_t%cg: no convergence, niter, error =', max_iter, error
  else if (verbose>0 .and. io%master) then
    usppt = 1e6*(wallclock()-wc)/product(shape(phi))
    print '(a,i7,i5,3x,2e12.3,f8.2)', &
     'poisson_t%cg: id, niter, initial/final error, us/pt =', &
      patch%id, iter, error0, error, usppt
  end if
  phi1 = phi
  if (detailed_timer) call trace%end (itimer)
contains
  !> =============================================================================
  !> Set scratch array outer boundary to zero
  !> =============================================================================
  subroutine zero_boundaries (f)
    real(8):: f(:,:,:)
    class(mesh_t), pointer:: m(:)
    !.............................................................................
    if (patch%is_periodic()) return
    m => patch%mesh
    if (size(f,1) > 1) then
      f(m(1)%lb:m(1)%lo,:,:) = 0d0
      f(m(1)%uo:m(1)%ub,:,:) = 0d0
    end if
    if (size(f,2) > 1) then
      f(:,m(2)%lb:m(2)%lo,:) = 0d0
      f(:,m(2)%uo:m(2)%ub,:) = 0d0
    end if
    if (size(f,3) > 1) then
      f(:,:,m(3)%lb:m(3)%lo) = 0d0
      f(:,:,m(3)%uo:m(3)%ub) = 0d0
    end if
  end subroutine zero_boundaries
  !> =============================================================================
  !> 2nd order Laplace operator
  !> =============================================================================
  subroutine laplace (ds, phi, out)
    real(8)                   :: ds(3)
    real(8), dimension(:,:,:) :: phi, out
    !.............................................................................
    integer                   :: i, j, k, im1, jm1, km1, ip1, jp1, kp1, li(3), n(3)
    real(8)                   :: c(3)
    !-----------------------------------------------------------------------------
    c  = 1d0/ds**2
    n  = patch%mesh%nc
    where (n == 1)
      c = 0.0
    end where
    li = patch%mesh%li
    !-----------------------------------------------------------------------------
    ! For periodic patches, wrap around
    !-----------------------------------------------------------------------------
    if (patch%is_periodic()) then
      do k=patch%mesh(3)%lb,patch%mesh(3)%ub
        kp1 = modulo(k-li(3)+1,n(3))+li(3)
        km1 = modulo(k-li(3)-1,n(3))+li(3)
        do j=patch%mesh(2)%lb,patch%mesh(2)%ub
          jp1 = modulo(j-li(2)+1,n(2))+li(2)
          jm1 = modulo(j-li(2)-1,n(2))+li(2)
          do i=patch%mesh(1)%lb,patch%mesh(1)%lo
            ip1 = modulo(i-li(1)+1,n(1))+li(1)
            im1 = modulo(i-li(1)-1,n(1))+li(1)
            out(i,j,k) = c(1)*(phi(ip1,j,k)+phi(im1,j,k)-2d0*phi(i,j,k)) &
                       + c(2)*(phi(i,jp1,k)+phi(i,jm1,k)-2d0*phi(i,j,k)) &
                       + c(3)*(phi(i,j,kp1)+phi(i,j,km1)-2d0*phi(i,j,k))
          end do
          do i=patch%mesh(1)%li,patch%mesh(1)%ui
            out(i,j,k) = c(1)*(phi(i+ione,j,k)+phi(i-ione,j,k)-2d0*phi(i,j,k)) &
                       + c(2)*(phi(i   ,jp1,k)+phi(i   ,jm1,k)-2d0*phi(i,j,k)) &
                       + c(3)*(phi(i   ,j,kp1)+phi(i   ,j,km1)-2d0*phi(i,j,k))
          end do
          do i=patch%mesh(1)%uo,patch%mesh(1)%ub
            ip1 = modulo(i-li(1)+1,n(1))+li(1)
            im1 = modulo(i-li(1)-1,n(1))+li(1)
            out(i,j,k) = c(1)*(phi(ip1,j,k)+phi(im1,j,k)-2d0*phi(i,j,k)) &
                       + c(2)*(phi(i,jp1,k)+phi(i,jm1,k)-2d0*phi(i,j,k)) &
                       + c(3)*(phi(i,j,kp1)+phi(i,j,km1)-2d0*phi(i,j,k))
          end do
        end do
      end do
    else
    !-----------------------------------------------------------------------------
    ! For non-periodic patches, loop from li to ui
    !-----------------------------------------------------------------------------
      call zero_boundaries (out)
      do k=patch%mesh(3)%li,patch%mesh(3)%ui
        do j=patch%mesh(2)%li,patch%mesh(2)%ui
          do i=patch%mesh(1)%li,patch%mesh(1)%ui
            out(i,j,k) = c(1)*(phi(i+ione,j,k)+phi(i-ione,j,k)-2d0*phi(i,j,k)) &
                       + c(2)*(phi(i,j+jone,k)+phi(i,j-jone,k)-2d0*phi(i,j,k)) &
                       + c(3)*(phi(i,j,k+kone)+phi(i,j,k-kone)-2d0*phi(i,j,k))
          end do
        end do
      end do
    end if
  end subroutine laplace
  !> =============================================================================
  !> Preconditioner for conjugate gradient method
  !> =============================================================================
  subroutine precond (res, out)
    real(8), dimension(:,:,:) :: res, out
    !.............................................................................
    integer                   :: i, j, k, im1, jm1, km1, ip1, jp1, kp1, li(3), n(3), &
                                 ndim
    real(8)                   :: a, b(3)
    !-----------------------------------------------------------------------------
    if (precondition) then
      n  = patch%mesh%nc
      ndim = sum(merge(1,0,n > 1))
      b = 0.5d0 / ndim
      a = 1d0 + 0.25d0 / ndim
      where (n == 1)
        b = 0.0d0
      end where
      li = patch%mesh%li
      !---------------------------------------------------------------------------
      ! For periodic patches, wrap around (cheap patch so no need to optimize)
      !---------------------------------------------------------------------------
      if (patch%is_periodic()) then
        do k=patch%mesh(3)%lb,patch%mesh(3)%ub
          kp1 = modulo(k-li(3)+1,n(3))+li(3)
          km1 = modulo(k-li(3)-1,n(3))+li(3)
          do j=patch%mesh(2)%lb,patch%mesh(2)%ub
            jp1 = modulo(j-li(2)+1,n(2))+li(2)
            jm1 = modulo(j-li(2)-1,n(2))+li(2)
            do i=patch%mesh(1)%lb,patch%mesh(1)%ub
              ip1 = modulo(i-li(1)+1,n(1))+li(1)
              im1 = modulo(i-li(1)-1,n(1))+li(1)
              out(i,j,k) = b(1)*(res(ip1,j,k)+res(im1,j,k))  &
                         + b(2)*(res(i,jp1,k)+res(i,jm1,k))  &
                         + b(3)*(res(i,j,kp1)+res(i,j,km1)) &
                         + a*res(i,j,k)
            end do
          end do
        end do
      !---------------------------------------------------------------------------
      ! For non- periodic patches, loop over internal cells
      !---------------------------------------------------------------------------
      else
        do k=patch%li(3),patch%ui(3)
          do j=patch%li(2),patch%ui(2)
            do i=patch%li(1),patch%ui(1)
              out(i,j,k) = b(1)*(res(i+ione,j,k)+res(i-ione,j,k))  &
                         + b(2)*(res(i,j+jone,k)+res(i,j-jone,k))  &
                         + b(3)*(res(i,j,k+kone)+res(i,j,k-kone)) &
                         + a*res(i,j,k)
            end do
          end do
        end do
        call zero_boundaries (out)
      end if
    else
      out = res
    end if
  end subroutine precond
END SUBROUTINE cg

!> =============================================================================
!> 2nd order successive over relaxation with red-black ordering
!> =============================================================================
SUBROUTINE sor (self, patch, phi1, d, d0, detailed_timer)
  class(poisson_t)       :: self
  class(patch_t)         :: patch
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: phi1
  real(kind=KindScalarVar), dimension(:,:,:), pointer, intent(in):: d
  real(8):: d0
  logical                :: detailed_timer
  !.............................................................................
  integer                :: i, j, k, sweep, ss, ione, jone, kone, m(3), iter, niter
  real                   :: error, error0
  real(8)                :: a, b, b1, b2, b3, res, numer, denom, wc, usppt
  integer, save          :: itimer=0
  real(8), dimension(:,:,:), allocatable:: phi, source

  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('poisson_t%sor_iter', itimer=itimer)

  associate(n=>patch%n,ds=>patch%ds,l=>patch%li,u=>patch%ui,gn=>patch%gn)
  allocate (phi(gn(1),gn(2),gn(3)), source(gn(1),gn(2),gn(3)))

  !-------------------------------------------------------------------
  ! Calculate optimal guess for the SOR prefactor depending on
  ! dimensionality, ds, and grid points n:
  !    delta = (1/sum(1/ds^2)) sum(1/ds^2 * cos(pi/n))
  !    omega_sor = 2 / (1 + (1 - delta^2)^0.5)
  ! For an even dimensionality grid
  !    delta = cos(pi / n)
  ! which for n > 10 is well approximated by
  !    delta ~ 1 - 1/2 (pi/n)**2
  ! and then
  !    omega_sor = 2 / (1 + pi/n)
  !-------------------------------------------------------------------
  numer = 0.0
  denom = 0.0
  do i=1,3
    if (patch%n(i) > 1) then
      numer = numer + cos(pi / patch%n(i))
      denom = denom + 1d0
    end if
  end do
  self%rjac2 = (numer / denom)**2
  self%omega_sor = 2d0 / (1d0 + sqrt(1d0 - self%rjac2))

  if (patch%is_periodic()) call patch%make_periodic (phi)
  ione = 1
  jone = 1
  kone = 1
  if (n(1) <= 1) ione = 0
  if (n(2) <= 1) jone = 0
  if (n(3) <= 1) kone = 0

  wc = wallclock()
  phi = phi1
  source = fourPiG * (d - d0)
  ! Calculate 1 / [(2/dx + 2/dy + 2/dz)] prefactor according to dimensionality
  a = 0.0
  if (n(1) > 1) a = a + 2. / ds(1)**2
  if (n(2) > 1) a = a + 2. / ds(2)**2
  if (n(3) > 1) a = a + 2. / ds(3)**2
  a = 1. / a
  b1 = a / ds(1)**2 * ione
  b2 = a / ds(2)**2 * jone
  b3 = a / ds(3)**2 * kone

  do iter=1,max_iter
    error = 0.0
    do sweep=0,1
      do k=l(3),u(3)
        do j=l(2),u(2)
          ss = modulo(j + k + sweep,2)
          do i=l(1)+ss,u(1),2
            !---------------------------------------------------------------------
            ! This form of the residual has no factor in front of phi, which makes
            ! it clear how to apply the over-relaxation factor
            !---------------------------------------------------------------------
            res = b1 * (phi(i+ione,j,k) + phi(i-ione,j,k)) &
                + b2 * (phi(i,j+jone,k) + phi(i,j-jone,k)) &
                + b3 * (phi(i,j,k+kone) + phi(i,j,k-kone)) &
                - phi(i,j,k)- a * source(i,j,k)
            phi(i,j,k) = phi(i,j,k) + self%omega_sor * res
            error = max(abs(res / (d(i,j,k) + floor)), error)
          end do
        end do
      end do
    end do
    !-----------------------------------------------------------------------------
    ! We are supposed to return the residual as (Laplace(phi)-source)/d, so we
    ! need to divide with the factor that source (FourPiGrho) was multiplied by
    !-----------------------------------------------------------------------------
    error = error / a
    if (iter==1) error0 = error
    if (verbose > 1) print *, 'poisson_t%sor: iter, error =', iter, error
    if (error < tolerance) exit
    !if (patch%is_periodic()) call patch%make_periodic (phi)
  end do

  if (iter==max_iter+1) then
    print *, 'poisson_t%sor: no convergence, niter, error =', max_iter, error
  else if (verbose > 0 .and. io%master) then
    usppt = 1e6 * (wallclock() - wc) / product(shape(phi))
    print '(a,i7,i5,3x,2e12.3,f8.2)', &
     'poisson_t%sor: id, niter, initial/final error, us/pt =', &
      patch%id, iter, error0, error, usppt
  end if
  niter = iter
  phi1 = phi
  deallocate (phi)

  end associate
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE sor

!> =============================================================================
!> 2nd order Laplace operator
!> =============================================================================
SUBROUTINE laplace (ds, phi, source, res)
  real(8)                   :: ds(3)
  real, dimension(:,:,:)    :: phi, source, res
  !.............................................................................
  integer                   :: ix, iy, iz, m(3), ione, jone, kone
  real(8)                   :: c(3)
  !-----------------------------------------------------------------------------
  ione = 1
  jone = 1
  kone = 1
  if (size(phi,1) <= 1) ione = 0
  if (size(phi,2) <= 1) jone = 0
  if (size(phi,3) <= 1) kone = 0

  c = 1d0/ds**2
  m = shape(phi)
  where (m == 1)
    c = 0.0d0
  end where
  do iz=2,m(3)-kone
  do iy=2,m(2)-jone
  do ix=2,m(1)-ione
    res(ix,iy,iz) = c(1)*(phi(ix+ione,iy,iz)+phi(ix-ione,iy,iz)-2d0*phi(ix,iy,iz)) &
                  + c(2)*(phi(ix,iy+jone,iz)+phi(ix,iy-jone,iz)-2d0*phi(ix,iy,iz)) &
                  + c(3)*(phi(ix,iy,iz+kone)+phi(ix,iy,iz-kone)-2d0*phi(ix,iy,iz)) &
                  - source(ix,iy,iz)
  end do
  end do
  end do
END SUBROUTINE laplace

END MODULE poisson_mod
