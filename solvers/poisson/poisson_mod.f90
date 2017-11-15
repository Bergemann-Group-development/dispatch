!> =============================================================================
!> Methods for solving the Poisson equation; derived from routines developed by
!> Troels Haugbølle for the course Computational Astrophysics, and from papers
!> on conjugate gradient methods cited in ramses/poisson/phi_fine_cg.f90:
!>   http://www.sciencedirect.com/science/article/pii/S0167819113000719
!>   http://www.vis.uni-stuttgart.de/~weiskopf/publications/pdp10.pdf
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
  contains
    procedure:: init
    procedure:: cg
    procedure, nopass:: laplace
  end type
  real,    save:: floor=1e-2, tolerance=1e-4, fourPiG=1.0
  integer, save:: verbose=0
  integer, save:: max_iter=500
  logical, save:: precondition=.true.
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
  real,              save :: omega_sor=0.0
  logical,           save :: first_time = .true.
  namelist /poisson_params/ tolerance, floor, max_iter, fourPiG, precondition, verbose
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
END SUBROUTINE init

!> =============================================================================
!> Conjugate gradient Poisson solver (Ghysels & Vanroose 2014, Algorithm 1)
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
  integer                :: iter, m(3)
  integer, save          :: itimer=0
  real(8), dimension(:,:,:), allocatable:: source, phi, r, p, s, u
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin('poisson_mod%cg', itimer=itimer)
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
    write (11) phi
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
    integer                   :: i, j, k, im1, jm1, km1, ip1, jp1, kp1, li(3), n(3)
    real(8)                   :: c(3)
    !-----------------------------------------------------------------------------
    c  = 1d0/ds**2
    n  = patch%mesh%nc
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
            out(i,j,k) = c(1)*(phi(i+1,j,k)+phi(i-1,j,k)-2d0*phi(i,j,k)) &
                       + c(2)*(phi(i,jp1,k)+phi(i,jm1,k)-2d0*phi(i,j,k)) &
                       + c(3)*(phi(i,j,kp1)+phi(i,j,km1)-2d0*phi(i,j,k))
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
            out(i,j,k) = c(1)*(phi(i+1,j,k)+phi(i-1,j,k)-2d0*phi(i,j,k)) &
                       + c(2)*(phi(i,j+1,k)+phi(i,j-1,k)-2d0*phi(i,j,k)) &
                       + c(3)*(phi(i,j,k+1)+phi(i,j,k-1)-2d0*phi(i,j,k))
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
    integer                   :: i, j, k, im1, jm1, km1, ip1, jp1, kp1, li(3), n(3)
    real(8)                   :: a, b, c
    !-----------------------------------------------------------------------------
    if (precondition) then
      a = 1d0 + 1d0/12d0
      b = 1d0/6d0
      n  = patch%mesh%nc
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
              out(i,j,k) = b*(res(ip1,j,k)+res(im1,j,k)  &
                         +    res(i,jp1,k)+res(i,jm1,k)  &
                         +    res(i,j,kp1)+res(i,j,km1)) &
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
              out(i,j,k) = b*(res(i+1,j,k)+res(i-1,j,k)  &
                         +    res(i,j+1,k)+res(i,j-1,k)  &
                         +    res(i,j,k+1)+res(i,j,k-1)) &
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

END MODULE poisson_mod
