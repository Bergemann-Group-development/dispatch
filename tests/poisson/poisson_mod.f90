!> =============================================================================
!> Poisson solver module, which can solve a Poisson problem on a patch in a
!> reasonably small number of iterations.  Larger problems must be solved
!> separately on each scale, imposing the correction from the larger scale
!> solution before iterating on the smaller scales.  In general, the patches
!> may have arbitrary location relative to each other, and shoul use the same
!> download procedure as is used in other parts of the evolution.
!>
!> As with the radiative transfer solution, the Poisson solver is called when
!> other solvers are ready to update; i.e., when guard zone values are known.
!> The solution generally consists of small corrections to previous solutions,
!> which are remember from the previous time step.  In the case where patches
!> are moving, it is still assumed that the local solution remains to some 
!> extend relevant from timestep to timestep, as is indeed the case if the 
!> motion of the patch follows the motion of its center of mass.
!>
!> In order to correcly correct the initial value of phi for larger scales,
!> one shold subtract off the larger scale component of phi before saving the
!> solution until the next time step.  The larger scale phi is then added to
!> phi, which means it is effectively corrected for the change of the larger
!> scale component, which may have a different center of mass motion.
!>
!> In principle the smaller scale solution should be fed back into the next
!> larger scale, by uploading the smaller scale phi, but we can also choose
!> to consider the larger scale solution, based on the mass resolved on that
!> scale, as being authoritative, and providing fixed boundary conditions for
!> the smaller scale solutions.
!>
!> It is by no means obvious that imposing a phi from the smaller scale solution
!> inside the larger scale one will result in a more precise solution. A case
!> in point is when solving for the gravitational potential of a central mass,
!> which is fully contained in the innermost patch. The correct solution in
!> larger scale patches then satisfies laplace(phi)=0. If one imposes the phi
!> resolved on smaller scales one should expect large truncation errors, which
!> means laplace(phi)=0 is NOT satisfied when part of the phi stencil rests on
!> the higher resolution part.
!>
!> The Rubik's Cube setup is a good test example!
!> =============================================================================
MODULE poisson_mod
  USE io_unit_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  implicit none
  private
  type, public :: poisson_t
    real    :: grav = 1.0
    real    :: sor = 1.0
    real    :: tolerance = 1e-4
    integer :: maxiter = 100
  contains
    procedure :: init
    procedure :: solve
    !procedure :: update => jacobi
    procedure :: update => gauss_seidel
    procedure :: laplace
  end type
  real, parameter :: fourpi=16.*atan(1.)
CONTAINS

!> -----------------------------------------------------------------------------
!> Choose solver and parameters
!> -----------------------------------------------------------------------------
SUBROUTINE init (self, mesh)
  class(poisson_t) :: self
  class(mesh_t), dimension(:), pointer, intent(in)  :: mesh
  real,    save   :: tolerance=1e-4, sor=1.0, grav=1.0
  integer, save   :: maxiter=100
  logical, save   :: first_time=.true.
  real, parameter :: pi=4.*atan(1.)
  real            :: delta, pref
  integer         :: i
  character(len=16) :: comment
  namelist /poisson_params/ grav, tolerance, maxiter, sor, comment
  !.............................................................................

  !$omp critical (poisson_cr)
  if (first_time) then
    first_time = .false.
    rewind (io_unit%input); read (io_unit%input, poisson_params)
    if (sor==0.) then
      delta = 0.0; pref = 0.0
      do i=1,3
        if (mesh(i)%n > 1) then
          pref  = pref + 1.
          delta = delta + cos(pi / mesh(i)%n)
        endif
      enddo
      delta = delta / pref
      sor = 2.0 / (1. + sqrt(1.0 - delta**2))
      comment =' (computed)'
    else
      comment =' (from namelist)'
    endif
    if (io_unit%master) write (*, poisson_params)
  end if
  !$omp end critical (poisson_cr)
  self%tolerance = tolerance
  self%maxiter   = maxiter
  self%sor       = sor
  self%grav      = grav
END SUBROUTINE init

!> -----------------------------------------------------------------------------
!> Given rho and the initial values of phi, iterate the solver method on the
!> internal part of phi.
!> -----------------------------------------------------------------------------
SUBROUTINE solve (self, mesh, rho, phi, iters)
  class(poisson_t) :: self
  class(mesh_t), dimension(:), pointer, intent(in)  :: mesh
  integer :: iters
  real, dimension(:,:,:),      pointer, intent(in)  :: rho
  real, dimension(:,:,:),      pointer, intent(out) :: phi
  real, dimension(:,:,:),      pointer              :: residual
  real    :: rho0, error
  integer :: iter, ix, iy, iz, l(3), u(3)
  integer, save :: itimer=0
  !.............................................................................
  call trace%begin ('poisson_mod::solve', itimer)
  allocate (residual(mesh(1)%gn,mesh(2)%gn,mesh(3)%gn))
  rho0 = maxval(rho)
  do iter=1,self%maxiter
    call self%laplace (mesh, phi, rho, residual)
    l = mesh%li
    u = mesh%ui
    error = maxval(abs(residual(l(1):u(1),l(2):u(2),l(3):u(3))))
    if (io%verbose > 2) then
      write (io_unit%log, '(a,i5,1p,e14.6)') 'poisson_t%solve: iter, error =', iter, error
      call flush (io_unit%log)
    end if
    if (error < self%tolerance*rho0) then
      iters = iter
      exit
    end if
    call self%update (mesh, phi, residual)
  end do
  deallocate (residual)
  call trace%end (itimer)
END SUBROUTINE solve

!> -----------------------------------------------------------------------------
!> Laplace operator, internal for now -- could be external
!> -----------------------------------------------------------------------------
SUBROUTINE laplace (self, mesh, phi, rho, out)
  class(poisson_t) :: self
  class(mesh_t), dimension(:), pointer, intent(in)  :: mesh
  real, dimension(:,:,:),      pointer, intent(in)  :: rho
  real, dimension(:,:,:),      pointer, intent(in)  :: phi
  real, dimension(:,:,:),      pointer, intent(out) :: out     
  real :: c(3), f
  integer :: ix, iy, iz
  !.............................................................................
  c = 1./mesh%d**2
  f = self%grav*fourpi
  do iz=mesh(3)%li,mesh(3)%ui
  do iy=mesh(2)%li,mesh(2)%ui
  do ix=mesh(1)%li,mesh(1)%ui
    out(ix,iy,iz) = &
      c(3)*((phi(ix  ,iy  ,iz-1)-phi(ix,iy,iz)) + (phi(ix  ,iy  ,iz+1)-phi(ix,iy,iz))) + &
      c(2)*((phi(ix  ,iy-1,iz  )-phi(ix,iy,iz)) + (phi(ix  ,iy+1,iz  )-phi(ix,iy,iz))) + &
      c(1)*((phi(ix-1,iy  ,iz  )-phi(ix,iy,iz)) + (phi(ix+1,iy  ,iz  )-phi(ix,iy,iz))) - &
      f*rho(ix,iy,iz)
  end do
  end do
  end do
END SUBROUTINE laplace

!> -----------------------------------------------------------------------------
!> Make one iteration on the Poisson problem
!> -----------------------------------------------------------------------------
SUBROUTINE jacobi (self, mesh, phi, residual)
  class(poisson_t) :: self
  class(mesh_t), dimension(:), pointer, intent(in)  :: mesh
  real, dimension(:,:,:),      pointer, intent(in)  :: phi
  real, dimension(:,:,:),      pointer              :: residual
  real    :: c
  integer :: ix, iy, iz
  !-----------------------------------------------------------------------------
  ! Assume that the mesh is isotropic
  !-----------------------------------------------------------------------------
  c = 1./sum(2./mesh%d**2)
  do iz=mesh(3)%li,mesh(3)%ui
  do iy=mesh(2)%li,mesh(2)%ui
  do ix=mesh(1)%li,mesh(1)%ui
    phi(ix,iy,iz) = phi(ix,iy,iz) + c*residual(ix,iy,iz)
    !phi(ix,iy,iz) = c*residual(ix,iy,iz)
  end do
  end do
  end do
END SUBROUTINE jacobi

!> -----------------------------------------------------------------------------
!> Make one iteration on the Poisson problem, using red-black Gauss-Seidel
!> -----------------------------------------------------------------------------
SUBROUTINE gauss_seidel (self, mesh, phi, residual)
  class(poisson_t) :: self
  class(mesh_t), dimension(:), pointer, intent(in)  :: mesh
  real, dimension(:,:,:),      pointer, intent(in)  :: phi
  real, dimension(:,:,:),      pointer              :: residual
  real    :: c
  integer :: ix, iy, iz, ss
  integer, save :: sweep=0
  !-----------------------------------------------------------------------------
  ! Assume that the mesh is isotropic
  !-----------------------------------------------------------------------------
  c = self%sor/sum(2./mesh%d**2)
  sweep = mod(sweep+1,2)
  do iz=mesh(3)%li,mesh(3)%ui
  do iy=mesh(2)%li,mesh(2)%ui
    ss = mod(iy+iz+sweep,2)
    do ix=mesh(1)%li+ss,mesh(1)%ui,2
      phi(ix,iy,iz) = phi(ix,iy,iz) + c*residual(ix,iy,iz)
    end do
  end do
  end do
END SUBROUTINE gauss_seidel

END MODULE poisson_mod
