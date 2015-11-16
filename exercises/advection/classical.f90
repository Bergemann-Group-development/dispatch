!> =============================================================================
!> $Id$
!> This version uses classical program & subroutine structure, with all the 
!> data carried in the main program, and subroutines written for the specific
!> purpose of the experiment.
!> =============================================================================
! ------------------------------------------------------------------------------
! The initial condition is a density Gaussian of width w
! ------------------------------------------------------------------------------
SUBROUTINE initial_condition (n, w, rho)
  implicit none
  integer:: n, ix, iy
  real:: w, x, y, rho(n,n)
  !.............................................................................
  do iy=1,n
  do ix=1,n
    x = ix-(n+1)*0.5
    y = iy-(n+1)*0.5
    rho(ix,iy) = exp(-(x**2+y**2)/w**2)
  end do
  end do
END SUBROUTINE initial_condition
! ------------------------------------------------------------------------------
! Derivatives evaluated with central differences
! ------------------------------------------------------------------------------
SUBROUTINE ddx (n, f, d)
  implicit none
  integer:: n
  real, dimension(n,n):: f, d
  !.............................................................................
  d(2:n-1,:) = (f(3:n,:)-f(1:n-2,:))*0.5
  d(1    ,:) = (f(2  ,:)-f(  n  ,:))*0.5
  d(  n  ,:) = (f(1  ,:)-f(  n-1,:))*0.5
END SUBROUTINE ddx
SUBROUTINE ddy (n, f, d)
  implicit none
  integer:: n
  real, dimension(n,n):: f, d
  !.............................................................................
  d(:,2:n-1) = (f(:,3:n)-f(:,1:n-2))*0.5
  d(:,1    ) = (f(:,2  )-f(:,  n  ))*0.5
  d(:,  n  ) = (f(:,1  )-f(:,  n-1))*0.5
END SUBROUTINE ddy
! ------------------------------------------------------------------------------
! Simplest possible time step -- first order in time
! ------------------------------------------------------------------------------
SUBROUTINE time_step (n, courant, dx, rho, ux, uy)
  implicit none
  integer:: n
  real:: courant, ux, uy, dx, dt, umax
  real, dimension(n,n):: rho, dfdx, dfdy
  !.............................................................................
  umax = sqrt(ux**2+uy**2)
  dt = courant*dx/umax
  call ddx (n, rho*ux, dfdx)
  call ddy (n, rho*uy, dfdy)
  rho = rho - (dt/dx)*(dfdx + dfdy)
END SUBROUTINE time_step
! ------------------------------------------------------------------------------
! This test should just move the Gaussian, without changing the shape, and w/o
! adding any erroneous new features.  Run it first for 100 steps, then 300, then
! 1000 steps.  Vary the Courant number (the fraction of a mesh the profiles is
! moved per step), and try to move the shape as far as possible (disregarding
! the number of time steps needed).
! ------------------------------------------------------------------------------
PROGRAM advection_test
  implicit none
  integer:: n=100, n_step=200, i
  real:: courant=0.2, w=10., ux=1.0, uy=0.5, dx=1.0
  real, allocatable:: rho(:,:)
  namelist /input/ n, n_step, courant, w, ux, uy
  !.............................................................................
  read (*,input); write (*,input)
  allocate (rho(n,n))
  call initial_condition (n, w, rho)
  open (1, file='classical.dat', form='unformatted', status='unknown')
  write (1) n
  do i=1,n_step
    call time_step (n, courant, dx, rho, ux, uy)
    if (modulo(i,10)==0) print*,i,minval(rho),maxval(rho)
    write (1) i, rho
  end do
  close (1)
END

