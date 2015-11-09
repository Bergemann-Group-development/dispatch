!> =============================================================================
!> $Id$
!> This version uses a Fortran module that isolates the experiment data from
!> the main program, which only is aware of two number: the current step number,
!> and the total number of steps (n_stp). The module also transfers information
!> about the calling sequences etc, which allows the use of array valued fns.
!> =============================================================================
MODULE experiment
  implicit none
  integer:: n=100, n_step=100
  real:: w=10, courant=0.2, dx=1.0, dt, ux=1.0, uy=0.5
  real, dimension(:,:), pointer:: rho
PRIVATE
PUBLIC init, update, output
CONTAINS
!> =============================================================================
!> =============================================================================
SUBROUTINE init(n_stp)
  integer:: ix, iy, n_stp
  real:: x, y
  namelist /input/ n, n_step, courant, w, ux, uy
  !.............................................................................
  read(*,input); write(*,input); n_stp=n_step
  allocate (rho(n,n))
  do iy=1,n
  do ix=1,n
    x = ix-(n+1)*0.5
    y = iy-(n+1)*0.5
    rho(ix,iy) = exp(-(x**2+y**2)/w**2)
  end do
  end do
  open (1, file='advection.dat', form='unformatted', status='unknown')
  write (1) n
END SUBROUTINE init
!> =============================================================================
!> Derivatives evaluated with central differences
!> =============================================================================
FUNCTION ddx(n,f)
  integer:: n
  real, dimension(n,n):: f, ddx
  !.............................................................................
  ddx = (cshift(f,dim=1,shift=1)-cshift(f,dim=1,shift=-1))*0.5/dx
END FUNCTION ddx
FUNCTION ddy(n,f)
  integer:: n
  real, dimension(n,n):: f, ddy
  !.............................................................................
  ddy = (cshift(f,dim=2,shift=1)-cshift(f,dim=2,shift=-1))*0.5/dx
END FUNCTION ddy
!> =============================================================================
!> Derivatives evaluated with central differences
!> =============================================================================
SUBROUTINE update(i)
  real:: dt, umax
  integer:: i
  !.............................................................................
  umax = sqrt(ux**2+uy**2)
  dt = courant*dx/umax
  rho = rho - dt*(ddx(n,rho*ux) + ddy(n,rho*uy))
  if (modulo(i,10)==0) print*,i,minval(rho),maxval(rho)
END SUBROUTINE update
!> =============================================================================
!> I/O routine
!> =============================================================================
SUBROUTINE output (i)
  integer:: i
  !.............................................................................
  write (1) i, rho
END SUBROUTINE output
END MODULE experiment

!> =============================================================================
!> This test should just move the Gaussian, without changing the shape, and w/o
!> adding any erroneous new features.  Run it first for 100 steps, then 300, then
!> 1000 steps.  Vary the Courant number (the fraction of a mesh the profiles is
!> moved per step), and try to move the shape as far as possible (disregarding
!> the number of time steps needed).
!> =============================================================================
PROGRAM advection_experiment
  USE experiment
  integer:: n_stp, i
  call init(n_stp)
  do i=1,n_stp
    call update(i)
    call output(i)
  end do
END PROGRAM

