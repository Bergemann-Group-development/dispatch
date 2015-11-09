!> =============================================================================
!> $Id$
!> =============================================================================
MODULE experiment
  implicit none
  type:: experiment_t
    integer:: n, n_step
    real:: w, courant, dx, dt, ux, uy
    real, dimension(:,:), pointer:: rho
    contains
    procedure:: init
    procedure:: update
    procedure:: output
  end type
PRIVATE
PUBLIC experiment_t
CONTAINS
!> =============================================================================
!> =============================================================================
SUBROUTINE init (self, n_stp)
  class(experiment_t):: self
  integer:: ix, iy, n=100, n_step=200, n_stp
  real:: dx=1.0, x, y, courant=0.2, w=10, ux=1.0, uy=0.5
  namelist /input/ n, n_step, courant, w, ux, uy
  !.............................................................................
  read(*,input); write(*,input)
  self%dx=dx; self%ux=ux; self%uy=uy; self%courant=courant; self%w=w; self%n=n
  self%n_step=n_step; self%n=n; n_stp=n_step
  allocate (self%rho(self%n,self%n))
  do iy=1,self%n
  do ix=1,self%n
    x = ix-(self%n+1)*0.5
    y = iy-(self%n+1)*0.5
    self%rho(ix,iy) = exp(-(x**2+y**2)/self%w**2)
  end do
  end do
  open (1, file='advection.dat', form='unformatted', status='unknown')
  write (1) n
END SUBROUTINE init
!> =============================================================================
!> Derivatives evaluated with central differences
!> =============================================================================
FUNCTION ddx (n, f)
  integer:: n
  real, dimension(n,n):: f, ddx
  !.............................................................................
  ddx = (cshift(f,dim=1,shift=1)-cshift(f,dim=1,shift=-1))*0.5
END FUNCTION ddx
FUNCTION ddy (n, f)
  integer:: n
  real, dimension(n,n):: f, ddy
  !.............................................................................
  ddy = (cshift(f,dim=2,shift=1)-cshift(f,dim=2,shift=-1))*0.5
END FUNCTION ddy
!> =============================================================================
!> Derivatives evaluated with central differences
!> =============================================================================
SUBROUTINE update (self, i)
  class(experiment_t):: self
  integer:: i
  real:: umax
  !.............................................................................
  umax = sqrt(self%ux**2+self%uy**2)
  self%dt = self%courant*self%dx/umax
  self%rho = self%rho - (self%dt/self%dx)* &
             (ddx(self%n, self%rho*self%ux) + ddy(self%n, self%rho*self%uy))
  if (modulo(i,10)==0) print*,i,minval(self%rho),maxval(self%rho)
END SUBROUTINE update
!> =============================================================================
!> I/O routine
!> =============================================================================
SUBROUTINE output (self, i)
  class(experiment_t):: self
  integer:: i
  !.............................................................................
  write (1) i, self%rho
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
  implicit none
  class(experiment_t), allocatable:: exp1, exp2
  integer:: n_step, i
  !.............................................................................
  allocate (exp1, exp2)
  call exp1%init(n_step)
  do i=1,n_step
    call exp1%update(i)
    call exp1%output(i)
  end do
END PROGRAM

