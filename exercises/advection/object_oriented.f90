!> =============================================================================
!> $Id$
!> =============================================================================
MODULE experiment
  implicit none
  type:: experiment_t
    integer:: n, ns, unit
    real:: w, co, dx, dt, ux, uy
    real, dimension(:,:), pointer:: rho
    type(experiment_t), pointer:: next => null()
    character(len=32):: filename
    contains
    procedure:: init
    procedure:: update
    procedure:: output
    procedure:: run
  end type
PRIVATE
PUBLIC experiment_t
CONTAINS
!> =============================================================================
!> =============================================================================
SUBROUTINE init (self, k)
  class(experiment_t):: self
  integer:: i, j, k
  integer:: n=100, n_step=200
  real:: dx=1.0, courant=0.2, w=10, ux=1.0, uy=0.5
  real:: x, y
  namelist /input/ n, n_step, courant, w, ux, uy, filename
  character(len=32):: filename
  !.............................................................................
  read (*,input)
  self%ns = n_step
  self%n  = n
  self%dx = dx
  self%ux = ux
  self%uy = uy
  self%w  = w
  self%co = courant

  write (filename,'(a,i2.2,a)') 'object_oriented_', k, '.dat'
  self%unit = 10+k
  self%filename = filename
  open (self%unit, file=filename, form='unformatted', status='unknown')
  write (self%unit) n
  write (*,input)

  allocate (self%rho(self%n,self%n))
  do j=1,self%n
  do i=1,self%n
    x = i-(self%n+1)*0.5
    y = j-(self%n+1)*0.5
    self%rho(i,j) = exp(-(x**2+y**2)/self%w**2)
  end do
  end do
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
SUBROUTINE update (self)
  class(experiment_t):: self
  real:: umax
  !.............................................................................
  umax = sqrt(self%ux**2+self%uy**2)
  self%dt = self%co*self%dx/umax
  self%rho = self%rho - (self%dt/self%dx)* &
             (ddx(self%n, self%rho*self%ux) + ddy(self%n, self%rho*self%uy))
END SUBROUTINE
!> =============================================================================
!> I/O routine
!> =============================================================================
SUBROUTINE output (self, istep)
  class(experiment_t):: self
  integer:: istep
  !.............................................................................
  write (self%unit) istep, self%rho
END SUBROUTINE
!> =============================================================================
!> Run the experiment
!> =============================================================================
SUBROUTINE run (self)
  class(experiment_t):: self
  integer:: istep
  !.............................................................................
  do istep=1,self%ns
    call self%update
    call self%output(istep)
  end do
  print *, self%filename, minval(self%rho), maxval(self%rho)
  deallocate (self%rho)
  close (self%unit)
END SUBROUTINE
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
  class(experiment_t), pointer, dimension(:) :: experiments
  integer:: n, k
  namelist /instances/ n
  !.............................................................................
  read (*,instances)
  allocate (experiments(n))

  do k=1,n
    call experiments(k)%init(k)
  end do

  !$omp parallel do
  do k=1,n
    call experiments(k)%run
  end do

END PROGRAM

