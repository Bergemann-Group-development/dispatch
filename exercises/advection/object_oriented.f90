!> =============================================================================
!> $Id$
!> This version uses a simple form of object-orientation, where everything
!> of concern to the experiments is stored inside the object(s), and which
!> illustrates the very significant advantage of such arrangements, in that
!> one can easily create multiple instances, which can !> optionally be
!> executed in parallel.
!>
!> In this example, which is constructed to just illustrate the principles,
!> there is no major advantage, compared to running the same code repeatedly
!> (or parallelizing with MPI), but soon we will use the same type of instancing
!> to create processes (tasks) that actually collaborate to advance the same
!> experiment, thus parallelizing with OpenMP on each node, without having to
!> parallelize each task over multiple OMP threads,
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
!> Initialize the experiment, reading parameters from an input file (with one
!> namelist input required per task).  Note that procedures that are part of
!> a derived type are called with an extra first argument, which is the
!> relevant instance of the derived data type (object).  This way, the
!> procedure automatically has access to all required information.  Here,
!> the argument is called "self", to indicate this.
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
!> Derivatives evaluated with central differences.  Note that these array-
!> valued functions are not made part of the experiment object (experiment_t),
!> they are just local procedures inside the module. In this isolated example,
!> one could let everithing except "init" and "run" be similarly local, but for
!> the purpose of illustrationm "update" and "output" are made parts of the
!> experiment_t object, as they need to be in more complex and realistic
!> applications.
!> =============================================================================
FUNCTION ddx (n, f)
  integer:: n
  real, dimension(n,n,n):: f, ddx
  !.............................................................................
  ddx(2:n-1,:,:) = (f(3:n,:,:)-f(1:n-2,:,:))*0.5
  ddx(1    ,:,:) = (f(2  ,:,:)-f(  n  ,:,:))*0.5
  ddx(  n  ,:,:) = (f(1  ,:,:)-f(  n-1,:,:))*0.5
END FUNCTION ddx
FUNCTION ddy (n, f)
  integer:: n
  real, dimension(n,n,n):: f, ddy
  !.............................................................................
  ddy(:,2:n-1,:) = (f(:,3:n,:)-f(:,1:n-2,:))*0.5
  ddy(:,1    ,:) = (f(:,2  ,:)-f(:,  n  ,:))*0.5
  ddy(:,  n  ,:) = (f(:,1  ,:)-f(:,  n-1,:))*0.5
END FUNCTION ddy
!> =============================================================================
!> Simple time stepping -- 1st order
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
!> I/O routines
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
!> adding any erroneous new features.  Run it with varying number of steps, and
!> varying Courant number (the fraction of a mesh the profiles is moved per step),
!>  and try to move the shape as far as possible (disregarding the number of
!> time steps needed).
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

  !-----------------------------------------------------------------------------
  ! Initialize the experiments, one after the other, reading the namelists in
  ! serial order.
  !-----------------------------------------------------------------------------
  do k=1,n
    call experiments(k)%init(k)
  end do

  !-----------------------------------------------------------------------------
  ! Run the experiments in parallel, on the available OMP threads.  If the
  ! number of experiments is much larger than the number of threads one can
  ! expect good load balance if the threads are assigned to experiments
  ! dynamically, rather than by just dividing the number of experiments
  ! equally over the number of available threads (the dedault). This can be
  ! achieved with options to the !$omp parallel directive, choosing
  ! schedule(dynamic,1), or schedule(guided).
  !-----------------------------------------------------------------------------
  !$omp parallel do schedule(dynamic,1)
  do k=1,n
    call experiments(k)%run
  end do

END PROGRAM

