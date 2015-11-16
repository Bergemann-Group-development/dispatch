!> =============================================================================
!> $Id$
!> This version uses a simple form of object-orientation, where everything
!> of concern to the experiments is stored inside the object(s), and which
!> illustrates the very significant advantage of such arrangements, in that
!> one can easily create multiple instances, which can optionally be
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
  private
  type, public:: experiment_t
    integer:: n, ns, unit
    real:: w, courant, dx, dt, velocity(3), position(3)=0.0, time=0.0
    real, dimension(:,:,:), pointer:: rho, rho0
    type(experiment_t), pointer:: next => null()
    character(len=32):: filename
    logical:: do_output=.false.
    contains
    procedure:: init
    procedure:: gaussian
    procedure:: update
    procedure:: output
    procedure:: print
    procedure:: run
  end type
CONTAINS
!> =============================================================================
!> Initialize the experiment, reading parameters from an input file (with one
!> namelist input required per task).  Note that procedures that are part of
!> a derived type are called with an extra first argument, which is the
!> relevant instance of the derived data type (object).  This way, the
!> procedure automatically has access to all required information.  Here,
!> the argument is called "self", to indicate this.
!> =============================================================================
SUBROUTINE init (self, m)
  class(experiment_t):: self
  integer:: i, j, k, m
  integer:: n=100, n_step=200
  real:: dx=1.0, courant=0.2, w=10, ux=1.0, uy=0.5
  logical:: do_output
  namelist /input/ n, n_step, courant, w, ux, uy, filename, do_output
  character(len=32):: filename
  !.............................................................................
  read (*,input)
  self%ns = n_step
  self%n  = n
  self%dx = dx
  self%w  = w
  self%courant   = courant
  self%velocity  = [ux,uy,0.0]
  self%do_output = do_output

  write (filename,'(a,i2.2,a)') 'objects_', m, '.dat'
  self%unit = 10+m
  if (self%do_output) then
    self%filename = filename
    print*, self%unit, filename
    open (self%unit, file=filename, form='unformatted', status='unknown')
    write (self%unit) n
  end if
  !write (*,input)
  call self%gaussian (self%rho)
END SUBROUTINE
!> =============================================================================
!> Construct a Gaussian centered on self%position, for initial conditions and
!> for comparisons.
!> =============================================================================
SUBROUTINE gaussian (self, data)
  class(experiment_t):: self
  real, dimension(:,:,:), pointer:: data
  integer:: i, j, k
  real:: x, y, position(3)
  !.............................................................................
  allocate (data (self%n, self%n, self%n))
  do k=1,self%n
  do j=1,self%n
  do i=1,self%n
    x = i-(self%n+1)*0.5 - self%position(1)
    y = j-(self%n+1)*0.5 - self%position(2)
    data(i,j,k) = exp(-(x**2+y**2)/self%w**2)
  end do
  end do
  end do
END SUBROUTINE
!> =============================================================================
!> Derivatives evaluated with central differences.
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
!> Simple time stepping -- 1st order:    f[n+1] = f[n] + dt*dfdt[n]
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  real:: umax
  !.............................................................................
  umax = sqrt(sum(self%velocity**2))
  self%dt = self%courant*self%dx/umax
  self%rho = self%rho - (self%dt/self%dx)* ( &
            ddx(self%n, self%rho*self%velocity(1)) &
          + ddy(self%n, self%rho*self%velocity(2)) )
  self%time = self%time + self%dt
  self%position = self%position + self%velocity*self%dt
END SUBROUTINE
!> =============================================================================
!> Print summary, with error measure
!> =============================================================================
SUBROUTINE print (self)
  class(experiment_t):: self
  integer:: istep, omp_get_thread_num
  !.............................................................................
  call self%gaussian (self%rho0)
  print'(2(a,i5,3x),a,f6.3,3x,a,1p,e10.2)', &
    'thread:', omp_get_thread_num(), &
    'n_step:', self%ns, &
    'Courant:', self%courant, &
    'rms:', sqrt(sum((self%rho-self%rho0)**2))/sqrt(sum(self%rho**2))
END SUBROUTINE
!> =============================================================================
!> Binary I/O to objects_nn.dat files (only if do_output=.true.)
!> =============================================================================
SUBROUTINE output (self)
  class(experiment_t):: self
  integer:: istep, omp_get_thread_num
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
  do istep=1,self%ns                            ! for self%ns time steps
    call self%update                            ! update
  end do
  call self%print                               ! print summary
  deallocate (self%rho, self%rho0)              ! deallocate memory
  if (self%do_output) then                      ! optionally ...
    call self%output                            ! write binary output
    close (self%unit)                           ! close binary file
  end if
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
  USE mpi_timer
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
  ! equally over the number of available threads (the default). This can be
  ! achieved with options to the !$omp parallel directive, choosing
  ! schedule(dynamic,1), or schedule(guided).
  !-----------------------------------------------------------------------------
  call tic                                      ! start Matlab style timer
  !$omp parallel do schedule(dynamic,1)
  do k=1,n
    call experiments(k)%run
  end do
  call toc ('run', n*experiments(1)%n**3)       ! stop MatLab style timer
END PROGRAM
