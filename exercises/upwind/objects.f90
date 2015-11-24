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
  USE mpi_mod
  USE mpi_coords_mod
  USE mpi_send_mod
  USE mpi_reduce_mod
  USE derivatives_mod
  USE upwind_mod
  USE central_mod
  implicit none
  private
  type, public:: experiment_t
    integer:: n(3), ns, unit
    real:: w, courant, dx, dt, velocity(3), position(3)=0.0, time=0.0
    real, dimension(:,:,:), pointer:: rho, rho0
    type(experiment_t), pointer:: next => null()
    character(len=32):: filename
    logical:: do_output=.false.
    logical:: do_upwind=.false.
    type(central_t) :: central
    type(upwind_t)  :: upwind
    class(derivatives_t), pointer :: d
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
  class(experiment_t),target:: self
  integer:: i, j, k, m, rank
  integer:: n=63, n_step=200
  real:: dx=1.0, courant=0.2, w=10, ux=1.0, uy=0.5, uz=0.8, upwind=1.0
  logical:: do_output=.true., do_upwind=.false.
  namelist /input/ n, n_step, courant, w, ux, uy, uz, filename, do_output, &
                   do_upwind, upwind
  character(len=32):: filename
  !.............................................................................
  read (1,input)

  call mpi_coords%init (dims=[n,n,n])
  call mpi_coords%print

  self%ns = n_step
  self%n  = n/mpi_coords%dims
  self%position = self%n*(mpi_coords%coords-0.5*(mpi_coords%dims-1))
  self%dx = dx
  self%w  = w
  self%courant   = courant
  self%velocity  = [ux,uy,uz]
  self%do_upwind = do_upwind

  ! Select which differential operator type to use
  if (do_upwind) then
    self%d => self%upwind
    self%upwind%factor = upwind
    self%upwind%velocity => self%velocity ! pass information about advection velocity
  else
    self%d => self%central
  endif
  self%d%n = self%n
  self%d%ds = self%dx
  self%d%dt = self%dt
    
  self%do_output = do_output
  if (mpi%rank==0) print*,'                              positions'
  call mpi%barrier (delay=10)
  do rank=0,mpi%size-1
    if (rank==mpi%rank) &
      print*, mpi%rank, self%position
    call mpi%barrier (delay=1)
  end do

  write (filename,'(a,i2.2,a)') 'objects_', m, '.dat'
  self%unit = 10+m
  if (self%do_output) then
    self%filename = filename
    print*, self%unit, filename
    open (self%unit, file=filename, form='unformatted', status='unknown')
    write (self%unit) n
  end if
  if (mpi%rank==0) write(*,input)
  call self%gaussian (self%rho)
END SUBROUTINE init
!> =============================================================================
!> Construct a Gaussian centered on self%position, for initial conditions and
!> for comparisons.
!> =============================================================================
SUBROUTINE gaussian (self, data)
  class(experiment_t):: self
  real, dimension(:,:,:), pointer:: data
  integer:: i, j, k
  real:: x, y, z, position(3)
  !.............................................................................
  allocate (data (self%n(1), self%n(2), self%n(3)))
  do k=1,self%n(3)
  do j=1,self%n(2)
  do i=1,self%n(1)
    x = i-(self%n(1)+1)*0.5 - self%position(1)
    y = j-(self%n(2)+1)*0.5 - self%position(2)
    z = k-(self%n(3)+1)*0.5 - self%position(3)
    data(i,j,k) = exp(-(x**2+y**2+z**2)/self%w**2)
  end do
  end do
  end do
END SUBROUTINE gaussian
!> =============================================================================
!> Simple time stepping -- 1st order:    f[n+1] = f[n] + dt*dfdt[n]
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  real:: umax
  !.............................................................................
  umax = sqrt(sum(self%velocity**2))
  self%dt = self%courant*self%dx/umax
  self%d%dt = self%dt
  self%rho = self%rho - (self%dt/self%dx)* ( &
            self%d%dx(self%rho*self%velocity(1)) &
          + self%d%dy(self%rho*self%velocity(2)) &
          + self%d%dz(self%rho*self%velocity(3)) )
  self%time = self%time + self%dt
  self%position = self%position + self%velocity*self%dt
END SUBROUTINE update
!> =============================================================================
!> Print summary, with error measure
!> =============================================================================
SUBROUTINE print (self)
  class(experiment_t):: self
  integer:: istep
  real:: sum1, sum2
  !.............................................................................
  call self%gaussian (self%rho0)
  call mpi%barrier (delay=2)
  sum1 = sum((self%rho-self%rho0)**2)
  sum2 = sum(self%rho**2)
  call mpi_reduce%sum_real (sum1)
  call mpi_reduce%sum_real (sum2)
  if (mpi%rank==0) &
   print'(1x,a,i5,3x,a,f6.3,3x,a,1p,e10.2)', &
    'n_step:', self%ns, &
    'Courant:', self%courant, &
    'rms:', sqrt(sum1/sum2)
END SUBROUTINE
!> =============================================================================
!> Binary I/O to objects_nn.dat files (only if do_output=.true.)
!> =============================================================================
SUBROUTINE output (self)
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
  USE mpi_timer_mod
  USE mpi_mod
  USE mpi_coords_mod
  USE mpi_send_mod
  USE mpi_reduce_mod
  implicit none
  class(experiment_t), pointer, dimension(:) :: experiments
  integer:: n, k
  namelist /instances/ n
  !.............................................................................
  call mpi%init
  call mpi_send%init
  call mpi_reduce%init
  open (1, file='input.nml', form='formatted')
  read (1,instances)
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
  if (mpi%rank==0) &
    call tic                                      ! start Matlab style timer
  do k=1,n
    call experiments(k)%run
  end do
  call mpi%barrier (delay=10)
  if (mpi%rank==0) &
    call toc ('run', n*product(experiments(1)%n)) ! stop MatLab style timer
  call mpi%barrier (delay=20)
  call mpi%end
END PROGRAM
