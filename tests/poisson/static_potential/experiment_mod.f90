!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE bits_mod
  USE solver_mod
  USE initial_mod
  USE pboundary_mod
  USE timestep_mod
  implicit none
  private

  type, public, extends(solver_t):: experiment_t
    real, dimension(:), allocatable:: error
  contains
    procedure:: init
    procedure:: update
    procedure:: has_finished
  end type
  logical, save:: first_time = .true.
  real, save:: ema_error=0.0, aepsilon=1e-3
CONTAINS

!> =============================================================================
!> Initialise the experiment data type.
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  real:: init_ema_error=0.0
  namelist / experiment_params / init_ema_error, aepsilon
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%init')

  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*, experiment_params)
    call timestep%init()
  end if
  !$omp end critical (input_cr)

  self%nv = 1
  self%nw = 2
  self%periodic = .false.
  
  call self%gpatch_t%init()
  call self%clear(bits%root_grid)
  self%kind = 'poisson'
  self%dtime = 32.0 * minval(self%ds)
  self%d0 = 0.0

  ! allocate error storage array
  allocate(self%error(self%nt))
  self%error(:) = 0.0
  ! initialise exponetial moving average (assumes 4piG = 1).
  if (init_ema_error == 0.0) then
    ema_error = max(ema_error, sum(self%mem(:,:,:,self%idx%d,1,1)))
  else
    ema_error = init_ema_error
  end if

  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Update procedure for this experiment.
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  real:: error
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update')
  
  call pboundaries%conditions(self, self%it)
  call self%gravity(error)
  self%error(self%it) = error
  ema_error = ema_error + 0.025 * (error - ema_error)

  print "('ID: ',i4.1,' EMA error = ',1pe11.3,' local error = ',1pe11.4,'; steps taken: ',i4.1)",&
    self%id, ema_error, error, self%istep

  call self%output
  call timestep%update(self%id, self%iit, self%dt, self%mem, self%mesh)

  call trace%end
END SUBROUTINE update

!> =============================================================================
!> Check to see if the test has converged.
!> =============================================================================
LOGICAL FUNCTION has_finished (self)
  class (experiment_t):: self
  !.............................................................................
  has_finished = .false.
  if (ema_error < aepsilon) then
    has_finished = .true.
    self%out_next = self%time
    call self%output
  end if
END FUNCTION has_finished

END MODULE experiment_mod
