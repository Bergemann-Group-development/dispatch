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
    !!!procedure:: has_finished
  end type
  logical, save:: first_time = .true.

CONTAINS

!> =============================================================================
!> Initialise the experiment data type.
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  real:: toler=5e-4
  namelist / experiment_params / toler
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

  call self%output
  call timestep%update(self%id, self%iit, self%dt, self%mem, self%mesh)
  call self%dnload  (only=self%idx%phi)
  call pboundaries%conditions(self, self%it)
  call self%gravity (error, only=self%idx%phi)

  call trace%end
END SUBROUTINE update

!> =============================================================================
!> Check to see if the current task has converged, defined here as meaning the
!> change in the error/residual since the previous iteration falls below a
!> user-defined residual.
!> =============================================================================
!LOGICAL FUNCTION has_finished (self)
!  class (experiment_t):: self
!  real:: rd
!  !.............................................................................
!  has_finished = .false.
!  if (self%error(self%iit(self%nt-1)) == 0.0) return
!  rd = abs(self%error(self%iit(self%nt)) - self%error(self%iit(self%nt-1))) &
!     / self%error(self%iit(self%nt-1))
!  if (rd < self%toler) has_finished = .true.

!END FUNCTION has_finished

END MODULE experiment_mod
