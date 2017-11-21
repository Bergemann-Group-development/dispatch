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
  contains
    procedure:: init
    procedure:: update
  end type
  logical, save:: first_time = .true.
CONTAINS

!> =============================================================================
!> Initialise the experiment data type.
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%init')

  if (first_time) then
    call timestep%init
    first_time = .false.
  end if

  self%nv = 1
  self%nw = 2
  self%periodic = .false.
  
  call self%gpatch_t%init()
  call self%clear(bits%root_grid)
  self%kind = 'poisson'
  self%dtime = 32.0 * minval(self%ds)
  self%d0 = 0.0

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
  call self%gravity

  call self%output
  call timestep%update(self%id, self%iit, self%dt, self%mem, self%mesh)

  call self%rotate
  call self%dnload (only=self%idx%d)
  call self%dnload (only=self%idx%phi)
  call pboundaries%conditions(self, self%it)
  call self%gravity

  call trace%end
END SUBROUTINE update

END MODULE experiment_mod
