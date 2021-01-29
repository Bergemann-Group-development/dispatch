!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE bits_mod
  USE solver_mod
  USE initial_mod
  USE timestep_mod
  USE pboundary_mod
  implicit none
  private

  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
    procedure:: update
  end type
  logical, save:: first_time = .true.
CONTAINS

!===============================================================================
!> Initialise the experiment data type.
!===============================================================================
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
  
  call self%idx%init (ienergy=2, .false.)
  call self%solver_t%init()
  self%level = 1
  self%kind = 'poisson'
  self%dtime = 1.0
  self%d0 = 0.0

  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Update procedure for this experiment.
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  real:: error
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update')
  call self%solver_t%update
  call pboundaries%conditions(self, self%it)
  call timestep%update(self%id, self%iit, self%dt, self%mem, self%mesh, self%lock)
  call self%counter_update
  call self%output
  call trace%end()
END SUBROUTINE update

END MODULE experiment_mod
