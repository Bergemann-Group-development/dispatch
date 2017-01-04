!> =============================================================================
!> Generic module, used by the task_mod module to initialize and update the
!> experiment
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE bits_mod
  USE rmhd_mod
  USE initial_mod
  USE timestep_mod
  USE task_mod
  USE rt_mod
  USE rmhd_mod
  implicit none
  private
  type, public, extends(rmhd_t):: experiment_t
    type(initial_t):: initial
  contains
    procedure init
    procedure update
    procedure boundary_conditions
  end type
  logical, save:: first_time=.true.
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  call trace%begin('experiment_mod::init', itimer=itimer)
  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  call self%initial%init
  self%nv  = self%initial%nv
  self%mhd = self%initial%mhd
  call self%rmhd_t%init
  self%i_d =1
  self%i_px=2
  self%i_py=3
  self%i_pz=4
  self%i_e =5
  self%i_tt=6
  self%i_bx=7
  self%i_by=8
  self%i_bz=9
  call self%set (bits%static)
  !----------------------------------------------------------------------------
  ! Initial conditions
  !----------------------------------------------------------------------------
  call self%initial%condition (self%mesh, self%mem(:,:,:,:,1,1))
  if (self%root_grid.or.self%is_set(bits%root_grid)) then
    call self%periodic_grid                     ! periodic root grid
  else
    call self%extrapolate_guards                ! time-extrapolate guard cells
  end if
  self%rt%mesh(2)%dsup => self%mesh(2)%dsup                             ! mesh spacing
  self%time = self%initial%time
  self%t(self%it) = self%initial%time
  call trace%end (itimer=itimer)
END SUBROUTINE init

!> =============================================================================
!> Compute time derivatives, and take a time step.  Note that this procedure is
!> only called if/when the task addet to the task list was an experiment_t task
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  type(timestep_t):: timestep
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update')
  if (first_time) then
    first_time = .false.
    call self%initial%cleanup
  end if
  call self%rmhd_t%update
  call self%boundary_conditions
  call self%courant_condition
  call self%output
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call trace%end
END SUBROUTINE update

!> =============================================================================
!> Apply boundary conditions
!> =============================================================================
SUBROUTINE boundary_conditions (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%boundary_conditions')
  associate (my => self%mesh(2)%gn)
  self%mem(:,  01:06,:,:,:,2) = 0.0
  self%mem(:,my-4:my,:,:,:,2) = 0.0
  end associate
  call trace%end
END SUBROUTINE boundary_conditions

END MODULE
