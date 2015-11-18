!> =============================================================================
!> $Id$
!> Prototype experiment, showing the steps and components needed:
!> 1) a "task_list%init" procedure, which adds new tasks to the task_list
!> 2) an "init" proc, which is executed by the task_list%init
!> 3) an "update" procedure, which is executed by task_list%update
!> 3) an "output" procedure, which does whatever is desired
!> =============================================================================

MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE mpi_coord_mod
  USE task_mod
  USE tasks_mod
  USE patch_mod
  USE task_list_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Task structure, overriding init, update, and output
  !-----------------------------------------------------------------------------
  type, public, extends(patch_t):: experiment_t
  contains
    procedure:: init
    procedure:: update
    procedure:: output
  end type
  !-----------------------------------------------------------------------------
  ! Task list structure, overriding only init
  !-----------------------------------------------------------------------------
  type, public, extends(task_list_t):: experiment_list_t
  contains
    procedure:: init => init_experiment
    procedure:: output => output_experiment
  end type
  logical:: first_time=.true.
  integer:: dims(3), n=64
  real(8):: end_time=1.0, grace=0.0
contains
!> =============================================================================
!> This is where to compute initial conditions
!> =============================================================================
subroutine init (self)
  class(experiment_t):: self
  !.............................................................................
  call trace_begin ('experiment_t%init')
  call self%task_t%init                       ! generic task initialization
  self%dtime = 0.01+(0.09*self%id)/n
  self%grace = grace*self%dtime
  self%end_time = end_time
  call trace_end
  first_time = .false.
end subroutine init
!> =============================================================================
!> Any special handling of a time step update, in addition to the standard one
!> =============================================================================
subroutine update (self)
  class(experiment_t):: self
  !.............................................................................
  call trace_begin ('experiment_t%update')
  call self%task_t%update
  print*, self%id, self%time, self%position
  call trace_end
end subroutine
!> =============================================================================
!> Output file handling here
!> =============================================================================
subroutine output (self)
  class(experiment_t):: self
  !.............................................................................
  call trace_begin ('experiment_t%output')
  call self%print                             ! just a generic print, for now
  call trace_end
end subroutine

!> =============================================================================
!> Initialize the experiment
!> =============================================================================
SUBROUTINE init_experiment (self, name)
  class(experiment_list_t):: self
  character(len=*):: name
  type(experiment_t), pointer:: exper
  class(task_t), pointer:: task
  integer:: i
  real(8):: size=1.0, patch_size=0.25, position(3)
  namelist /experiment_params/ end_time, grace, dims, size, patch_size
  !.............................................................................
  call trace_begin ('experiment_list_t%init_experiment')
  !-----------------------------------------------------------------------------
  ! Read parameters
  !-----------------------------------------------------------------------------
  rewind (io%input)
  read (io%input, experiment_params)
  if (io%master.and.first_time) write (*, experiment_params)
  call mpi_coord%init (dims=dims)
  !-----------------------------------------------------------------------------
  ! Initialize a task list, and add a few tasks to it
  !-----------------------------------------------------------------------------
  self%name = name
  self%size = size
  position = 0.5*patch_size
  do i=1,n
    allocate (exper)                            ! allocate task structure
    exper%size = patch_size
    exper%position = position
    exper%id = i
    task => exper                               ! class pointer
    call self%append (task)                     ! append to the list
    position(1) = position(1) + patch_size
    if (position(1) > size) then
      position(1) = 0.5*patch_size
      position(2) = position(2) + patch_size
    end if
    if (position(2) > size) then
      position(2) = 0.5*patch_size
      position(3) = position(3) + patch_size
    end if
  end do
  call self%task_list_t%init (name)             ! initial_conditions and queue
  !-----------------------------------------------------------------------------
  ! task_list speficic initialization
  !-----------------------------------------------------------------------------
  call self%print                               ! print the list (if verbose)
  call trace_end
END SUBROUTINE init_experiment

SUBROUTINE output_experiment (self)
  class(experiment_list_t):: self
  print '(a,3(5x,a,i5))', 'output_experiment:', &
    'n_ready:', self%n_ready, &
    'n_fail:' , self%n_fail, &
    'n_nbor:' , self%n_nbor
END SUBROUTINE output_experiment

END MODULE
