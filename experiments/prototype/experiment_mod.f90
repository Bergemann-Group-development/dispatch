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
  USE item_mod
  USE task_mod
  USE task_list_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Task structure, overriding init, update, and output
  !-----------------------------------------------------------------------------
  type, public, extends(tasks_t):: experiment_t
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
  end type
  logical:: first_time=.true.
  integer:: dims(3)
  integer:: id=0
contains
!> =============================================================================
!> This is where to compute initial conditions
!> =============================================================================
subroutine init (self)
  class(experiment_t):: self
  !.............................................................................
  call trace_begin ('experiment_t%init')
  call self%item_t%init                       ! generic item initialization
  call trace_end
  first_time = .false.
end subroutine init
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
!> Any special handling of a time step update, in addition to the standard one
!> =============================================================================
subroutine update (self)
  class(experiment_t):: self
  !.............................................................................
  call trace_begin ('experiment_t%update')
  call self%task_t%update
  call trace_end
end subroutine
!> =============================================================================
!> Initialize the experiment
!> =============================================================================
SUBROUTINE init_experiment (list, name)
  class(experiment_list_t):: list
  character(len=*):: name
  type(experiment_t), pointer:: task
  class(item_t), pointer:: item
  integer:: i, n=5
  namelist /experiment_params/ dims
  !.............................................................................
  call trace_begin ('experiment_t%init_experiment')
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
  call list%list_t%init (name)                  ! generic initialization
  do i=1,n
    allocate (task)                             ! allocate task structure
    call task%init                              ! compute initial conditions
    item => task
    call list%append (item)                     ! append to the list
  end do
  !-----------------------------------------------------------------------------
  ! task_list speficic initialization 
  !-----------------------------------------------------------------------------
  call list%task_list_t%init (name)             ! specific initialization
  call list%print                               ! print the list (if verbose)
  call trace_end
END SUBROUTINE
END MODULE
