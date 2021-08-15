!===============================================================================
!> Generic DISPATCH main program.  A selection of components may be added; the
!> first one that defines a task list defines the job.
!===============================================================================
PROGRAM dispatch
  USE setup_mod
  USE restart_mod
  USE dispatcher_mod
  USE cartesian_vb_mod
  USE task_list_mod
  type(cartesian_vb_t):: cartesian              ! Use Cartesian patch arrangement
  class(task_list_t), pointer:: task_list => null()
  !.............................................................................
  call setup%init                               ! Generic setup (MPI, I/O, ...)
  call dispatcher%init                          ! Initialize the dispatcher
  call restart%init (task_list)                 ! Possibly initialize task_list
  call cartesian%init (task_list)               ! .. or else here ...
  call dispatcher%execute (task_list)           ! Run dispatcher on the task_list
  call setup%end                                ! End setup
!===============================================================================
END PROGRAM dispatch
