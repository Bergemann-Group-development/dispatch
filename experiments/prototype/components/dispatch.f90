!> =============================================================================
!> The definition module experiment_mod.o in the specific experiment directory,
!> needs to have an init method, and needs to have a task list; that's all.
!> 
!> The experiment_mod layer is only needed to anonymize the specifics that go
!> into the creation and population of the experiment%task_list.
!> =============================================================================
PROGRAM dispatch
  USE experiment_mod
  USE mpi_mod
  USE io_mod
  implicit none
  type(experiment_t), pointer:: experiment
  !.............................................................................
  call mpi%init
  call io%init
  allocate (experiment)
  call experiment%init
  do while (experiment%task_list%n > 0)
    call experiment%task_list%update
  end do
  call mpi%end
END PROGRAM
