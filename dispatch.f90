!> =============================================================================
!> $Id$
!> =============================================================================
PROGRAM dispatch
  USE trace_mod
  USE task_list_mod
  USE mesg_list_mod
  USE experiment_mod
  USE mpi_mod
  implicit none
  type(task_list_t):: task_list, ready_queue
  type(mesg_list_t):: mesg_list, recv_list
  type(experiment_t):: experiment
  logical:: running = .true.
  !....................................................................
  call trace_begin ('dispatch')
  call mpi%init                                     ! start up MPI
  call io%init                                      ! start up I/O
  call task_list%init ('task_list')                 ! initialize tasks_list
  call ready_queue%init ('ready_queue')             ! initialize ready_queue
  !$omp parallel
  call trace('init_experiment')
  call init_experiment (task_list)                  ! populate lists with tasks
  do while (task_list%size > 0)                     ! as long as there are tasks
    call mesg_list%probe                            ! probe for new messages
    call mesg_list%receive (recv_list)              ! collect received messages
    call recv_list%unpack (task_list)               ! unpack received messages
    call ready_queue%head%item%execute              ! execute the oldest task
  end do
  !$omp end parallel
  call mpi%end                                      ! finalize
END PROGRAM
