!> =============================================================================
!> $Id$
!> =============================================================================
PROGRAM dispatch
!  USE list_mod
  USE task_list_mod
  USE mesg_list_mod
  USE experiment_mod
  implicit none
  type(task_list_t), pointer:: task_list, ready_queue
  type(mesg_list_t), pointer:: mesg_list, recv_list
  type(experiment_t):: experiment
  logical:: running = .true.
  !....................................................................
  call init_mpi                                     ! start up MPI
  !$omp parallel
  call task_list%init ('task_list')                 ! initialize tasks_list
  call ready_queue%init ('ready_queue')             ! initialize ready_queue
  call init_experiment (task_list)                  ! populate lists with tasks
  do while (task_list%size > 0)                     ! as long as there are tasks
    call mesg_list%probe                            ! probe for new messages
    call mesg_list%receive (recv_list)              ! collect received messages
    call recv_list%unpack (task_list)               ! unpack received messages
    call ready_queue%head%item%execute              ! execute the oldest task
  end do
  !$omp end parallel
END PROGRAM
