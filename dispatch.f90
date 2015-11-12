!> =============================================================================
!> $Id$
!> =============================================================================
PROGRAM dispatch
  USE list_mod
  !USE task_mod
  USE mesg_list_mod
  USE experiment_mod
  implicit none
  type(list_t), pointer:: task_list, ready_queue
  type(mesg_list_t), pointer:: mesg_list, recv_list
  type(experiment_t):: experiment
  logical:: running = .true.
  !....................................................................
  call init_mpi                                 ! start up MPI
  !$omp parallel
  call experiment%init (task_list, ready_queue) ! populate lists with tasks
  do while (task_list%size > 0)                 ! as long as there are tasks
    call mesg_list%probe                        ! probe for new messages
    call mesg_list%receive (recv_list)          ! collect received messages
    call recv_list%unpack (task_list)           ! unpack received messages
    call ready_queue%head%item%execute          ! execute the oldest task
  end do
  !$omp end parallel
END PROGRAM
