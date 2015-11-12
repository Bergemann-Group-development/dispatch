!> =============================================================================
!> $Id$
!> =============================================================================
PROGRAM dispatch
  USE list_mod
  USE task_mod
  USE experiment_mod
  implicit none
  type(list_t), pointer:: task_list, ready_queue
  type(list_t), pointer:: send_list, recv_list, mesg_list
  type(experiment_t):: experiment
  logical:: running = .true.
  !....................................................................
  call experiment%init (ready_queue)        ! populate ready_queue with tasks
  do while (task_list%size > 0)             ! as long as there are tasks
    call recv_list%check                    ! check incoming MPI
    call recv_list%move_to (mesg_list)      ! move completed messages
    call mesg_list%unpack_to (task_list)    ! unpack completed messages
    call ready_queue%head%execute           ! execute the oldest task
  end do
END PROGRAM
