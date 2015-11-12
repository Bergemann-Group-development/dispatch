!> =============================================================================
!> $Id$
!> =============================================================================
PROGRAM dispatch
  USE list_mod
  USE task_mod
  USE mesg_list_mod
  USE experiment_mod
  implicit none
  type(list_t), pointer:: task_list, ready_queue
  type(mesg_list_t), pointer:: recv_list, mesg_list
  type(experiment_t):: experiment
  logical:: running = .true.
  !....................................................................
  call experiment%init (task_list, ready_queue) ! populate lists with tasks
  do while (task_list%size > 0)                 ! as long as there are tasks
    call recv_list%probe                        ! check incoming MPI
    call mesg_list%receive (recv_list)          ! move completed messages
    call mesg_list%unpack (task_list)           ! unpack completed messages
    call ready_queue%head%item%execute          ! execute the oldest task
  end do
END PROGRAM
