!===============================================================================
!> Generic DISPATCH main program
!===============================================================================
PROGRAM dispatch
  USE io_mod
  USE mpi_mod
  USE scaling_mod
  USE timer_mod
  USE download_mod
  USE dispatcher_mod
  USE cartesian_mod
  USE task_list_mod
  implicit none
  type(cartesian_t):: cartesian
  type(task_list_t), pointer:: task_list
  !.............................................................................
  call mpi%init                                 ! MPI startup
  call io%init                                  ! I/O system
  call scaling%init                             ! Scaling to code units
  call download%init                            ! Downloading and interpolation parameters
  call dispatcher%init                          ! Initialize the dispatcher
  call cartesian%init
  task_list => cartesian%task_list
  call setup (task_list)                        ! Populate task list
  call dispatcher%execute (task_list)           ! Run the dispatcher
  call mpi%end
END PROGRAM dispatch

!===============================================================================
!> Add a few sink tasks
!===============================================================================
SUBROUTINE setup (task_list)
  USE task_list_mod
  USE cartesian_mod
  USE sinks_mod
  USE task_mod
  implicit none
  type(task_list_t):: task_list
  !.............................................................................
  class(sinks_t), pointer:: sink
  class(task_t), pointer:: task
  integer:: i, n_sink=1
  !-----------------------------------------------------------------------------
  do i=1,n_sink
    allocate (sink)
    call sink%init
    sink%position = (i-(n_sink+1)/2.0)*0.25
    task => sink
    call task_list%append_task (task)
  end do
END SUBROUTINE setup
