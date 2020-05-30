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
  call cartesian%init                           ! Cartesian patches
  task_list => cartesian%task_list              ! Initial task list
  call setup (task_list)                        ! Add sinks
  call dispatcher%execute (task_list)           ! Run the dispatcher
  call mpi%end                                  ! MPI finalize
END PROGRAM dispatch

!===============================================================================
!> Add a few sink tasks
!===============================================================================
SUBROUTINE setup (task_list)
  USE task_list_mod
  USE cartesian_mod
  USE sinks_mod
  USE task_mod
  USE bits_mod
  USE io_mod
  USE link_mod
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
    sink%ds    = 0d0
    sink%size  = task_list%head%task%size*0.5
    sink%level = task_list%head%task%level
    sink%position = (i-(n_sink+1)/2.0)*0.25 + task_list%position
    sink%velocity = [0.0,0.0,0.0]
    call sink%set (bits%internal)
    print '(a,2i6,3(2x,3f6.3))', &
      'setup:', i, sink%id, sink%size, sink%position, task_list%position
    allocate (sink%link)
    call sink%link%init
    sink%link%task => sink
    call task_list%append (sink%link)
  end do
  call task_list%init_all_nbors
END SUBROUTINE setup
