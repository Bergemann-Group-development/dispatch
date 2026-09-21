!===============================================================================
!> Generic DISPATCH main program.  A selection of components may be added; the
!> first one that defines a task list defines the job.
!>
!> The call hierarchy is
!>
!> dispatch
!>   mpi%init                   ! start up MPI
!>   io%init                    ! assign unit numbers and open log files
!>   experiment%begin           ! generic startup
!>     code_begin               ! code startup
!>   cartesian%init             ! intialize a Cartesian patch arrangement
!>     experiment%init          ! initialize an experiment patch
!>       patch_t%init           ! allocate memory, mesh, etc
!>         task_t%init          ! assign unique task id
!>       code_init              ! initialize a code patch
!>     task_list%append         ! append patch to task list
!>   task_list%execute          ! execute the task list
!>     task_list%update         ! update the task list
!>       experiment%update      ! generic task update
!>         code_update          ! external code time step
!>   experiment%end             ! generic shut down
!>     code_end                 ! external code shut down
!>   mpi%end                    ! close down MPI
!>
!===============================================================================
PROGRAM dispatch
  USE mpi_mod
  USE os_mod
  USE setup_mod
  !USE restart_mod
  USE dispatcher_mod
  !USE cartesian_mod
  USE task_list_mod
  USE experiment_mod
  USE read_external_mod
  USE timer_mod
  USE mpi_debug_mod
  !type(cartesian_t):: cartesian                 ! Use Cartesian patch arrangement
  class(task_list_t), pointer:: task_list => null()
  !integer :: i
  !.............................................................................
  call setup%init                               ! Generic setup (MPI, I/O, ...)
  call dispatcher%init                          ! Initialize the dispatcher
  call check_mpi_debug()
  call read_external%init (task_list)           ! Possibly initialize task_list
  !call timer%print() 
  !call restart%init (task_list)                 ! Possibly initialize task_list
  !call cartesian%init (task_list)               ! .. or else here ...
  !call read_snapshot%init (task_list)           ! Construct ICs
  call dispatcher%execute (task_list)           ! Run dispatcher on the task_list
  if (timer%verbose > 0) call timer%print() 
  call mpi%barrier()
  write(*,'(a, I3, a, F7.2, a)') 'maximum memory on rank', mpi%rank, ' was', os%max_mem, 'Gb'
  write(*,'(a, I3, a, F7.2, a)') 'OpenMULTI3D terminated normally.'
  !call os%print_tracked_mem()
  !i = findloc(os%labels(1:os%nlabels), 'rt_patch%mem', 1)
  !print *, os%n_allocated(i), os%labelled_mem(i)
  !write(*,'(a, I3, a, F7.2, a)') 'still allocated memory on rank', mpi%rank, ' is', os%gb, 'Gb'
  call setup%end                                ! End setup
!===============================================================================
END PROGRAM dispatch
  