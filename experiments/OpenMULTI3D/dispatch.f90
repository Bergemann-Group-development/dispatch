!===============================================================================
!> OpenMULTI3D main program.
!>
!> The call hierarchy is
!>
!> dispatch
!>   setup%init                 ! MPI, I/O and log files
!>   dispatcher%init            ! select and initialize the dispatcher
!>   read_external%init         ! read the atmosphere and build the task list
!>   dispatcher%execute         ! run the task list
!>   write_linelist_out         ! line-list diagnostics, if collected
!>   setup%end                  ! shut down
!>
!===============================================================================
PROGRAM dispatch
  USE mpi_mod
  USE os_mod
  USE setup_mod
  USE console_mod
  USE dispatcher_mod
  USE task_list_mod
  USE experiment_mod
  USE read_external_mod
  USE timer_mod
  USE mpi_debug_mod
  USE m3d_patch_mod, only: write_linelist_out
  class(task_list_t), pointer:: task_list => null()
  !.............................................................................
  call setup%init                               ! Generic setup (MPI, I/O, ...)
  call dispatcher%init                          ! Initialize the dispatcher
  call check_mpi_debug()
  call read_external%init (task_list)           ! Initialize task_list
  call dispatcher%execute (task_list)           ! Run dispatcher on the task_list
  call write_linelist_out()                     ! Line-list diagnostics, if collected
  if (timer%verbose > 0) call timer%print() 
  call mpi%barrier()
  if (console_on(ch_m3d, 0, all_ranks=.true.)) &
    write(console_unit,'(a, I3, a, F7.2, a)') 'maximum memory on rank', mpi%rank, ' was', os%max_mem, 'Gb'
  if (console_on(ch_m3d, 0)) write(console_unit,'(a)') 'OpenMULTI3D terminated normally.'
  call setup%end                                ! End setup
!===============================================================================
END PROGRAM dispatch
  