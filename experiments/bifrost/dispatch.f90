!> =============================================================================
!> Generic DISPATCH main program for Cartesian meshes.  Note that this level
!> of the code is NOT code specific; it would / could be the same for any 
!> DISPATCH experiment.  Calls to the external code are ONLY made in experiment_t.
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
!> =============================================================================
PROGRAM dispatch_prog
  USE dispatch_mod
  USE extern_mod
  USE cartesian_mod
  type (cartesian_t):: cartesian
  !.............................................................................
  call dispatch%init
  call extern%init
  call cartesian%init
  call cartesian%task_list%execute
  call extern%end
  call dispatch%end
!> =============================================================================
END PROGRAM dispatch_prog
