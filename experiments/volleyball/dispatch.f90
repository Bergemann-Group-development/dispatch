!===============================================================================
!> Generic DISPATCH main program for Cartesian meshes
!===============================================================================
PROGRAM dispatch
  USE setup_mod
  USE restart_mod
  USE dispatcher_mod
  USE volleyball_mod
  type (volleyball_t):: volleyball               ! Volleyball patch arrangement
  !.............................................................................
  call setup%init                                ! Standard setup
  call dispatcher%init                           ! Initialize the dispatcher
  call volleyball%init                           ! Initialize the task list
  call dispatcher%execute (volleyball%task_list) ! Run dispatcher on task_list
  call setup%end                                 ! End run
!===============================================================================
END PROGRAM dispatch
