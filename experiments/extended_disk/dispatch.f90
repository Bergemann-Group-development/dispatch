!===============================================================================
!> Generic DISPATCH main program for Cartesian meshes
!===============================================================================
PROGRAM dispatch
  USE io_unit_mod
  USE setup_mod
  USE cylinder_mod
  USE dispatcher_mod
  type (cylinder_t):: cylinder
  !.............................................................................
  call setup%init                               ! Standard setup (MPI, I/O, scaling, ...)
  call dispatcher%init                          ! Set dispatcher method
  call cylinder%init                            ! Initialize the patch list and patches
  call cylinder%task_list%initialize            ! Make task list ready to execute
  call dispatcher%execute (cylinder%task_list)  ! Call dispatcher
  write (io_unit%mpi,*) 'calling setup%end'
  flush (io_unit%mpi)
  call setup%end
!===============================================================================
END PROGRAM dispatch
