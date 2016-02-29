PROGRAM dispatch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set of Cartesian patches on a task list, and update until finished
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE mpi_mod
  USE rubiks_cube_mod
  type (rubiks_cube_t):: planet
  real:: position(3)=0.0, size(3)=1.0
  !.............................................................................
  call io%init
 !call mpi%init
  call planet%init('planet_atmoshpere') ! initialize the patch list and patches
  call planet%task_list%execute         ! execute the list = update member patches
  print*,'stopping: all tasks finished'
!> =============================================================================
END PROGRAM dispatch
