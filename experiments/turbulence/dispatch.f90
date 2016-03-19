PROGRAM dispatch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set of Cartesian patches on a task list, and update until finished
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE mpi_mod
  USE cartesian_mod
  type (cartesian_t):: cartesian
  type (task_list_t), pointer:: task_list
  real:: position(3)=0.0, size(3)=1.0
  !.............................................................................
  call io%init('turbulence')
 !call mpi%init
  call cartesian%init                   ! Initialize the patch list and patches
  task_list => cartesian%task_list      ! for now, just short hand -- could assemble
  task_list%n_tasks = product(cartesian%dims)
  call task_list%execute                ! update the list = update member patches
  print*,'stopping: all tasks finished'
!> =============================================================================
END PROGRAM dispatch
