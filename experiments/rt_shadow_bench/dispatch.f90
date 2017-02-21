PROGRAM dispatch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set of Cartesian patches on a task list, and update until finished
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE mpi_mod
  USE cartesian_mod
  USE radiative_transfer_mod
  USE eos_mod
  type (cartesian_t):: cartesian
  type (task_list_t), pointer:: task_list
  !.............................................................................
  call mpi%init
  call io%init('shadow benchmark')
  call eos%init()
  call cartesian%init                   ! Initialize the patch list and patches
  task_list => cartesian%task_list
  call rt_setup (task_list)             ! add RT tasks to the task list
  call task_list%execute                ! update the list = update member patches
  call mpi%end
!> =============================================================================
END PROGRAM dispatch
