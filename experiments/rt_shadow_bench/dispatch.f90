PROGRAM dispatch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set of Cartesian patches on a task list, and update until finished
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE mpi_mod
  USE cartesian_mod
  USE radiative_transfer_mod
  type (cartesian_t):: cartesian
  type (task_list_t), pointer:: task_list
  real:: position(3)=0.0, size(3)=1.0
  !.............................................................................
  call mpi%init
  call io%init('shadow benchmark')
  call cartesian%init                   ! Initialize the patch list and patches
  call rt_setup (cartesian%task_list)           ! Setup the RT solver
  call cartesian%task_list%execute      ! update the list = update member patches
  call mpi%end
!> =============================================================================
END PROGRAM dispatch
