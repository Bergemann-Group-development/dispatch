PROGRAM dispatch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a single patch with a density tophat and call the Poisson solver
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE mpi_mod
  USE cartesian_mod
  type (cartesian_t):: cartesian
  !.............................................................................
  call io%init
  !call mpi%init
  call cartesian%init('orszag-tang')    ! Initialize the patch list and patches
  call cartesian%task_list%execute      ! update the list = update member patches
  print*,'stopping: all tasks finished'
!> =============================================================================
END PROGRAM dispatch
