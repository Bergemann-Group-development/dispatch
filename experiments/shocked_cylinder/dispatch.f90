!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Setup a 2-D shock wave in Cartesian coords. that interacts with a solid
!> cylinder. The cylinder is surrounded by an envelope of cylindrical grids.
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
PROGRAM dispatch
  USE setup_mod
  USE mpi_mod
  USE io_mod
  implicit none
  !.............................................................................
  call mpi%init
  call setup%init
  call setup%task_list%execute
  call mpi%end
  print*,'All tasks finished. Goodbye.'
!===============================================================================
END PROGRAM dispatch
