!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Set up a cylindrical annulus at a particular radius in a PPD and embed a 
!> small, spherical planetesimal in it. The spherical planetesimal is surrounded
!> by (a) spherical mesh(es). It is assumed that the the radius of the planetesimal
!> R_p < Bondi radius < Hill radius, and thus the planetesimal's gravity is not
!> important.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
PROGRAM dispatch
  USE setup_mod
  USE mpi_mod
  USE io_mod
  implicit none
  !.............................................................................
  call mpi%init
  call io%init ('planetesimal_annulus')
  call setup%init
  call setup%task_list%execute
  call mpi%end
  print*,'All tasks finished. Goodbye.'
!> =============================================================================
END PROGRAM dispatch
