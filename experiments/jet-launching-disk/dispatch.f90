!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Setup a jet launching disk in cylindrical or spherical coords., including a
!> gravitating point mass for the central star.
!> Current setup follows Suriano et al. (2017, 2018).
!> Particles are added on top via "pebble patches".
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
PROGRAM dispatch
  USE setup_mod
  USE dispatcher_mod
  USE experiment_setup_mod
  implicit none
  !.............................................................................
  call setup%init
  call dispatcher%init
  call experiment_setup%init
  call dispatcher%execute(experiment_setup%task_list)
  call setup%end
  print*,'All tasks finished. Goodbye.'
!===============================================================================
END PROGRAM dispatch
