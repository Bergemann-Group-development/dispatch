PROGRAM experiment
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a single star with oribiting particles
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE star_mod
  USE particles_mod
  USE omp_timer_mod
  implicit none
  type (particles_t):: particles
  !.............................................................................
  call particles%init
  call tic
  do while (particles%time < particles%end_time)
    call particles%update
  end do
  call toc ('particle update', particles%n_step)
!> =============================================================================
END PROGRAM experiment
