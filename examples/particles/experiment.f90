PROGRAM experiment
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a single star with oribiting particles
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE star_mod
  USE particles_mod
  USE omp_timer_mod
  implicit none
  type (particle_task_t):: particles
  !.............................................................................
  call particles%init
  call tic
  call particles%update
  !call particles%update_cache
  call toc ('particle update', particles%n_step)

  particles%dt_aver = particles%dt_aver/particles%n_step
  particles%dt_harmonic = particles%n_step/particles%dt_harmonic
  print'(4(a,1p,e10.2/))', &
   '     dt_min:', particles%dt_min, &
   '    dt_aver:', particles%dt_aver, &
   'dt_harmonic:', particles%dt_harmonic, &
   '     dt_max:', particles%dt_max

  call particles%init
  particles%n_step = 0
  call tic
  call particles%update_cache
  call toc ('particle update_cache', particles%n_step)
!> =============================================================================
END PROGRAM experiment
