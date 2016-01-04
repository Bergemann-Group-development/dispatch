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
  do while (particles%time < particles%end_time)
    call particles%update
  end do
  particles%dt_aver = particles%dt_aver/particles%n_step
  particles%dt_harmonic = particles%n_step/particles%dt_harmonic
  print'(4(a,1p,e10.2/))', &
   '     dt_min:', particles%dt_min, &
   '    dt_aver:', particles%dt_aver, &
   'dt_harmonic:', particles%dt_harmonic, &
   '     dt_max:', particles%dt_max
  call toc ('particle update', particles%n_step)

  call particles%init
  call particles%update_cache

  !call particles%init
  !call particles%update
  !do while (particles%time < particles%end_time)
  !  call particles%update_bags
  !end do
!> =============================================================================
END PROGRAM experiment
