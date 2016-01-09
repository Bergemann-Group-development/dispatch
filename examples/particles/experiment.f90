PROGRAM experiment
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a single star with oribiting particles
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  !USE star_mod
  !USE standard_disk_mod
  USE particles_mod
  USE omp_timer_mod
  implicit none
  type (particle_task_t):: particles
  !.............................................................................
  !call star%init
  !call standard_disk%init
  call particles%init
  !
  call tic
  do while (particles%time < particles%end_time)
    particles%time = particles%time + particles%dtime
    call particles%update_list
  end do
  call toc (' particle update_list', particles%n_step)

  particles%dt_aver = particles%dt_aver/particles%n_step
  particles%dt_harmonic = particles%n_step/particles%dt_harmonic
  particles%e_kin = particles%e_kin/particles%n_step
  particles%e_aver = particles%e_aver/particles%n_step
  particles%e2_aver = particles%e2_aver/particles%n_step
  print'((a,1p,e10.2))', &
   '     n_step:', real(particles%n_step), &
   '     dt_min:', particles%dt_min, &
   '    dt_aver:', particles%dt_aver, &
   'dt_harmonic:', particles%dt_harmonic, &
   '     dt_max:', particles%dt_max, &
   '     e_aver:', particles%e_aver, &
   '      e_kin:', particles%e_kin, &
   '      e_rms:', sqrt(particles%e_rms/particles%n_step), &
   '  e_max_rms:', particles%e_max_rms

  call particles%init
  particles%n_step = 0
  call tic
  do while (particles%time < particles%end_time)
    particles%time = particles%time + particles%dtime
    call particles%update_cache
  end do
  call toc ('particle update_cache', particles%n_step)
!> =============================================================================
END PROGRAM experiment
