PROGRAM experiment
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a single star with orbiting particles
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE particles_mod
  USE omp_timer_mod
  USE standard_disk_mod
  implicit none
  type (particle_task_t):: particles
  !.............................................................................
  call io%init
  call particles%init
  call particles%stat
  call tic
  !-----------------------------------------------------------------------------
  ! Advance time in steps "dtime", to mimic conditions when patch tasks are also
  ! being updated. "dtime" represents typical patch update times.
  !-----------------------------------------------------------------------------
  do while (particles%time < particles%end_time)
    particles%time = particles%time + particles%dtime
    call particles%update_kdk
    call particles%output('kdk')
  end do
  call toc (' particle update', particles%n_step)
  !-----------------------------------------------------------------------------
  ! Calculate and print statistics
  !-----------------------------------------------------------------------------
  particles%dt_aver = particles%dt_aver/particles%n_step
  particles%dt_harmonic = particles%n_step/particles%dt_harmonic
  particles%e_kin = particles%e_kin/particles%n_step
  particles%e_aver = particles%e_aver/particles%n_step
  particles%e2_aver = particles%e2_aver/particles%n_step
  call particles%stat
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
  call particles%stat
  particles%time = 0.0
  call tic
  !-----------------------------------------------------------------------------
  ! Advance time in steps "dtime", to mimic conditions when patch tasks are also
  ! being updated. "dtime" represents typical patch update times.
  !-----------------------------------------------------------------------------
  do while (particles%time < particles%end_time)
    particles%time = particles%time + particles%dtime
    call particles%update_rk3
    call particles%output('rk3')
  end do
  call toc (' particle update', particles%n_step)
  !-----------------------------------------------------------------------------
  ! Calculate and print statistics
  !-----------------------------------------------------------------------------
  particles%dt_aver = particles%dt_aver/particles%n_step
  particles%dt_harmonic = particles%n_step/particles%dt_harmonic
  particles%e_kin = particles%e_kin/particles%n_step
  particles%e_aver = particles%e_aver/particles%n_step
  particles%e2_aver = particles%e2_aver/particles%n_step
  call particles%stat
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

!> =============================================================================
END PROGRAM experiment
