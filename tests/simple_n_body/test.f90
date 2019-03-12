!===============================================================================
!> Program to test functionality and speed of particles, when using an array
!> of structures
!===============================================================================
PROGRAM test
  USE particle_mod
  USE particles_mod
  USE particle_history_mod
  USE omp_timer_mod
  USE timer_mod
  implicit none
  type (particles_t):: particles
  integer:: n
  real(8):: wc

  call timer%init
  n = 3000
  read *,n
  call particles%init (n)
  print *, particles%n
  call particles%randomize
  !call particles%print
  wc = wallclock()
  call particles%move_past (1.0)
  wc = wallclock()-wc
  call particles%print
  print *,wc, wc/particles%moves
  call timer%print
END PROGRAM

