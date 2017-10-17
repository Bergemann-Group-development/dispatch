PROGRAM dispatch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set of Cartesian patches on a task list, and update until finished
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE eos_mod
  USE setup_mod
  USE cartesian_mod
  USE radiative_transfer_mod
  implicit none
  type (cartesian_t):: cartesian
  real:: position(3)=0.0, size(3)=1.0
  !.............................................................................
  call setup%init
  call cartesian%init                   ! Initialize the patch list and patches
  if (eos%do_rt) &
    call rt_setup (cartesian%task_list)
  call cartesian%task_list%execute      ! update the list = update member patches
  call setup%end
!> =============================================================================
END PROGRAM dispatch
