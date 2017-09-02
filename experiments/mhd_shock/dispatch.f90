PROGRAM dispatch
  USE setup_mod
  USE cartesian_mod
  type (cartesian_t):: cartesian
  !.............................................................................
  call setup%init
  call cartesian%init ('mhd_shock')
  call cartesian%task_list%execute     ! update the list = update member patches
  call setup%end
!> =============================================================================
END PROGRAM dispatch
