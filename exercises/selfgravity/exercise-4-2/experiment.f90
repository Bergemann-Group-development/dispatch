PROGRAM experiment
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a single patch with a density tophat and call the Poisson solver
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE pde_mod
  !USE patch_list_mod
  !type (patch_list_t):: list
  type(pde_t):: BE
  !.............................................................................
  call BE%init
  call BE%solve
  call BE%output
  !call list%init ('patch list')         ! Initialize the patch list and patches
  !do while (list%n > 0)                 ! As long as there are patches in the list
  !  call list%update                    ! update the list = update member patches
  !end do
  !print*,'stopping: all tasks finished'
!> =============================================================================
END PROGRAM experiment
