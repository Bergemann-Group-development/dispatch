PROGRAM dispatch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Run the DISPATCH code, by initializing a patch list and updating it
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE patch_list_mod
  type (patch_list_t):: list
  !.............................................................................
  call list%init ('patch list')         ! Initialize the patch list and patches
  do while (list%n > 0)                 ! As long as there are patches in the list
    call list%update                    ! update the list = update member patches
  end do
  print*,'stopping: all tasks finished'
!> =============================================================================
END PROGRAM dispatch
