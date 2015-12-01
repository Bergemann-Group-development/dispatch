PROGRAM experiment
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Make list of links with pointers to patches, and with procedures to initialize
!> and update the patches.  Keep updating until the list is empty
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE patch_list_mod
  type (patch_list_t):: list
  !.............................................................................
  call io%init                          ! Initialize input (input.nml) / output
  call list%init ('patch list')         ! Initialize the patch list and patches
!  call list%split (list%head)           ! try splitting a single patch
  do while (list%n > 0)                 ! As long as there are patches in the list
    call list%update                    ! update the list = update member patches
  end do
  call list%print ('after split of first patch')
!> =============================================================================
END PROGRAM experiment
