PROGRAM experiment
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Make list of links with pointers to patches, and with procedures to initialize
!> and update the patches.  Keep updating until the list is empty
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE lists_mod
  type (list_t):: list
  !.............................................................................
  call io%init                          ! Initialize input (input.nml) / output
  call list%init                        ! Initialize the patches
  do while (list%n > 0)                 ! As long as there are patches in the list
    call list%update                    ! update the list = update member patches
  end do
!> =============================================================================
END PROGRAM experiment
