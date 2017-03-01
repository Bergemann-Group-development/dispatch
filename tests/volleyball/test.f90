!> =============================================================================
!> Program to test and demonstrate the volleyball patch set generator.
!> =============================================================================
PROGRAM test
  USE io_mod
  USE task_mod
  USE patch_mod
  USE volleyball_mod
  USE link_mod
  implicit none

  type(volleyball_t):: volleyball
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call io%init
  call volleyball%init

  link => volleyball%task_list%head
  do while (associated(link))
    task => link%task
    select type (task)
    class is (patch_t)
      call task%output_legacy()
    end select
    link => link%next
  end do

END PROGRAM
