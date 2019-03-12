!===============================================================================
!===============================================================================
MODULE pebbles_mod
  USE task_list_mod
  USE pebble_mod
  USE trace_mod
  implicit none
  private
  type, public:: pebbles_t
    type(task_list_t), pointer:: task_list
  contains
    procedure:: init
  end type
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self)
  class(pebbles_t):: self
  type(pebble_t), dimension(:), pointer:: pebble
  class(task_t), pointer:: task
  integer:: i, npebble=100
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  call trace_begin ('pebbles_t%init')
  allocate (self%task_list)
  call self%task_list%list_t%init ('pebbles')
  allocate (pebble(npebble))
  do i=1,npebble
    call pebble(i)%init
    task => pebble(i)
    call self%task_list%append (task)
  end do
  call trace_end
END SUBROUTINE init

END MODULE pebbles_mod
