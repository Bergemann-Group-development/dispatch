!> =============================================================================
!> =============================================================================
MODULE planets_mod
  USE atmosphere_mod
  USE planet_mod
  USE particle_mod
  USE task_mod
  USE task_list_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  type, public:: planets_t
    type(task_list_t), pointer:: task_list
  contains
    procedure:: init
  end type
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self)
  class(planets_t):: self
  type(planet_t), dimension(:), pointer:: planet
  type(atmosphere_t), dimension(:), pointer:: atmosphere
  class(task_t), pointer:: task
  integer:: i, nplanet=2
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  call trace_begin ('planets_t%init')
  allocate (self%task_list)
  call self%task_list%list_t%init ('planets')
  allocate (planet(nplanet), atmosphere(nplanet))
  do i=1,nplanet
    call planet(i)%init
    task => planet(i)
    call self%task_list%append (task)
    call atmosphere(i)%init (planet(i))
    call self%task_list%append_list (atmosphere(i)%task_list)
  end do
  call trace_end
END SUBROUTINE init

END MODULE planets_mod
