!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Atmosphere component, with a number of tasks centered around a planet
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE atmosphere_mod
  USE planet_mod
  USE task_mod
  USE patch_mod
  USE task_list_mod
  implicit none
  private
  type, public:: atmosphere_t
    class(task_list_t), pointer:: task_list
    integer:: npatch=100
  contains
    procedure:: init
  end type
CONTAINS

!> =============================================================================
!> Initialize a number of atmosphere patches around a planet
!> =============================================================================
SUBROUTINE init (self, planet)
  class(atmosphere_t):: self
  type(planet_t):: planet
  class(patch_t), pointer:: patch
  class(task_t), pointer:: task
  integer:: i
  !.............................................................................
  allocate (self%task_list)
  do i=1,self%npatch
    allocate (patch)
    call patch%init
    patch%position = patch%position + planet%position
    task => patch
    call self%task_list%append (task)
  end do
END SUBROUTINE init
END MODULE atmosphere_mod
