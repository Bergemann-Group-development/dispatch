MODULE setup_mod
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set a planet at a given position, and a set of Rubik's Cube patches
!> surrounding it, representing its extended atmosphere.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE omp_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  USE planet_mod
  USE atmosphere_mod
  USE rubiks_mod
  implicit none
  private
  type, public, extends(rubiks_t):: setup_t
  contains
    procedure init
  end type
  type(setup_t), public:: setup
CONTAINS

!> =============================================================================
!> Setup a general Rubik hierarchy: NxNzN patches per level, and zoom factor M
!> Rubik's cube has N=3, M=3
!> Rubik's revenge has N=4, M=2
!> =============================================================================
SUBROUTINE init (self, label)
  class(setup_t):: self
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  call trace_begin ('setup_t%init')
  call io%init
  call scaling%init
  call planet%init
  call atmosphere%init
  self%size = planet%radius*2.2
  print*, 'Rubik hierarchy size:', self%size*scaling%l/cgs%r_earth
  planet%position = 0.0
  self%position = planet%position
  call self%rubiks_t%init_hierarchy (label)
  call trace_end
END SUBROUTINE init

END MODULE setup_mod
