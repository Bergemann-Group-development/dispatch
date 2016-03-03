MODULE setup_mod
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set a planet at a given position, and a set of Rubik's Cube patches
!> surrounding it, representing its extended atmosphere.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  USE planet_mod
  USE atmosphere_mod
  USE rubiks_cube_mod
  implicit none
  private
  type, public, extends(rubiks_cube_t):: setup_t
  contains
    procedure init
  end type
  type(setup_t), public:: setup
!> =============================================================================
CONTAINS

SUBROUTINE init (self, label)
  class(setup_t):: self
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  call trace_begin ('setup_t%init')
  call scaling%init
  call planet%init
  call atmosphere%init
  self%size = planet%radius*1.1
  print*, 'size1', self%size*scaling%l/cgs%r_earth
  self%position = planet%position
  call self%rubiks_cube_t%init (label)
  call trace_end
END SUBROUTINE init

END MODULE setup_mod
