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
  integer, save:: n_patch=3, m_zoom=3
  integer:: iostat
  namelist /rubiks_params/ n_patch, m_zoom
  !-----------------------------------------------------------------------------
  call trace_begin ('setup_t%init')
  call io%init
  rewind (io%input)
  read (io%input, rubiks_params, iostat=iostat)
  if (io%master) write (*, rubiks_params)
  call scaling%init
  call planet%init
  call atmosphere%init
  self%size = planet%radius*2.2
  print*, 'Rubik hierarchy size:', self%size*scaling%l/cgs%r_earth
  self%position = planet%position
  call self%rubiks_t%init_hierarchy (n_patch, m_zoom, label)
  call trace_end
END SUBROUTINE init

END MODULE setup_mod
