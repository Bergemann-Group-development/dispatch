MODULE setup_mod
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a blob and set of Rubik's Cube patches
!> surrounding it. Cast a high temperature radiation from one boundary.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE bits_mod
  USE omp_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  USE experiment_mod
  USE rubiks_mod
  USE radiative_transfer_mod
  implicit none
  private
  type, public, extends(rubiks_t):: setup_t
  contains
    procedure:: init
    procedure, nopass:: close
  end type
  type(setup_t), public:: setup
CONTAINS

!> =============================================================================
!> Setup a general Rubik hierarchy: NxNzN patches per level, and zoom factor M
!>
!> Rubik's cube has N=3, M=3
!> Rubik's revenge has N=4, M=2
!>
!> We set the position of the planet as the origin of the coordinate system,
!> which is natural when treating the gravitational influence of the Sun by going
!> into a comoving and rotating coordinate system.  The basic size unit is the
!> size of the central patch, which is taken to be 10% larger than the planet.
!> =============================================================================
SUBROUTINE init (self, label)
  class(setup_t):: self
  character(len=*), optional:: label
  !.............................................................................
  !-----------------------------------------------------------------------------
  call trace_begin ('setup_t%init')
  call io%init (label)
  call scaling%init
  self%size = 0.14
  self%position = 0.5
  call self%rubiks_t%init_hierarchy (label)
  call self%task_list%remove_parents
  call rt_setup (self%task_list)                                        ! add RT tasks to the task list
  call trace_end
END SUBROUTINE init

SUBROUTINE close
END SUBROUTINE close

END MODULE setup_mod
