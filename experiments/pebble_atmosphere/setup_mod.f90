MODULE setup_mod
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set a planet at a given position, and a set of Rubik's Cube patches
!> surrounding it, representing its extended atmosphere.
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE bits_mod
  USE omp_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  USE planet_mod
  USE atmosphere_mod
  USE eos_mod
  USE experiment_mod
  USE pebble_trace_mod
  !USE radiative_transfer_mod
  USE rubiks_mod
  USE shared_mod
  USE list_mod
  USE task_mod
  implicit none
  private
  type, public, extends(rubiks_t):: setup_t
  contains
    procedure:: init
    procedure, nopass:: close
  end type
  type(setup_t), public:: setup
CONTAINS

!===============================================================================
!> Setup a general Rubik hierarchy: NxNzN patches per level, and zoom factor M
!>
!> Rubik's cube has N=3, M=3
!> Rubik's revenge has N=4, M=2
!>
!> We set the position of the planet as the origin of the coordinate system,
!> which is natural when treating the gravitational influence of the Sun by going
!> into a comoving and rotating coordinate system.  The basic size unit is the
!> size of the central patch, which is taken to be 10% larger than the planet.
!===============================================================================
SUBROUTINE init (self, label)
  class(setup_t):: self
  character(len=*), optional:: label
  !.............................................................................
  class(link_t), pointer:: link
  class(task_t), pointer :: task
  class(experiment_t), pointer:: patch
  !-----------------------------------------------------------------------------
  call trace_begin ('setup_t%init')
  call planet%init
  call atmosphere%init
  self%size = planet%radius*2.4
  planet%position = 0.0
  self%position = planet%position
  call self%rubiks_t%init_hierarchy (label)
  !if (eos%do_rt) then
    !call rt_setup (self%task_list,2)
  !end if
  print*, 'task_list%n =', self%task_list%n
  call trace_end
END SUBROUTINE init

SUBROUTINE close
  call pebble_trace%close_all
END SUBROUTINE close

END MODULE setup_mod
