!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Pebble data type
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE pebble_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! A single pebble in a linked list
  !-----------------------------------------------------------------------------
  type, public:: pebble_t
    integer:: id                        ! pebble nr
    integer:: status=0                  ! status bits
    real::  p(3)                        ! position
    real::  v(3)                        ! velocity
    real:: vg(3)                        ! gas velocity
    real:: ag(3)                        ! acceleration of gravity
    real:: ad(3)                        ! drag acceleration
    real:: m                            ! mass
    real:: s                            ! size
    real:: d                            ! solid density
    real:: dg                           ! gas density
    real:: c                            ! drag acceleration rate
    real:: dt                           ! time step
    real(8):: t                         ! particle time
    type(pebble_t), pointer:: prev=>null(), next=>null()
  contains
    procedure:: init
  end type
  integer, save:: id = 0
CONTAINS

!> =============================================================================
!> Give each pebble a uniq id
!> =============================================================================
SUBROUTINE init (self)
  class(pebble_t):: self
  !$omp atomic
  id = id+1
  self%id = id
END SUBROUTINE init

END MODULE
