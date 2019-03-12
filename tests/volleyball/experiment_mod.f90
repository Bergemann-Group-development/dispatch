!===============================================================================
!> Dummy experiment module for the volleyball test.
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE mhd_mod
  implicit none
  private

  type, public, extends(mhd_t):: experiment_t
  contains
  end type

CONTAINS

END MODULE experiment_mod
