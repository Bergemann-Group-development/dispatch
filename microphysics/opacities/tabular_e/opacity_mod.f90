!> =============================================================================
!> Equation of state module for any sort of tables, provided by a reader
!> =============================================================================
MODULE opacity_mod
  USE io_mod
  USE trace_mod
  USE opacity_reader_mod
  USE eos_mod
  implicit none
  private
  type, public:: opacity_t
  contains
    procedure, nopass:: init
  end type
  type(opacity_t), public:: opacity
CONTAINS

SUBROUTINE init
  call opacity_reader%init
END SUBROUTINE init

END MODULE opacity_mod
