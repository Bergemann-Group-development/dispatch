!===============================================================================
!> Program to test and demonstrate the volleyball patch set generator.
!===============================================================================
PROGRAM test
  USE io_mod
  USE volleyball_layers_mod
  implicit none
  type(volleyball_layers_t):: volleyball_layers
  !-----------------------------------------------------------------------------
  call io%init
  call volleyball_layers%init
END PROGRAM
