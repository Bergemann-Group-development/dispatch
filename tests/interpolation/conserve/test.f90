!===============================================================================
!> Test code for conserve_mod
!===============================================================================
PROGRAM test
  USE mpi_mod
  USE io_mod
  USE conserve_mod
  !-----------------------------------------------------------------------------
  call mpi%init
  call io%init
  call conserve%test
END PROGRAM
