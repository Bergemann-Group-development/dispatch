PROGRAM test
  USE io_mod
  USE mpi_mod
  USE rt_sc_mod
  implicit none
  type(rt_sc_t):: rt_sc
  !----------------------------------------------------------------------------
  call io%init
  call mpi%init
  call rt_sc%init
  call rt_sc%test
  call rt_sc%update
END PROGRAM test
