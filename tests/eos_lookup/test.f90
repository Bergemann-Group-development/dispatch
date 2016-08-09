!> =============================================================================
!> =============================================================================
PROGRAM test_eos
  USE io_mod
  USE mpi_mod
  USE omp_timer_mod
  USE eos_mod
  implicit none
  real, pointer, dimension(:,:,:):: rho, ee, pg
  integer:: n=4
  !-----------------------------------------------------------------------------
  call mpi%init
  call io%init ('eos')
  allocate (rho(n,n,n), ee(n,n,n), pg(n,n,n))
  rho=1.0
  ee=4.0
  call eos%lookup (rho, ee, pg=pg)
  print *,'pg =',pg(1:n,1,1)
END PROGRAM
