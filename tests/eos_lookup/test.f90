!> =============================================================================
!> Test program for EOS lookup module eos_mod.f90
!> =============================================================================
PROGRAM test_eos
  USE io_mod
  USE mpi_mod
  USE omp_timer_mod
  USE eos_mod
  implicit none
  real, save, pointer, dimension(:,:,:):: rho, ee, pg, tt, lnrk, ne
  real, save, pointer, dimension(:,:,:,:):: lnsource
  integer:: n=100, i
  !$omp threadprivate(pg,tt,lnrk,ne,lnsource)
  !-----------------------------------------------------------------------------
  call mpi%init
  call io%init ('lookup')
  call eos%init (n,n,n)
  allocate (rho(n,n,n), ee(n,n,n))
  rho=1.0
  ee=5.5
  !$omp parallel default(shared)
  allocate (pg(n,n,n), tt(n,n,n), lnrk(n,n,n), ne(n,n,n))
  allocate (lnsource(n,n,n,4))
  call tic
  call eos%lookup (rho, ee, pg=pg)
  call toc ('pg',n**3)
  !$omp barrier
  call tic
  call eos%lookup (rho, ee, tt=tt, ne=ne)
  call toc ('tt,ne',n**3)
  !$omp barrier
  call tic
  call eos%lookup (rho, ee, lnrk=lnrk, lnsource=lnsource)
  call toc ('lnrk,lns',n**3)
  !$omp barrier
  call tic
  call eos%lookup (rho, ee, pg=pg, tt=tt, ne=ne, lnrk=lnrk, lnsource=lnsource)
  call toc ('all 8',n**3)
  !$omp end parallel
  print *,'  pg =',  pg(1:4,1,1)
  print *,'  tt =',  tt(1:4,1,1)
  print *,'  ne =',  ne(1:4,1,1)
  print *,'lnrk =',lnrk(1:4,1,1)
  do i=1,4
    print *,' lns =',lnsource(1:4,1,1,i)
  end do
END PROGRAM
