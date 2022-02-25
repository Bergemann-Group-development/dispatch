!===============================================================================
!> Test program for EOS lookup module eos_mod.f90
!===============================================================================
PROGRAM test_eos
  USE io_mod
  USE mpi_mod
  USE omp_timer_mod
  USE eos_mod
  implicit none
  real, save, pointer, dimension(:,:,:):: rho, ee, pg, tt, rk, ne
  real, save, pointer, dimension(:,:,:,:):: src
  integer:: n=200, i
  !$omp threadprivate(pg,tt,rk,ne,src)
  !-----------------------------------------------------------------------------
  call io%init ('lookup')
  call eos%init (n,n,n)
  allocate (rho(n,n,n), ee(n,n,n))
  rho=1.0
  ee=5.5
  !$omp parallel default(shared)
  allocate (pg(n,n,n), tt(n,n,n), rk(n,n,n), ne(n,n,n))
  allocate (src(n,n,n,4))
  call tic
  call eos%lookup (rho, ee)
  call timer%toc ('  void',n**3)
  !$omp barrier
  call tic
  call eos%lookup (rho, ee, pg=pg)
  call timer%toc ('    pg',n**3)
  !$omp barrier
  call tic
  call eos%lookup (rho, ee, tt=tt, ne=ne)
  call timer%toc (' tt,ne',n**3)
  !$omp barrier
  call tic
  call eos%lookup (rho, ee, rk=rk, src=src)
  call timer%toc ('rk,src',n**3)
  !$omp barrier
  call tic
  call eos%lookup (rho, ee, pg=pg, tt=tt, ne=ne, rk=rk, src=src)
  call timer%toc (' all 8',n**3)
  !$omp end parallel
  print *,'  pg =', pg(1:4,1,1)
  print *,'  tt =', tt(1:4,1,1)
  print *,'  ne =', ne(1:4,1,1)
  print *,'  rk =', rk(1:4,1,1)
  do i=1,4
    print *,' src =',src(1:4,1,1,i)
  end do
END PROGRAM
