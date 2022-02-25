!===============================================================================
!> Test program for EOS lookup module eos_mod.f90
!===============================================================================
PROGRAM test_eos
  USE io_mod
  USE mpi_mod
  USE timer_mod
  USE eos_mod
  implicit none
  real, save, allocatable, target, dimension(:,:,:):: d, ee, lne, pg, tt, lnd
  integer:: n=200, i, dim(3), ni
  !$omp threadprivate(pg,tt)
  !-----------------------------------------------------------------------------
  call io%init ('lookup')
  call eos%init
  allocate (d(n,n,n), lnd(n,n,n), ee(n,n,n), lne(n,n,n))
  d=1.0
  ee=7.423690
  !$omp parallel default(shared)
  allocate (pg(n,n,n), tt(n,n,n))
  dim = [n,n,n]
  call tic
  call eos%lookup (dim, d=d, ee=ee)
  call timer%toc ('  void',n**3)
  !$omp barrier
  call tic
  call eos%lookup (dim, d=d, ee=ee, pg=pg)
  call timer%toc ('    pg',n**3)
  !$omp barrier
  print *,'using d and ee'
  ni=6
  do i=1,ni
    call tic
    call eos%lookup (dim, d=d, ee=ee, tt=tt, pg=pg, iskip=i)
    call timer%toc (' tt,pg',n**3)
    !$omp barrier
  end do
  lne = log(ee)
  print *,'using d and lne'
  do i=1,ni
    call tic
    call eos%lookup (dim, d=d, lne=lne, tt=tt, pg=pg, iskip=i)
    call timer%toc (' tt,pg',n**3)
    !$omp barrier
  end do
  print *,'using lnd and ee'
  do i=1,ni
    call tic
    call eos%lookup (dim, lnd=lnd, ee=ee, tt=tt, pg=pg, iskip=i)
    call timer%toc (' tt,pg',n**3)
    !$omp barrier
  end do
  lnd = log(d)
  print *,'using lnd and lne'
  do i=1,ni
    call tic
    call eos%lookup (dim, lnd=lnd, lne=lne, tt=tt, pg=pg, iskip=i)
    call timer%toc (' tt,pg',n**3)
    !$omp barrier
  end do
  print *,'using lnd and lne'
  do i=1,ni
    call tic
    call eos%lookup (dim, lnd=lnd, lne=lne, tt=tt, pg=pg, iskip=i)
    call timer%toc (' tt,pg',n**3)
    !$omp barrier
  end do
  !$omp end parallel
  print *,'  pg =', pg(1:4,1,1)
  print *,'  tt =', tt(1:4,1,1)
END PROGRAM
