!===============================================================================
!> Test program for EOS lookup module eos_mod.f90
!===============================================================================
PROGRAM test_eos
  USE io_mod
  USE mpi_mod
  USE timer_mod
  USE itable_mod
  implicit none
  real, pointer, contiguous:: lnr(:,:,:), ee(:,:,:), v(:,:,:)
  type(itable_t):: eos
  integer:: n=300, i, iv, i1, i2, i3, nv
  !-----------------------------------------------------------------------------
  call mpi%init
  call io%init ('lookup')
  call timer%init
  call eos%init ('dispatch.eos')
  nv = eos%nv
  allocate (lnr(n,n,n), ee(n,n,n), v(n,n,n))
  do i3=1,n
  do i2=1,n
  do i1=1,n
    lnr(i1,i2,i3) = 1.0 + 2.0*(mod(i1*10,n)/n-0.5)
    ee (i1,i2,i3) = 6.5 + 2.0*(mod(i2*20,n)/n-0.5)
  end do
  end do
  end do
  do i1=1,5
    call timer%tic
    do iv=1,nv
      call eos%lookup (iv, lnr, ee, v)
    end do
    call timer%toc ('itable',n**3)
  end do
  deallocate (lnr,ee,v)
END PROGRAM
