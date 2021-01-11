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
  real, pointer, contiguous:: vv(:,:,:,:)
  type(itable_t):: eos
  integer:: n, i, iv, i1, i2, i3, nv
  !-----------------------------------------------------------------------------
  call mpi%init
  call io%init ('lookup')
  call timer%init
  call eos%init ('test.eos')
  do n=16,32,16
    nv = eos%nv
    write (stdout,'(a,i2,a)') 'with n=',n,' and separate lookup of variables'
    allocate (lnr(n,n,n), ee(n,n,n), v(n,n,n))
    do i3=1,n
    do i2=1,n
    do i1=1,n
      lnr(i1,i2,i3) = 1.5 + 2.0*(mod(i1*11,n)/real(n))
      ee (i1,i2,i3) = 6.5 + 2.0*(mod(i2* 7,n)/real(n))
    end do
    end do
    end do
    do i1=1,5
      call timer%tic
      do iv=1,nv
        call eos%lookup (iv, lnr, ee, v)
      end do
      call timer%toc ('itable',nv*n**3)
    end do
    deallocate (v)
    !---------------------------------------------------------------------------
    write (stdout,'(2(a,i2),a)') 'with n=',n,' ,nv=',nv,' and v(n,n,n,nv) ordering'
    allocate (vv(n,n,n,nv))
    do i1=1,5
      call timer%tic
      call eos%lookup (lnr, ee, vv)
      call timer%toc ('itable',nv*n**3)
    end do
    deallocate (vv)
    !---------------------------------------------------------------------------
    write (stdout,'(2(a,i2),a)') 'with n=',n,' ,nv=',nv,' and v(nv,n,n,n) ordering'
    allocate (vv(nv,n,n,n))
    do i1=1,5
      call timer%tic
      call eos%lookup (lnr, ee, vv)
      call timer%toc ('itable',nv*n**3)
    end do
    deallocate (vv)
    !---------------------------------------------------------------------------
    nv = 6
    write (stdout,'(2(a,i2),a)') 'with n=',n,' ,nv=',nv,' and v(n,n,n,nv) ordering'
    allocate (vv(n,n,n,nv))
    do i1=1,5
      call timer%tic
      call eos%lookup (lnr, ee, vv)
      call timer%toc ('itable',nv*n**3)
    end do
    deallocate (vv)
    !---------------------------------------------------------------------------
    write (stdout,'(2(a,i2),a)') 'with n=',n,' ,nv=',nv,' and v(nv,n,n,n) ordering'
    allocate (vv(nv,n,n,n))
    do i1=1,5
      call timer%tic
      call eos%lookup (lnr, ee, vv)
      call timer%toc ('itable',nv*n**3)
    end do
    deallocate (vv)
    deallocate (lnr,ee)
  end do
END PROGRAM
