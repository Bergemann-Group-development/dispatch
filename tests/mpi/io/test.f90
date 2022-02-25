!===============================================================================
!> Test program for MPI parallel I/O.  gsize3 is the global size of the array,
!> lsize3 is the rank-thread-local size of the array. Each thread handles a
!> separate array f, with an offset (in this test) that increases with MPI rank
!> in the y-direction and with OMP thread in the z-direction
!===============================================================================
PROGRAM test_mpi_io
  USE io_mod
  USE io_unit_mod
  USE mpi_mod
  USE mpi_io_mod
  USE omp_mod
  USE omp_timer_mod
  USE timer_mod
  implicit none
  include 'mpif.h'
  type(mpi_file_t):: file1, file2
  type(mpi_io_t), save:: io_f, io_h, io_1
  integer:: nrec=10
  integer, parameter:: n_patch=16
  integer, parameter:: n=32, nh=20
  real, dimension(n,n,n), save:: f
  real, dimension(nh), save:: h
  integer, save:: verbose=1
  !$omp threadprivate (f,h,io_f,io_h)
  !-----------------------------------------------------------------------------
  call mpi%init
  call io%init ('test')
  io%verbose = 0
  !$omp parallel default(shared)
  if (io_unit%master) print*,'test single_thread mode'
  call test_single_thread
  !$omp end parallel
  call mpi%end
  call check_files
contains

!===============================================================================
!> This subroutine is called by each thread, so variable declared here become
!> private
!===============================================================================
subroutine test_single_thread
  integer:: rec
  integer, dimension(1):: gsize1, lsize1
  integer, dimension(3):: gsize3, lsize3
  !-----------------------------------------------------------------------------
  ! Set thread local windows to the file
  !-----------------------------------------------------------------------------
  gsize3 = [n, mpi%size*n, omp%nthreads*n]
  lsize3 = [n, n, n]
  lsize1 = shape(h)
  gsize1 = lsize1*mpi%size*omp%nthreads
  !-----------------------------------------------------------------------------
  ! Open the two files -- note that omp single implies barrier at end
  !-----------------------------------------------------------------------------
  call file1%openw ('data/f.dat')
  call file2%openw ('data/h.dat')
  call io_f%use (file1)
  call io_h%use (file2)
  !-----------------------------------------------------------------------------
  ! Do blocking write
  !-----------------------------------------------------------------------------
  call io_h%init (product(gsize1), product(lsize1))
  call io_f%init (product(gsize3), product(lsize3))
  call io_f%check_init
  do rec=1,nrec
    f=rec*100+mpi%rank*10+omp%thread
    h=rec*100+mpi%rank*10+omp%thread
    call tic
    !if (mpi%rank /= 1) then
      call io_f%write (f, rec, omp%thread+1)
      if (verbose>1) &
        print '(a,3i5,2f12.0)', ' write: rank, thread, rec, min, max', &
        mpi%rank, omp%thread, rec, minval(f), maxval(f)
    !end if
    call io_h%write (h, rec, omp%thread+1)
    call timer%toc('write',n**3)
  end do
  !-----------------------------------------------------------------------------
  ! Do non-blocking write
  !-----------------------------------------------------------------------------
  !$omp barrier
  do rec=1,nrec
    f=rec*100+mpi%rank*10+omp%thread
    h=rec*100+mpi%rank*10+omp%thread
    call tic
    if (mpi%rank /= 1) then
      call io_f%iwrite (f, rec, omp%thread+1)
      if (verbose>1) &
        print '(a,3i5,2f12.0)', ' write: rank, thread, rec, min, max', &
        mpi%rank, omp%thread, rec, minval(f), maxval(f)
    end if
    call io_h%iwrite (h, rec, omp%thread+1)
    call timer%toc('iwrite',n**3)
  end do
  !-----------------------------------------------------------------------------
  ! Test and wait
  !-----------------------------------------------------------------------------
  do rec=1,nrec
    call tic
    if (.not.io_f%test()) call io_f%wait
    if (.not.io_h%test()) call io_h%wait
    call timer%toc('wait',n**3)
  end do
  !-----------------------------------------------------------------------------
  ! Close and re-open -- with implied MPI barriers.  Note omp barrier before
  !-----------------------------------------------------------------------------
  !$omp barrier
  call file1%close
  call file2%close
  call file1%openr ('data/f.dat')
  call file2%openr ('data/h.dat')
  call io_f%use (file1)
  call io_h%use (file2)
  !-----------------------------------------------------------------------------
  ! Do blocking read (should really be to separate buffers)
  !-----------------------------------------------------------------------------
  call io_h%init (product(gsize1), product(lsize1))
  call io_f%init (product(gsize3), product(lsize3))
  call io_f%check_init
  do rec=1,nrec
    call tic
    if (mpi%rank /= 1) then
      call io_f%read (f, rec, omp%thread+1)
      if (verbose>1) &
        print '(a,3i5,2f12.0)', ' read: rank, thread, rec, min, max', &
        mpi%rank, omp%thread, rec, minval(f), maxval(f)
      call io_h%read (h, rec, omp%thread+1)
    end if
    call timer%toc('read',n**3)
  end do
  !-----------------------------------------------------------------------------
  ! Wait for all threads, then close
  !-----------------------------------------------------------------------------
  !$omp barrier
  call file1%close
  call file2%close
end subroutine

!===============================================================================
!===============================================================================
subroutine check_files
  real, dimension(:,:,:), allocatable:: tmp
  integer:: rec, recl
  !-----------------------------------------------------------------------------
  recl = 4*n**3
  open (1, file='data/f.dat', access='direct', status='old', recl=recl)
  allocate (tmp(n,n,n))
  do rec=1,nrec
    read (1, rec=rec) tmp
    print *, rec, minval(tmp), maxval(tmp)
  end do
  close (1)
end subroutine check_files

END PROGRAM
