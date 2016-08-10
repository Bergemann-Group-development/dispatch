!> =============================================================================
!> Test program for MPI parallel I/O.  gsize3 is the global size of the array,
!> lsize3 is the rank-thread-local size of the array, and start3 is the offset
!> into the global array.  Each thread handles a separate array f, with an 
!> offset (in this test) that increases with MPI rank in the y-direction and 
!> with OMP thread in the z-direction
!> =============================================================================
PROGRAM test_mpi_io
  USE io_mod
  USE io_unit_mod
  USE mpi_mod
  USE mpi_io_mod
  USE omp_mod
  USE omp_timer_mod
  implicit none
  include 'mpif.h'
  type(mpi_io_t), save:: io_f, io_h
  integer, dimension(3):: gsize3, lsize3, start3
  integer, dimension(1):: gsize1, lsize1, start1
  integer:: rec, nrec=10
  integer, parameter:: n_patch=16
  integer, parameter:: n=32, nh=20
  real, dimension(n,n,n), save:: f
  real, dimension(nh), save:: h
  !$omp threadprivate (f,h)
  !-----------------------------------------------------------------------------
  call mpi%init
  call io%init ('test')
  io%verbose = 0
  !$omp parallel default(shared) private(gsize3,lsize3,start3,gsize1,lsize1,start1)
  if (mpi%mode == MPI_THREAD_MULTIPLE+9999) then
    if (io_unit%master) print*,'test multiple_thread mode'
    call test_multiple_thread
  else
    if (io_unit%master) print*,'test single_thread mode'
    call test_single_thread
  end if
  !$omp end parallel
  call mpi%end
contains

!> =============================================================================
!> Multiple thread parallel I/O, seems to not be a working option (even with
!> thread-private io_f and io_h).
!> =============================================================================
subroutine test_multiple_thread
  type(mpi_io_t), save:: io_f, io_h
  !$omp threadprivate (io_f,io_h)
  gsize3 = [n, mpi%size*n, omp%nthreads*n]
  lsize3 = shape(f)
  start3 = [0, mpi%rank*n, omp%thread*n]
  write (io_unit%log,*) 'gsize3:',gsize3
  write (io_unit%log,*) 'lsize3:',lsize3
  write (io_unit%log,*) 'start3:',start3
  !$omp critical (open_cr)
  call io_f%init (gsize3, lsize3, start3)
  call io_f%openw ('data/f.dat')
  !$omp end critical (open_cr)
  !$omp barrier
  !
  lsize1 = shape(h)
  gsize1 = lsize1*mpi%size*omp%nthreads
  start1 = lsize1*(omp%thread+omp%nthreads*mpi%rank)
  write (io_unit%log,*) 'gsize1:',gsize1
  write (io_unit%log,*) 'lsize1:',lsize1
  write (io_unit%log,*) 'start1:',start1
  call io_h%init (gsize1, lsize1, start1)
  call io_h%openw ('data/h.dat')
  !$omp barrier
  do rec=1,nrec
    f=rec
    h=rec
    call tic
    call io_f%write (f, rec)
    call io_h%write (h, rec)
    call toc('write',n**3)
  end do
  do rec=1,nrec
    call tic
    if (.not.io_f%test()) call io_f%wait
    if (.not.io_h%test()) call io_h%wait
    call toc('wait',n**3)
  end do
  !$omp barrier
  call io_f%close
  call io_h%close
  !
  !$omp critical (open_cr)
  call io_f%openr ('data/f.dat')
  call io_h%openr ('data/h.dat')
  !$omp end critical (open_cr)
  !$omp barrier
  do rec=1,nrec
    call io_f%read (f, rec)
    call io_h%read (h, rec)
    if (.not.io_f%test()) call io_f%wait
    if (.not.io_h%test()) call io_h%wait
  end do
  !$omp barrier
  !$omp critical (close_cr)
  call io_f%close
  call io_h%close
  !$omp end critical (close_cr)
end subroutine

!> =============================================================================
!> Single thread parallel I/O, seems to be the only working option
!> =============================================================================
subroutine test_single_thread
  !-----------------------------------------------------------------------------
  ! Set thread local windows to the file
  !-----------------------------------------------------------------------------
  gsize3 = [n, mpi%size*n, omp%nthreads*n]
  lsize3 = shape(f)
  start3 = [0, mpi%rank*n, omp%thread*n]
  lsize1 = shape(h)
  gsize1 = lsize1*mpi%size*omp%nthreads
  start1 = lsize1*(omp%thread+omp%nthreads*mpi%rank)
  !-----------------------------------------------------------------------------
  ! Open the two files -- note that omp single implies barrier at end
  !-----------------------------------------------------------------------------
  !$omp single
  write (io_unit%log,*) 'gsize3:',gsize3
  write (io_unit%log,*) 'lsize3:',lsize3
  write (io_unit%log,*) 'start3:',start3
  call io_f%init (gsize3, lsize3, start3)
  call io_f%openw ('data/f.dat')
  write (io_unit%log,*) 'gsize1:',gsize1
  write (io_unit%log,*) 'lsize1:',lsize1
  write (io_unit%log,*) 'start1:',start1
  call io_h%init (gsize1, lsize1, start1)
  call io_h%openw ('data/h.dat')
  !$omp end single
  !-----------------------------------------------------------------------------
  ! Do non-blocking write
  !-----------------------------------------------------------------------------
  !$omp critical (write_cr)
  do rec=1,nrec
    f=rec
    h=rec
    call tic
    call io_f%write (f, rec)
    call io_h%write (h, rec)
    call toc('write',n**3)
  end do
  !-----------------------------------------------------------------------------
  ! Test and wait
  !-----------------------------------------------------------------------------
  do rec=1,nrec
    call tic
    if (.not.io_f%test()) call io_f%wait
    if (.not.io_h%test()) call io_h%wait
    call toc('wait',n**3)
  end do
  !$omp end critical (write_cr)
  !$omp barrier
  !$omp single
  !-----------------------------------------------------------------------------
  ! Close and re-open -- with implied MPI barriers.  Note omp barrier before
  !-----------------------------------------------------------------------------
  call io_f%close
  call io_h%close
  !!call mpi%barrier
  call io_f%openr ('data/f.dat')
  call io_h%openr ('data/h.dat')
  !$omp end single
  !$omp critical (read_cr)
  !-----------------------------------------------------------------------------
  ! Do non-blocking read (should really be to separate buffers), then wait
  !-----------------------------------------------------------------------------
  do rec=1,nrec
    call tic
    call io_f%read (f, rec)
    call io_h%read (h, rec)
    call toc('read',n**3)
  end do
  do rec=1,nrec
    call tic
    if (.not.io_f%test()) call io_f%wait
    if (.not.io_h%test()) call io_h%wait
    call toc('wait',n**3)
  end do
  !$omp end critical (read_cr)
  !$omp barrier
  !-----------------------------------------------------------------------------
  ! Wait for all threads, then close
  !-----------------------------------------------------------------------------
  !$omp single
  call io_f%close
  call io_h%close
  !$omp end single
end subroutine

END PROGRAM
