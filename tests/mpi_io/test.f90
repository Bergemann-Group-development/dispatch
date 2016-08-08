!> =============================================================================
!> =============================================================================
PROGRAM test_mpi_io
  USE io_mod
  USE mpi_mod
  USE mpi_io_mod
  USE omp_timer_mod
  implicit none
  type(mpi_io_t):: mpi_io
  integer, dimension(3):: gsize, lsize, start
  integer:: rec, nrec=10
  integer, parameter:: n=16
  real, dimension(n,n,n):: f
  !-----------------------------------------------------------------------------

  call mpi%init
  io%verbose = 1
  gsize = shape(f)
  lsize = 2
  lsize = gsize
  start = 0
  if (io%master) then
    print*,'gsize:',gsize
    print*,'lsize:',lsize
    print*,'start:',start
  end if
  call mpi_io%init (gsize, lsize, start)
  call mpi_io%openw ('data/test.dat')
  do rec=1,nrec
    f=rec
    call mpi_io%write (f, rec)
    if (.not.mpi_io%test()) call mpi_io%wait
  end do
  call mpi_io%close

  call mpi_io%openr ('data/test.dat')
  do rec=1,nrec
    call mpi_io%read (f, rec)
    if (.not.mpi_io%test()) call mpi_io%wait
  end do
  call mpi_io%close
END PROGRAM
