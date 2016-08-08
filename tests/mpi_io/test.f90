!> =============================================================================
!> =============================================================================
PROGRAM test_mpi_io
  USE io_mod
  USE mpi_mod
  USE mpi_io_mod
  USE omp_timer_mod
  implicit none
  type(mpi_io_t):: io_f, io_h
  integer, dimension(3):: gsize3, lsize3, start3
  integer, dimension(1):: gsize1, lsize1, start1
  integer:: rec, nrec=10
  integer, parameter:: n=16, nh=20
  real, dimension(n,n,n):: f
  real, dimension(nh):: h
  !-----------------------------------------------------------------------------

  call mpi%init
  io%verbose = 1
  gsize3 = shape(f)
  lsize3 = gsize3
  start3 = 0
  if (io%master) then
    print*,'gsize3:',gsize3
    print*,'lsize3:',lsize3
    print*,'start3:',start3
  end if
  call io_f%init (gsize3, lsize3, start3)
  call io_f%openw ('data/f.dat')
  !
  gsize1 = shape(h)
  lsize1 = gsize1
  start1 = 0
  if (io%master) then
    print*,'gsize1:',gsize1
    print*,'lsize1:',lsize1
    print*,'start1:',start1
  end if
  call io_h%init (gsize1, lsize1, start1)
  call io_h%openw ('data/h.dat')
  do rec=1,nrec
    f=rec
    h=rec
    call io_f%write (f, rec)
    call io_h%write (h, rec)
    if (.not.io_f%test()) call io_f%wait
    if (.not.io_h%test()) call io_h%wait
  end do
  call io_f%close
  call io_h%close

  call io_f%openr ('data/test.dat')
  call io_h%openr ('data/test.dat')
  do rec=1,nrec
    call io_f%read (f, rec)
    call io_h%read (h, rec)
    if (.not.io_f%test()) call io_f%wait
    if (.not.io_h%test()) call io_h%wait
  end do
  call io_f%close
END PROGRAM
