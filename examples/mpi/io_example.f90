PROGRAM io_example
  USE mpi_coords
  USE mpi_par_io
  implicit none
  integer, dimension(3):: dd, d, offset, n_mpi
  real, allocatable, dimension(:,:,:):: f
  logical:: ok=.true.
  integer:: rec, rank
  real:: f1
  type(mpi_par_io_t):: mp
!-------------------------------------------------------------------------------
! Initialize MPI with a cartesian MPI arrangement
!-------------------------------------------------------------------------------
  call mp%init (mpi_dims=n_mpi)
  dd = [32,32,1]                                                                ! global dimensions
  call mp%init (dims=dd)
  n_mpi = mp%dims
!-------------------------------------------------------------------------------
! Open file for write-read and write unique values
!-------------------------------------------------------------------------------
  d = dd/n_mpi                                                                  ! our size
  offset = mp%coord*d                                                           ! our offset
  allocate (f(d(1),d(2),d(3)))                                                  ! allocate array
  call mp%file_openw ('tmp.dat', dd, d, offset)                                 ! open for writing
  do rec=0,9                                                                    ! write 10 records
    f = mp%rank + 100*rec                                                       ! unique values
    call mp%file_write (f, d, rec)                                              ! collective write
  end do
  call mp%file_close
!-------------------------------------------------------------------------------
! Re-open file for reading and check unique values
!-------------------------------------------------------------------------------
  call mp%file_openr ('tmp.dat', dd, d, offset)                                 ! open the same file
  do rec=0,9                                                                    ! 10 records
    call mp%file_read (f, d, rec)                                               ! read one record
    f1 = mp%rank + 100*rec                                                      ! ok value
    ok = ok .and. f(1,1,1)==f1                                                  ! check if ok values
   !print'(a,i4,2f9.1,l3)','mpi_rank, f =', mpi_rank, f1, f(1,1,1), ok          ! check unique values
  end do
  print'(a,i4,4x,a)', ' mp%rank', mp%rank, merge('OK ','NOK',ok)
  call mp%file_close
!-------------------------------------------------------------------------------
  deallocate (f)
  call mp%end
END PROGRAM io_example
