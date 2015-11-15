PROGRAM io_example
  USE mpi_mod
  USE mpi_coord_mod
  USE mpi_par_io_mod
  implicit none
  integer, dimension(3):: dd, d, offset, n_mpi
  real, allocatable, dimension(:,:,:):: f
  logical:: ok=.true.
  integer:: rec, rank
  real:: f1
!-------------------------------------------------------------------------------
! Initialize MPI with a cartesian MPI arrangement
!-------------------------------------------------------------------------------
  call mpi_par_io%init (mpi_dims=n_mpi)
  dd = [32,32,1]                                                                ! global dimensions
  call mpi_par_io%init (dims=dd)
  n_mpi = mpi_coord%dims
!-------------------------------------------------------------------------------
! Open file for write-read and write unique values
!-------------------------------------------------------------------------------
  d = dd/n_mpi                                                                  ! our size
  offset = mpi_par_io%coord*d                                                   ! our offset
  allocate (f(d(1),d(2),d(3)))                                                  ! allocate array
  call mpi_par_io%file_openw ('tmp.dat', dd, d, offset)                         ! open for writing
  do rec=0,9                                                                    ! write 10 records
    f = mpi%rank + 100*rec                                                      ! unique values
    call mpi_par_io%file_write (f, d, rec)                                      ! collective write
  end do
  call mpi_par_io%file_close
!-------------------------------------------------------------------------------
! Re-open file for reading and check unique values
!-------------------------------------------------------------------------------
  call mpi_par_io%file_openr ('tmp.dat', dd, d, offset)                         ! open the same file
  do rec=0,9                                                                    ! 10 records
    call mpi_par_io%file_read (f, d, rec)                                       ! read one record
    f1 = mpi%rank + 100*rec                                                     ! ok value
    ok = ok .and. f(1,1,1)==f1                                                  ! check if ok values
   !print'(a,i4,2f9.1,l3)','mpi_rank, f =', mpi_rank, f1, f(1,1,1), ok          ! check unique values
  end do
  print'(a,i4,4x,a)', ' mpi%rank', mpi%rank, merge('OK ','NOK',ok)
  call mpi_par_io%file_close
!-------------------------------------------------------------------------------
  deallocate (f)
  call mpi_par_io%end
END PROGRAM io_example
