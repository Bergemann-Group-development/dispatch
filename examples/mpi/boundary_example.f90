!*******************************************************************************
PROGRAM boundary_example
  USE mpi_coord_mod
  implicit none
  integer, save:: mpi_dims(3) = [1,1,1], dims(3) = [64,64,1]
  integer:: rank
  type(mpi_coord_t):: mpi
!-------------------------------------------------------------------------------
  call mpi%mpi_t%init
  if ((mpi%size/2)*2 /= mpi%size) &
    call mpi%end ('ERROR: MPI size must be multiple of 2')
  mpi_dims(1) = 2
  mpi_dims(2) = mpi%size/mpi_dims(1)
  call mpi%init (mpi_dims=mpi_dims)
  if (mpi%ok) then
    if (mpi%rank==0) print*,'Find the nearest neighbors from mpi%dn & mpi%up:'
    call mpi%print
  end if
  call mpi%end
END PROGRAM boundary_example
