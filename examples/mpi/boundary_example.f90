!*******************************************************************************
PROGRAM boundary_example
  USE mpi_coords
  implicit none
  integer, save:: n_mpi(3) = (/1,1,1/), gn(3) = (/64,64,1/)
  integer:: rank
  type(mpi_coord_t):: mp
!-------------------------------------------------------------------------------
  n_mpi(1) = 2
  n_mpi(2) = mp%size/n_mpi(1)
  call mp%init

  if (mp%rank==0) print*,'Finding the nearest neighbors from mpi_dn & mpi_up:'

  call mp%barrier('1',delay=1)
  do rank=0,mp%size-1
    if (rank==mp%rank) &
      print'(a,i4,3(5x,a,3i4))', &
        'rank:', mp%rank, &
        'pos:' , mp%coord, &
        'dn:'  , mp%dn, &
        'up:'  , mp%up
    call mp%barrier('2',delay=1)
  end do
  call mp%end
END PROGRAM boundary_example
