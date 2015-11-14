!*******************************************************************************
PROGRAM all_example
  USE mpi_all
  USE io_mod
  implicit none
  integer:: i, rank, ndata=100
  real, allocatable, dimension(:):: send1, recv1
  real, allocatable, dimension(:,:):: send2, recv2
  character(len=mch), save:: id= &
    'all_example.f90 $Id$'
 type(mpi_t):: mp
!-------------------------------------------------------------------------------
  call mp%init
  call print_id (id)
!-------------------------------------------------------------------------------
! Allocate send and receive buffers, and fill the first send buffer
!-------------------------------------------------------------------------------
  allocate (send1(ndata), recv1(ndata))
  allocate (send2(ndata, mp%size), recv2(ndata, mp%size))
  do i=1,ndata
  do rank=1,mp%size
     send2(i,rank) = i+100*rank
  end do
  end do
!-------------------------------------------------------------------------------
! Send ndata values from the rank 0 to all others
!-------------------------------------------------------------------------------
  call scatter_reals_mpi (mp, send2, recv1, ndata, 0, mp%size)
  print *,'recv1:',mp%rank, recv1(1:ndata:ndata-1)
  send1 = recv1 + 1
!-------------------------------------------------------------------------------
! Gather ndata values from all ranks to rank 0
!-------------------------------------------------------------------------------
  call  gather_reals_mpi (mp, send1, recv2, ndata, 0, mp%size)
  if (mp%rank==0) then
    do rank=1,mp%size
      print *,'recv2:',mp%rank,rank,recv2(1:ndata:ndata-1,rank)
    end do
  end if
!-------------------------------------------------------------------------------
! Send different data from each rank to all other ranks
!-------------------------------------------------------------------------------
  do rank=1,mp%size
  do i=1,ndata
    send2(i,rank) = i + rank*2000
  end do
  end do
  call allgather_reals_mpi (mp, send2, recv2, ndata, 0, mp%size)
  do rank=1,mp%size
  do i=1,ndata
    print *,'allgather recv2:',mp%rank,rank,recv2(1:ndata:ndata-1,rank)
  end do
  end do
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  call  alltoall_reals_mpi (mp, send2, recv2, ndata, mp%size)
  do rank=1,mp%size
  do i=1,ndata
    print *,'alltoall recv2:',mp%rank,rank,recv2(1:ndata:ndata-1,rank)
  end do
  end do
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  call alltoallv_reals_mpi (send, send_count, send_offset, &
!                            recv, recv_count, recv_offset, &
!                            count)
  deallocate (send1, recv1, send2, recv2)
  call mp%end
END PROGRAM all_example
