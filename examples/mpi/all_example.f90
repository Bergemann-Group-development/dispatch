!*******************************************************************************
PROGRAM all_example
  USE mpi_mod
  USE mpi_all_mod
  USE io_mod
  implicit none
  integer:: i, rank, ndata=10
  real, allocatable, dimension(:):: send1, recv1
  real, allocatable, dimension(:,:):: send2, recv2
  character(len=120), save:: id= &
    'all_example.f90 $Id$'
!-------------------------------------------------------------------------------
  call mpi%init
  call print_id (id)
  call mpi%print
!-------------------------------------------------------------------------------
! Allocate send and receive buffers, and fill the first send buffer
!-------------------------------------------------------------------------------
  allocate (send1(ndata), recv1(ndata))
  allocate (send2(ndata, mpi%size), recv2(ndata, mpi%size))
  print*,'send/recv buffer size:', sizeof(send2)
  do i=1,ndata
    send1(i) = i + 100*mpi%rank
    do rank=1,mpi%size
      send2(i,rank) = i+10*rank + 100*mpi%rank
    end do
  end do
!-------------------------------------------------------------------------------
! Send ndata values from rank 0 to all others
!-------------------------------------------------------------------------------
  call mpi_all%scatter_reals (send2, recv1, ndata, 0, mpi%size)
  print *,'recv1:',mpi%rank, recv1(1:ndata:ndata-1)
  send1 = recv1 + 1
!-------------------------------------------------------------------------------
! Gather ndata values from all ranks to rank 0
!-------------------------------------------------------------------------------
  call mpi_all%gather_reals (send1, recv2, ndata, 0, mpi%size)
  if (mpi%rank==0) then
    do rank=1,mpi%size
      print *,'recv2:',mpi%rank,recv2(1:ndata:ndata-1,rank),rank
    end do
  end if
!-------------------------------------------------------------------------------
! Send different data from each rank to all other ranks
!-------------------------------------------------------------------------------
  do rank=1,mpi%size
  do i=1,ndata
    send2(i,rank) = i + rank*2000
  end do
  end do
  call mpi_all%allgather_reals (send2, recv2, ndata, 0, mpi%size)
  do rank=1,mpi%size
  do i=1,ndata
    print *,'allgather recv2:',mpi%rank,rank,recv2(1:ndata:ndata-1,rank)
  end do
  end do
!-------------------------------------------------------------------------------
! Send data from all ranks to all ranks
!-------------------------------------------------------------------------------
  call mpi_all%alltoall_reals (send2, recv2, ndata, mpi%size)
  do rank=1,mpi%size
  do i=1,ndata
    print *,'alltoall recv2:',mpi%rank,rank,recv2(1:ndata:ndata-1,rank)
  end do
  end do
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  call alltoallv_reals_mpi (send, send_count, send_offset, &
!                            recv, recv_count, recv_offset, &
!                            count)
  deallocate (send1, recv1, send2, recv2)
  call mpi%end
END PROGRAM all_example
