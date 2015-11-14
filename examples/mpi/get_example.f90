!*******************************************************************************
PROGRAM get_example
  USE mpi_rma
  USE io_mod
  implicit none
  integer, parameter:: nw=100
  real, dimension(nw):: a, b
  integer i, win, rank, offset
  character(len=mch), save:: id= &
    'rma_example.f90 $Id$'
  type(mpi_rma_t):: mp
!-------------------------------------------------------------------------------
  call mp%init
  call print_id (id)
  call mp%barrier ('init')
!-------------------------------------------------------------------------------
  call mp%win_create_real (win, a, nw)                                          ! make a available
  do i=1,nw
    a(i) = i + nw*mp%rank
  end do
!-------------------------------------------------------------------------------
  rank = mod(mp%rank+1, mp%size)                                                ! get from one rank up
  offset = 10
  call mp%win_fence (win)                                                       ! <-  NO change to ..
  !! FIXME call mp%win_get_reals (win, b, nw-offset, offset, rank)                           !   | .. a in this ..
  call mp%win_fence (win)                                                       ! <-  .. interval
  print *,mp%rank,b(1)                                                          ! check results
!-------------------------------------------------------------------------------
  call mp%win_free (win)
  call mp%end
END PROGRAM get_example
