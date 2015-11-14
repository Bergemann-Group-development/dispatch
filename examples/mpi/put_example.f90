!*******************************************************************************
PROGRAM put_example
  USE mpi_rma
  USE io_mod
  implicit none
  integer, parameter:: nw=100
  real, dimension(nw):: a, b
  integer i, win, rank, offset
  character(len=mch), save:: id= &
    'put_example.f90 $Id$'
  type(mpi_rma_t):: mp
!-------------------------------------------------------------------------------
  call mp%init
  call print_id (id)
  call mp%barrier('init')
!-------------------------------------------------------------------------------
  call mp%win_create_real (win, b, nw)                                          ! make b available
  do i=1,nw
    a(i) = i + nw*mp%rank
  end do
!-------------------------------------------------------------------------------
  rank = modulo(mp%rank-1, mp%size)                                             ! put to one rank down
  offset = 10
  call mp%win_fence (win)                                                       ! <-  NO access to ..
  !!! FIXME call mp%put_reals (win, a(1+offset), nw-offset, offset, rank)                 !   | .. a&b in this ..
  call mp%win_fence (win)                                                       ! <-  .. interval
  print *,mp%rank,b(1+offset)                                                   ! check results
!-------------------------------------------------------------------------------
  call mp%win_free (win)
  call mp%end
END PROGRAM put_example
