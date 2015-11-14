!*******************************************************************************
PROGRAM send_example
  USE mpi_send_recv
  USE io_mod
  implicit none
  integer, parameter:: nw=100
  real, dimension(nw):: a, b
  integer i, to, from, offset, req(2)
  character(len=mch), save:: id= &
    'rma_example.f90 $Id$'
  type(mpi_t):: mp
!-------------------------------------------------------------------------------
  call mp%init
  call print_id (id)
  call mp%barrier('init')
!-------------------------------------------------------------------------------
  do i=1,nw
    a(i) = i + nw*mp%rank
  end do
!-------------------------------------------------------------------------------
  to   = modulo(mp%rank-1, mp%size)                                           ! send to one rank down
  from = modulo(mp%rank+1, mp%size)                                           ! recv from one rank up
  offset = 10
  call isend_reals_mpi (mp, a(1+offset), nw-offset,   to, req(1))                   ! <-  NO change to ..
  call irecv_reals_mpi (mp, b          , nw-offset, from, req(2))                   !   | .. a in this ..
  call waitall_mpi (mp, req, 2)                                                     ! <-  .. interval
  print *,mp%rank,b(1)                                                         ! check results
!-------------------------------------------------------------------------------
  call mp%end
END PROGRAM send_example
