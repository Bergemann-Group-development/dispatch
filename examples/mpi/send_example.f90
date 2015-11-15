!*******************************************************************************
PROGRAM send_example
  USE mpi_send_mod
  USE io_mod
  implicit none
  integer, parameter:: nw=100
  real, dimension(nw):: a, b
  integer i, to, from, offset, req(2)
  character(len=120), save:: id= &
    'rma_example.f90 $Id$'
!-------------------------------------------------------------------------------
  call mpi_send%init
  call print_id (id)
  call mpi_send%barrier('init')
!-------------------------------------------------------------------------------
  do i=1,nw
    a(i) = i + nw*mpi_send%rank
  end do
!-------------------------------------------------------------------------------
  to   = modulo(mpi_send%rank-1, mpi_send%size)                                 ! send to one rank down
  from = modulo(mpi_send%rank+1, mpi_send%size)                                 ! recv from one rank up
  offset = 10
  call mpi_send%isend_reals (a(1+offset), nw-offset,   to, req(1))              ! <-  NO change to ..
  call mpi_send%irecv_reals (b          , nw-offset, from, req(2))              !   | .. a in this ..
  call mpi_send%waitall (req, 2)                                                ! <-  .. interval
  print *,mpi_send%rank,b(1)                                                    ! check results
!-------------------------------------------------------------------------------
  call mpi_send%end
END PROGRAM send_example
