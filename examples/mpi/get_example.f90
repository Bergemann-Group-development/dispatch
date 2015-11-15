!*******************************************************************************
PROGRAM get_example
  USE mpi_rma_mod
  USE io_mod
  implicit none
  integer, parameter:: nw=100
  real, dimension(nw):: a, b
  integer i, win, rank, offset
  logical:: ok
  character(len=120), save:: id= &
    'get_example.f90 $Id$'
!-------------------------------------------------------------------------------
  call mpi_rma%init
  call print_id (id)
  call mpi_rma%barrier('init',delay=1)
!-------------------------------------------------------------------------------
  call mpi_rma%win_create_real (win, b, nw)                                     ! make b available
  do i=1,nw
    b(i) = i + nw*mpi_rma%rank                                                  ! put results in window
  end do
!-------------------------------------------------------------------------------
  rank = modulo(mpi_rma%rank+1, mpi_rma%size)                                   ! get from one rank up

  offset = 0
  call mpi_rma%win_fence (win)                                                  ! <-  NO access to ..
  call mpi_rma%win_get_reals (win, a(1+offset), nw-offset, offset, rank)        !   | .. a&b in this ..
  call mpi_rma%win_fence (win)                                                  ! <-  .. interval
  ok = a(offset+1) == offset+1 + nw*modulo(mpi_rma%rank+1, mpi_rma%size)
  print *,mpi_rma%rank,b(offset+1:offset+5)                                     ! check results
  print *,'The result is '//merge('ok ','nok',ok)
  call mpi_rma%barrier(delay=1)

  offset = 10
  call mpi_rma%win_fence (win)                                                  ! <-  .. interval
  call mpi_rma%win_get_reals (win, a(1+offset), nw-offset, offset, rank)        !   | .. a&b in this ..
  call mpi_rma%win_fence (win)                                                  ! <-  .. interval
  ok = a(offset+1) == offset+1 + nw*modulo(mpi_rma%rank+1, mpi_rma%size)
  print *,mpi_rma%rank,b(offset+1:offset+5)                                     ! check results
  print *,'The result is '//merge('ok ','nok',ok)
  call mpi_rma%barrier(delay=1)
!-------------------------------------------------------------------------------
  call mpi_rma%win_free (win)
  call mpi_rma%end
END PROGRAM get_example
