!*******************************************************************************
PROGRAM main
  USE mpi_mod
  USE mpi_reduce_mod
  USE io_mod
  implicit none
  character(len=120), save:: id= &
    'main.f90 $Id$'
  integer rank, i, m
  integer, parameter:: na=10
  real, dimension(na):: a
!-------------------------------------------------------------------------------
  call mpi%init
  call print_id(id)
  do i=1,na
    a(i) = i+mpi%rank
  end do

  call mpi%barrier('before print')
  if (mpi%master) print*,'just print'
  m = mpi%rank
  call mpi_reduce%max_int (m)
  call mpi_reduce%max_reals (a, na)
  print *,mpi%rank,mpi%size,m,a(na)
  call mpi%barrier('after print')

  if (mpi%master) print*,'ordered print'
  call mpi%barrier('before loop')
  do rank=0,mpi%size-1
    m = rank
    call mpi_reduce%max_int (m)
    call mpi_reduce%max_reals (a, na)
    if (rank==mpi%rank) print*,mpi%rank,mpi%size,m,a(na)
    call mpi%barrier('in loop',delay=1)
  end do

  call mpi%end
END PROGRAM main
!*******************************************************************************
