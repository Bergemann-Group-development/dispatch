!*******************************************************************************
PROGRAM main
  USE mpi_basic
  USE mpi_reduction
  USE io_mod
  implicit none
  character(len=mch), save:: id= &
    'main.f90 $Id$'
  integer rank, i, m
  integer, parameter:: na=10
  real, dimension(na):: a
  type(mpi_reduction_t):: mp
!-------------------------------------------------------------------------------
  call mp%init
  call print_id(id)
  do i=1,na
    a(i) = i+mp%rank
  end do

  call mp%barrier('before print')
  if (master) print*,'just print'
  m = mp%rank
  call mp%max_int (m)
  call mp%max_reals (a, na)
  print *,mp%rank,mp%size,m,a(na)
  call mp%barrier('after print')

  if (master) print*,'ordered print'
  call mp%barrier('before loop')
  do rank=0,mp%size-1
    m = rank
    call mp%max_int (m)
    call mp%max_reals (a, na)
    if (rank==mp%rank) print*,mp%rank,mp%size,m,a(na)
    call mp%barrier('in loop',delay=1)
  end do

  call mp%end
END PROGRAM main
!*******************************************************************************
