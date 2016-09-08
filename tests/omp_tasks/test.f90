!> *****************************************************************************
!> Module with list handling for generic class task_t objects
PROGRAM test
  USE mpi_mod
  USE omp_mod
  USE omp_timer_mod
  implicit none
  integer:: i, j
  !-----------------------------------------------------------------------------
  call mpi%init
  call omp%init
  !$omp parallel
  !$omp single
  do i=0,10
    if (i==mpi%rank) then
      j = i
      !$omp task
      print*,mpi%rank,omp%thread,j,' start'
      call delay (real(mpi%rank+2))
      print*,mpi%rank,omp%thread,j,' end'
      !$omp end task
    end if
  end do
  !$omp end single
  print*,mpi%rank,omp%thread,' done'
  !$omp end parallel
  call mpi%end
CONTAINS
  !-----------------------------------------------------------------------------
  subroutine delay (s)
  real(4):: s
  real(8):: wt
  !-----------------------------------------------------------------------------
  wt = wallclock()
  do while ((wallclock()-wt) < s)
  end do
  end subroutine
END PROGRAM
