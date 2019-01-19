!> =============================================================================
!> Test program for remote memory access MPI_Get_accumulate
!> =============================================================================
PROGRAM test
  USE io_mod
  USE mpi_mod
  USE mpi_rma_mod
  USE omp_mod
  USE omp_timer_mod
  implicit none
  integer:: w
  type(mpi_rma_t):: win
  integer, parameter:: ni=1, nt=10
  integer:: i(ni), j(ni), m(nt), target, it, rank
  real(8):: wc
  !-----------------------------------------------------------------------------
  ! Initialize MPI, create a window and lock (shared)
  !-----------------------------------------------------------------------------
  call mpi%init
  call mpi%barrier (delay=0.1)
  target = 0
  i = (100+mpi%rank)
  j = 1
  call win%create_int (w, i, ni)
  call win%lock       (w, target)
  !-----------------------------------------------------------------------------
  ! Repeatedly add 1 to the window variable, from all ranks
  !-----------------------------------------------------------------------------
  do it=1,nt
    call win%add_int  (w, target, j, i)
    m(it) = i(1)
  end do
  !-----------------------------------------------------------------------------
  ! Unlock, free the window, and close MPI
  !-----------------------------------------------------------------------------
  call win%unlock     (w, target)
  call win%free       (w)
  !-----------------------------------------------------------------------------
  ! Sync printing
  !-----------------------------------------------------------------------------
  call mpi%barrier
  do rank=0,mpi%size-1
    if (rank==mpi%rank) then
      do it=1,nt
        print *, rank, m(it)
      end do
    end if
    call mpi%barrier
  end do
  call mpi%end
END PROGRAM
