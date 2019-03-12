!===============================================================================
!> Test program for remote memory access MPI_Get_accumulate
!===============================================================================
PROGRAM test
  USE io_mod
  USE mpi_mod
  USE mpi_rma_mod
  implicit none
  type(mpi_rma_t):: incr, store
  integer:: i(1)=1, nt=10, master=0, it, target
  !-----------------------------------------------------------------------------
  ! Initialize MPI, create a window for accumulation and lock it for shared 
  ! access.  Create another window for storing the results, for later printing
  !-----------------------------------------------------------------------------
  call mpi%init
  call incr%create (ni=1)
  call incr%lock (master)
  call store%create (ni=nt)
  !-----------------------------------------------------------------------------
  ! Repeatedly add 1 to the master window value, from all ranks
  !-----------------------------------------------------------------------------
  incr%i(1)=100
  do it=1,nt
    call incr%accumulate (master, i=i)
    store%i(it) = incr%i(1)
  end do
  !-----------------------------------------------------------------------------
  ! Get values on master and print.  A barrier is needed here, to ensure that
  ! the target ranks are finished with the previous step
  !-----------------------------------------------------------------------------
  call mpi%barrier
  if (mpi%rank==0) then
    do target=0,mpi%size-1
      call store%lock (target)
      call store%get (target)
      call store%unlock (target)
      do it=1,nt
        print *, target, it, store%i(it)
      end do
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Unlock, free the window, and close MPI. A barrier is needed here to ensure
  ! that the target does not free the window.
  !-----------------------------------------------------------------------------
  call mpi%barrier
  call incr%unlock (master)
  call incr%free
  call store%free
  call mpi%end
END PROGRAM
