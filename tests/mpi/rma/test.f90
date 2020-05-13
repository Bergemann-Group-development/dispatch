!===============================================================================
!> Test program for remote memory access via MPI_Put, MPI_Get
!===============================================================================
PROGRAM test_mpi_io
  USE io_mod
  USE mpi_mod
  USE mpi_rma_mod
  USE omp_mod
  USE omp_timer_mod
  implicit none
  integer:: w
  type(mpi_rma_t):: win
  integer, parameter:: nr=2, ni=4, master=0
  integer:: i(ni), iget(ni), target, k, n
  real, save:: r(nr)=0.
  real(8):: wc
  !.............................................................................
  r = (30+mpi%rank)
  !-----------------------------------------------------------------------------
  ! Initialize MPI and create a window
  !-----------------------------------------------------------------------------
  call mpi%init
  target = mod(mpi%rank+1,mpi%size)
  call win%create (nr=nr)

  !-----------------------------------------------------------------------------
  ! Set a value in the window variable, make sure all ranks are in sync.
  ! Access the value from the next rank.  In this example a barrier is needed to
  ! make sure the target of the get operation has actually stored the value before
  ! the get call.  In a real application, getting the instantaneous value is
  ! exactly the point of use MPI_Get
  !-----------------------------------------------------------------------------
  call win%lock_x    (mpi%rank)
  win%r = (40+mpi%rank)
  call win%unlock    (mpi%rank)
  call win%sync
  call mpi%barrier
  call win%lock      (target)
  call win%get       (target)
  call win%unlock    (target)
  print '(a,i3,2x,8f6.1)', 'a', mpi%rank, win%r
  call mpi%barrier (delay=0.25)

  !-----------------------------------------------------------------------------
  ! Same thing, but w/o the barrier, to show that some values were changed before
  ! the get, others were not
  !-----------------------------------------------------------------------------
  call win%lock_x    (mpi%rank)
  r = (140+mpi%rank)
  call win%unlock    (mpi%rank)
  call win%lock      (target)
  call win%get       (target)
  call win%unlock    (target)
  print '(a,i3,2x,8f6.1)', 'b', mpi%rank, win%r
  call mpi%barrier (delay=0.25)

  !-----------------------------------------------------------------------------
  ! Same thing, but w/o the barrier, and with a delay of one rank
  !-----------------------------------------------------------------------------
  if (mpi%rank==1) then
    call omp%delay (2.)
    call win%lock_x    (mpi%rank)
    r = (240+mpi%rank)
    call win%unlock    (mpi%rank)
  end if
  call win%lock      (target)
  call win%get       (target)
  call win%unlock    (target)
  print '(a,i3,2x,8f6.1)', 'c', mpi%rank, win%r
  call mpi%barrier (delay=0.25)

  !-----------------------------------------------------------------------------
  ! Write to, and get, the value from the next rank. Note that none of the
  ! RMA calls are collective and blocking, hence, in this example, the only
  ! way to make sure the RMA call is complete, before the print statement, is
  ! to add a barrier.
  !-----------------------------------------------------------------------------
  win%r = [0.,1.] + mpi%rank*3
  call mpi%barrier
  call win%lock_x (target)
  call win%put    (target)
  call win%unlock (target)
  !call win%sync
  call mpi%barrier
  print '(a,i3,2x,8f6.1)', 'd', mpi%rank, win%r
  call mpi%barrier (delay=0.25)

  !-----------------------------------------------------------------------------
  ! Write to, and get, the value from the next rank, but only from rank 1,
  ! demonstrating that lock/get/put/unlock are not collective calls. The flush
  ! call makes sure the values to be printed are flushed back to main memory,
  ! from possibly having been updated only in cache / registers.
  !-----------------------------------------------------------------------------
  win%r = [0.,1.] + mpi%rank*3 +50
  call mpi%barrier
  if (mpi%rank==1) then
    call win%lock_x (target)
    call win%put    (target)
    call win%unlock (target)
  end if
  call mpi%barrier
  print '(a,i3,2x,8f6.1)', 'e', mpi%rank, win%r
  call mpi%barrier (delay=0.25)

  !-----------------------------------------------------------------------------
  ! Speed measurement of lock/unlock pair.
  !-----------------------------------------------------------------------------
  n = 1000000
  wc = wallclock()
  if (mpi%rank==1) then
    do k=1,n
      call win%lock    (target)
      call win%unlock  (target)
    end do
  end if
  call win%lock   (target)
  call win%get    (target)
  call win%unlock (target)
  print '(a,i3,2x,8f6.1)', 'f', mpi%rank, (wallclock()-wc)*1e6/n
  call mpi%barrier (delay=0.25)

  !-----------------------------------------------------------------------------
  ! Free the window and close MPI
  !-----------------------------------------------------------------------------
  call win%free
  call mpi%end
END PROGRAM
