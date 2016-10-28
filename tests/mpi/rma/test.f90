!> =============================================================================
!> Test program for remote memory access via MPI_Put, MPI_Get
!> =============================================================================
PROGRAM test_mpi_io
  USE io_mod
  USE mpi_mod
  USE mpi_rma_mod
  USE omp_mod
  USE omp_timer_mod
  implicit none
  integer:: w
  type(mpi_rma_t):: win
  integer, parameter:: nr=2, ni=4
  integer:: i(ni), iget(ni), target, k, n
  real, save:: r(nr)=0., rget(nr)=0., rput(nr)=0.
  real(8):: wc
  !.............................................................................
  r = (30+mpi%rank)
  !-----------------------------------------------------------------------------
  ! Initialize MPI and create a window
  !-----------------------------------------------------------------------------
  call mpi%init
  target = mod(mpi%rank+1,mpi%size)
  call win%create_real (w, r, nr)
  call mpi%barrier (delay=0.25)

  !-----------------------------------------------------------------------------
  ! Set a value in the window variable, make sure all ranks are in sync.
  ! Access the value from the next rank.  In this example a barrier is needed to
  ! make sure the target of the get operation has actually stored the value before
  ! the get call.  In a real application, getting the instantaneous value is
  ! exactly the point of use MPI_Get
  !-----------------------------------------------------------------------------
  call win%lock_x    (w, mpi%rank)
  r = (40+mpi%rank)
  call win%unlock    (w, mpi%rank)
  call win%sync      (w)
  call mpi%barrier
  call win%lock      (w, target)
  call win%get_real  (w, target, rget, nr, 0)
  call win%unlock    (w, target)
  print '(a,i3,2x,8f6.1)', 'a', mpi%rank, r, rget
  call mpi%barrier (delay=0.25)

  !-----------------------------------------------------------------------------
  ! Same thing, but w/o the barrier, to show that some values were changed before
  ! the get, others were not
  !-----------------------------------------------------------------------------
  call win%lock_x    (w, mpi%rank)
  r = (140+mpi%rank)
  call win%unlock    (w, mpi%rank)
  call win%lock      (w, target)
  call win%get_real  (w, target, rget, nr, 0)
  call win%unlock    (w, target)
  print '(a,i3,2x,8f6.1)', 'b', mpi%rank, r, rget
  !flush (io%output)
  call mpi%barrier (delay=0.25)

  !-----------------------------------------------------------------------------
  ! Same thing, but w/o the barrier, and with a delay of one rank
  !-----------------------------------------------------------------------------
  if (mpi%rank==1) then
    call omp_timer%delay (2.)
    call win%lock_x    (w, mpi%rank)
    r = (240+mpi%rank)
    call win%unlock    (w, mpi%rank)
  end if
  call win%lock      (w, target)
  call win%get_real  (w, target, rget, nr, 0)
  call win%unlock    (w, target)
  print '(a,i3,2x,8f6.1)', 'c', mpi%rank, r, rget
  !flush (io%output)
  call mpi%barrier (delay=0.25)

  !-----------------------------------------------------------------------------
  ! Write to, and get, the value from the next rank. Note that none of the
  ! RMA calls are collective and blocking, hence, in this example, the only
  ! way to make sure the RMA call is complete, before the print statement, is
  ! to add a barrier.
  !-----------------------------------------------------------------------------
  rput = [0.,1.] + mpi%rank*3
  call win%lock      (w, target)
  call win%put_real  (w, target, rput, nr, 0)
  call win%get_real  (w, target, rget, nr, 0)
  call win%unlock    (w, target)
  call mpi%barrier
  call win%sync      (w)
  print '(a,i3,2x,8f6.1)', 'd', mpi%rank, r, rget, rput
  call mpi%barrier (delay=0.25)

  !-----------------------------------------------------------------------------
  ! Write to, and get, the value from the next rank, but only from rank 1,
  ! demonstrating that lock/get/put/unlock are not collective calls. The sync
  ! call makes sure the values to be printed are flushed back to main memory,
  ! from possibly having been updated only in cache / registers.
  !-----------------------------------------------------------------------------
  rput = [0.,1.] + mpi%rank*3 +50
  if (mpi%rank==1) then
    call win%lock      (w, target)
    call win%put_real  (w, target, rput, nr, 0)
    call win%get_real  (w, target, rget, nr, 0)
    call win%unlock    (w, target)
  end if
  call mpi%barrier
  if (mpi%rank==1) then
    call win%sync      (w)
  end if
  print '(a,i3,2x,8f6.1)', 'e', mpi%rank, r, rget, rput
  call mpi%barrier (delay=0.25)

  !-----------------------------------------------------------------------------
  ! Speed measurement of lock/unlock pair.
  !-----------------------------------------------------------------------------
  n = 100000
  wc = wallclock()
  if (mpi%rank==1) then
    do k=1,n
      call win%lock    (w, target)
      call win%unlock  (w, target)
    end do
  end if
  call win%lock      (w, target)
  call win%put_real  (w, target, rput, nr, 0)
  call win%get_real  (w, target, rget, nr, 0)
  call win%unlock    (w, target)
  print '(a,i3,2x,8f6.1)', 'f', mpi%rank, (wallclock()-wc)*1e6/n
  call mpi%barrier (delay=0.25)

  !-----------------------------------------------------------------------------
  ! Free the window and close MPI
  !-----------------------------------------------------------------------------
  call win%free      (w)
  call mpi%end
END PROGRAM
