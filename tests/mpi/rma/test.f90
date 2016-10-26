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
  integer:: i(ni), iget(ni), rank
  real, save:: r(nr)=0., rget(nr)=0., rput(nr)=0.
  !.............................................................................
  r = (30+mpi%rank)
  !-----------------------------------------------------------------------------
  ! Initialize MPI and create a window
  !-----------------------------------------------------------------------------
  call mpi%init
  call win%create_real (w, r, nr)
  !-----------------------------------------------------------------------------
  ! Set a value in the window variable, make sure it is synced
  !-----------------------------------------------------------------------------
  rank = mod(mpi%rank+1,mpi%size)
  r = (40+mpi%rank)
  call win%sync      (w)
  !-----------------------------------------------------------------------------
  ! Access the value from the next rank
  !-----------------------------------------------------------------------------
  call win%lock      (w, rank)
  call win%get_real  (w, rank, rget, nr, 0)
  call win%unlock    (w, rank)
  print '(a,i3,2x,8f5.1)', 'a', mpi%rank, r, rget
  call mpi%barrier (delay=0.2)
  !-----------------------------------------------------------------------------
  ! Write to, and get, the value from the next rank
  !-----------------------------------------------------------------------------
  rput = [0.,1.] + mpi%rank*3
  call win%lock      (w, rank)
  call win%put_real  (w, rank, rput, nr, 0)
  call win%get_real  (w, rank, rget, nr, 0)
  call win%sync      (w)
  call win%unlock    (w, rank)
  call win%sync      (w)
  print '(a,i3,2x,8f5.1)', 'b', mpi%rank, r, rget, rput
  !-----------------------------------------------------------------------------
  ! Free the window and close MPI
  !-----------------------------------------------------------------------------
  call win%free      (w)
  call mpi%end
END PROGRAM
