!> =============================================================================
!> Test program for mpi_counter.  This test works fine with 64 ranks and 6
!> threads per rank, but fails (the ranks never finish) for 128 ranks (and
!> also for 80 ranks).  It shows that -- until we understand why and can make
!> the mpi_counter_t data type entirely robust, we cannot rely on it.
!> =============================================================================
PROGRAM test_mpi_counter
  USE omp_lib
  USE io_mod
  !USE io_unit_mod
  USE mpi_mod
  USE omp_mod
  USE omp_timer_mod
  implicit none
  type(mpi_counter_t):: mpi_counter, remain
  integer(4):: ni, nc, i, j, j1, k, thread
  real(8):: time
  character(16):: file
  logical:: ok
  namelist /test_params/ ni
  !-----------------------------------------------------------------------------
  ! Initialize MPI and create a counter
  !-----------------------------------------------------------------------------
  call mpi%init
  call omp%init
  call io%init
  ni = 100
  rewind (io%input)
  read (io%input,test_params)
  call remain%init (0)
  call mpi%barrier ('begin')
  j = remain%update(1)
  write (stderr,*) wallclock(), mpi%rank, 'begins', j
  !$omp parallel default(none) shared(stderr,k,ni,mpi,mpi_counter,remain,nc) &
  !$omp private(i,j,file)
  !$omp single
  call mpi_counter%init (0)
  write (file,'("tmp.",i4.4)') 10+mpi%rank
  open (unit=10+mpi%rank,file=file,form='formatted',status='unknown')
  nc = omp_get_num_threads()
  !$omp end single
  !-----------------------------------------------------------------------------
  ! Each thread updates the counter ni time
  !-----------------------------------------------------------------------------
  do i=1,ni
    j = mpi_counter%update(1)
    write (10+mpi%rank,*) wallclock(), omp%thread, i, j
    if (mpi%rank==0 .and. omp%thread==0 .and. mod(i,ni/10)==0) then
      print *,i,omp_get_thread_num(),omp%thread
    end if
  end do
  !$omp end parallel
  close (10+mpi%rank)
  j = remain%update(-1)
  write (stderr,*) wallclock(), mpi%rank, 'finished', j
  call mpi%barrier ('end')
  if (mpi%rank==0) then
    write (stderr,*) 'OMP threads:', nc
    write (stderr,*) '  MPI ranks:', mpi%size
    write (stderr,*) ' iterations:', ni
    write (stderr,*) 'analyzing ...'
    ok = .true.
    call system ('cat tmp.* | sort -n -k 4 -u > out.dat')
    open (unit=9,file='out.dat',form='formatted',status='old')
    j1 = 1
    do k=1,ni*mpi%size*nc
      read (9,*) time, thread, i, j
      if (j /= j1) then
        ok = .false.
      else
        j1 = j1+1
      end if
      !write (stderr,*) time, thread, i, j, j1, k, ok
    end do
    if (ok) then
      write (stderr,*) 'no errors found'
      call system ('rm out.dat')
    end if
    call system ('rm tmp.*')
  end if
  call mpi%end
END PROGRAM
