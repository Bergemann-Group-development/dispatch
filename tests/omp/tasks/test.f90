!> *****************************************************************************
PROGRAM test
  USE omp_lib
  USE mpi_mod
  USE omp_mod
  USE omp_timer_mod
  USE timer_mod
  implicit none
  integer:: i, j, n, parent, child, task
  !-----------------------------------------------------------------------------
  call mpi%init
  call omp%init
  !-----------------------------------------------------------------------------
  ! Each of the OMP_NUM_THREAD tasks makes 10 iterations, each time spawning a
  ! new task, which is costly.  The total number of tasks started is thus 10*
  ! omp_num_threads(). When the omp_num_thread() tasks are finished with the 
  ! loop (takes very little time), they wait for their 10 sub-tasks, which means
  ! their thread is suspended, and may in fact take on any of the sub-tasks.
  ! When all 10 sub-tasks are done, the initial thread wakes up and finishes.
  !-----------------------------------------------------------------------------
  n = 0
  call tic
  !$omp parallel private (parent, child)
  do i=1,10
    parent = omp_get_thread_num()
    !$omp task
    !$omp critical
    n = n+1
    task = n
    !$omp end critical
    child = omp_get_thread_num()
    !print*,mpi%rank,parent,' starts', child,' on task',task
    call delay (real(mpi%rank+2))
    child = omp_get_thread_num()
    !print*,mpi%rank,parent,'   ends', child,' on task',task
    !$omp end task
  end do
  print*,wallclock(),mpi%rank,parent,' loop done', n
  !$omp taskwait
  print*,wallclock(),mpi%rank,parent,' wait done', n
  !$omp end parallel
  call toc ('took')
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
