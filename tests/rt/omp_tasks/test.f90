!===============================================================================
!> Demonstrate and test OMP task creation. The program generates a number of
!> independent task, which are then executed by all threads, including the 
!> master thread, which created the tasks.
!===============================================================================
MODULE task_mod
  USE omp_lib
  USE omp_mod
  USE omp_timer_mod
  implicit none
  type:: task_t
    integer:: id
    integer(8):: n
  contains
    procedure:: init
    procedure:: update
  end type
  real, save:: rn=3e8
  character(len=64), save:: format='(i10,1p,e14.4,2x,a20,i5,e14.4)'
CONTAINS

!===============================================================================
!> Initialize a task
!===============================================================================
SUBROUTINE init (self, id)
  class(task_t):: self
  integer:: id
  !-----------------------------------------------------------------------------
  self%id = id
  write (*,format) omp_get_thread_num(),wallclock(),'   init: task', self%id
END SUBROUTINE

!===============================================================================
!> Update a task
!===============================================================================
SUBROUTINE update (self)
  class(task_t):: self
  integer(8):: i, n
  !-----------------------------------------------------------------------------
  write (*,format) omp_get_thread_num(),wallclock(),' update begin: task', self%id
  self%n = 0
  n = rn
  do i=1,n
    self%n = self%n+i
  end do
  write (*,format) omp_get_thread_num(),wallclock(),' update   end: task', self%id, real(self%n)
END SUBROUTINE update
END MODULE task_mod

!===============================================================================
!> Test program
!===============================================================================
PROGRAM test
  USE task_mod
  implicit none
  type(task_t), pointer:: task
  integer:: i, ntask, ncore=2
  real(8):: wc
  !-----------------------------------------------------------------------------
  call omp%init
  !-----------------------------------------------------------------------------
  ! Start a parallel region, and start by using only the master thread
  !-----------------------------------------------------------------------------
  ntask = 50
  print '(a,1pe7.0)','loop length =',rn
  !$omp parallel
  !$omp master
  wc = wallclock()
  print *,'number of   cores =', ncore
  print *,'number of threads =', omp_get_num_threads()
  !-----------------------------------------------------------------------------
  ! Start ntask tasks -- one for each instance of task_t
  !-----------------------------------------------------------------------------
  do i=1,ntask
     write (*,format) omp_get_thread_num(),wallclock(),'creating task', i
     !$omp task
     allocate (task)
     call task%init (i)
     call task%update
     !$omp end task
  end do
  !-----------------------------------------------------------------------------
  ! Wait for all tasks to finish, which releases the master thread, so it can
  ! participate in the work
  !-----------------------------------------------------------------------------
  write (*,format) omp_get_thread_num(),wallclock(),' before taskwait'
  !$omp taskwait
  !-----------------------------------------------------------------------------
  ! After all work is done, we are back to using the master thread only
  !-----------------------------------------------------------------------------
  write (*,format) omp_get_thread_num(),wallclock(),'  after taskwait'
  wc = wallclock()-wc
  print '(a20,f7.3)','total time =', wc
  print '(a20,f7.3)','nanosec/iteration =', 1e9*wc*min(ncore,omp_get_num_threads())/(rn*ntask)
  !$omp end master
  !$omp end parallel
END PROGRAM

