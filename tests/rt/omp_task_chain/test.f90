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
    integer:: ix, iy
    integer:: suspend, resume
    integer(8):: n
    real(8):: time=0.0, dtime=1.0
    type(task_t), pointer:: child => null()
    type(task_t), pointer:: next => null()
  contains
    procedure:: init
    procedure:: update
  end type
  type(task_t), pointer, save:: head => null()
  type:: thread_t
    integer(8):: nupdate=0
  end type
  type(thread_t), dimension(:), pointer:: threads
  integer, save:: nx=320, ny=320, ntask
  integer(8):: n
  real, save:: rn=3e8
  character(len=64), save:: format='(i10,f12.6,2x,a20,3i7,3x,2f6.2)'
CONTAINS

!===============================================================================
!> Initialize a task
!===============================================================================
SUBROUTINE init (self, ix, iy)
  class(task_t):: self
  integer, save:: id=0
  integer:: ix, iy
  !-----------------------------------------------------------------------------
  !$omp critical
  id = id + 1
  self%id = id
  !$omp end critical
  self%ix = ix
  self%iy = iy
END SUBROUTINE

!===============================================================================
!> Update a task. For each task, task%child points to the next task in a chain
!> to be executed by the same thread.
!===============================================================================
SUBROUTINE update (self)
  class(task_t):: self
  integer(8):: i, thread
  !-----------------------------------------------------------------------------
  self%n = 0
  n = rn*5**2/nx**2
  do i=1,n
    self%n = self%n+i
  end do
  self%time = self%time + self%dtime
  thread = omp_get_thread_num()
  !write (*,format) thread,wallclock(),' update   end: task', &
    !self%id, self%ix, self%iy, self%time
  !$omp atomic
  ntask = ntask+1
  threads(thread)%nupdate = threads(thread)%nupdate+1
  if (associated(self%child)) then
    call self%child%update
  end if
END SUBROUTINE update
END MODULE task_mod

!===============================================================================
!> Test program
!===============================================================================
PROGRAM test
  USE task_mod
  implicit none
  class(task_t), pointer:: task, prev => null(), tmp, child
  integer:: ncore=2, ix, iy, iter
  real(8):: wc
  !-----------------------------------------------------------------------------
  call omp%init
  !-----------------------------------------------------------------------------
  ! Start a parallel region, and start by using only the master thread
  !-----------------------------------------------------------------------------
  ntask = 0
  print '(a,1pe7.0)','loop length =',rn
  !$omp parallel
  !$omp master
  wc = wallclock()
  print *,'number of   cores =', ncore
  print *,'number of threads =', omp_get_num_threads()
  allocate (threads(0:omp_get_num_threads()-1))
  !-----------------------------------------------------------------------------
  ! Initialize task chain, and nbor relations.  head is a pointer to the first
  ! task, and task%next points to the next one, while task%child points to the
  ! next task in the y-direction.
  !-----------------------------------------------------------------------------
  do ix=1,nx
  do iy=1,ny
     allocate (task)
     if (ix==1 .and. iy==1) head => task
     if (associated(prev)) prev%next => task
     call task%init (ix,iy)
     if (iy > 1) then
       prev%child => task
     end if
     prev => task
  end do
  end do
  write (*,format) omp_get_thread_num(),wallclock(),'created task', task%id, task%ix, task%iy
  !-----------------------------------------------------------------------------
  ! Start an OMP task chain for each boundary task_t instance, five times.
  ! Note that is NECESSARY that new OMP task created in the do loop creates
  ! a private copy of the task pointer, since otherwise it will become modified
  ! in the next iteration, possibly / likely even before the task actually starts.
  !-----------------------------------------------------------------------------
  do iter=1,5
    task => head
    do while (associated(task))
      if (task%iy==1) then
        !$omp task firstprivate (task)
        call task%update
        !$omp end task
      end if
      task => task%next
    end do
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
  print '(a20,f7.3)','nanosec/iteration =', 1e9*wc*min(ncore,omp_get_num_threads())/(n*ntask)
  print '(a20,i7)'  ,'number of tasks   =', ntask
  print '(a20,f7.3)','millisec per task =', 1e3*wc*min(ncore,omp_get_num_threads())/ntask
  print '(a20,20i8)'  ,'number of updates =', threads%nupdate
  !$omp end master
  !$omp end parallel
END PROGRAM

