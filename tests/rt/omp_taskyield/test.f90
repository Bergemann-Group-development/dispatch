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
    integer(8):: n, nyield
    real(8):: time=0.0, dtime=1.0
    type(task_t), pointer:: nbor => null()
    type(task_t), pointer:: next => null()
  contains
    procedure:: init
    procedure:: update
  end type
  type(task_t), pointer, save:: head => null()
  type:: thread_t
    integer(8):: nyield=0, nupdate=0
  end type
  type(thread_t), dimension(:), pointer:: threads
  integer, save:: nx=5, ny=5, ntask
  integer(8), save:: nyield=0
  real, save:: rn=3e8
  character(len=64), save:: format='(i10,f12.6,2x,a20,3i5,3x,2f6.2)'
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
  write (*,format) omp_get_thread_num(),wallclock(),'   init: task', self%id
END SUBROUTINE

!===============================================================================
!> Update a task.  Tasks with iy=1 have no nbor, and hence can update at once.
!> Other tasks have to wait until their nbor (at iy-1) has updated, before they
!> can update.
!===============================================================================
SUBROUTINE update (self)
  class(task_t):: self
  integer(8):: i, n, thread
  !-----------------------------------------------------------------------------
  if (associated(self%nbor)) then
    write (*,format) omp_get_thread_num(),wallclock(),' update begin: task', &
      self%id, self%nbor%id, self%iy, self%time, self%nbor%time
    self%nyield = 0
    self%suspend = omp_get_thread_num()
    thread = omp_get_thread_num()
    do while (self%time >= self%nbor%time) 
      !write (*,format) omp_get_thread_num(),wallclock(),' update yield: task', self%id
      !$omp atomic
      self%nyield = self%nyield+1
      !$omp atomic
      threads(thread)%nyield = threads(thread)%nyield+1
      !$omp taskyield
    end do
    self%resume = omp_get_thread_num()
    !if (self%resume /= self%suspend) &
    !  print *,'task',self%id,' suspended by',self%suspend,' resumed by',self%resume, &
    !  ' after',self%nyield
  else
    write (*,format) omp_get_thread_num(),wallclock(),' update begin: task', self%id
  end if
  self%n = 0
  n = rn
  do i=1,n
    self%n = self%n+i
  end do
  self%time = self%time + self%dtime
  thread = omp_get_thread_num()
  write (*,format) thread,wallclock(),' update   end: task', &
    self%id, self%ix, self%iy, self%time
  !$omp atomic
  ntask = ntask+1
  threads(thread)%nupdate = threads(thread)%nupdate+1
END SUBROUTINE update
END MODULE task_mod

!===============================================================================
!> Test program
!===============================================================================
PROGRAM test
  USE task_mod
  implicit none
  type(task_t), pointer:: task, prev => null(), tmp
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
  ! task, and task%next points to the next one.
  !-----------------------------------------------------------------------------
  do ix=1,nx
  do iy=1,ny
     allocate (task)
     if (ix==1 .and. iy==1) head => task
     if (associated(prev)) prev%next => task
     call task%init (ix,iy)
     write (*,format) omp_get_thread_num(),wallclock(),'created task', task%id, task%ix, task%iy
     if (iy > 1) then
       task%nbor => prev
     end if
     prev => task
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Start an OMP task for each task_t instance
  !-----------------------------------------------------------------------------
  do iter=1,5
    task => head
    do while (associated(task))
       tmp => task
       !$omp task firstprivate (tmp)
       call tmp%update
       !$omp end task
       task => task%next
    end do
    !xomp taskwait
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
  print '(a20,i7)'  ,'number of tasks   =', ntask
  print '(a20,1p,10e10.3)'  ,'number of yields  =', real(threads%nyield)
  print '(a20,20i5)'  ,'number of updates =', threads%nupdate
  !$omp end master
  !$omp end parallel
END PROGRAM

