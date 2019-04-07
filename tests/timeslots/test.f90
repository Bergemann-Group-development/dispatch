!===============================================================================
!> This test demonstrates how to access the time slot information atomically,
!> and how -- apparently -- one does not need to protect changing the %next
!> pointer in a task list (such an assignment also cannot be preceeded by an
!> OMP atomic directive).
!>
!> The main point of test is to demonstrate that once an nbor has acquired the
!> time slot index array and the times, the most recent part of that remains
!> valid, and accessible without locking, until it is eventially overwritten.
!===============================================================================
MODULE task_mod
  implicit none
  private
  type, public:: task_t
    integer:: id=0
    integer:: nt=5
    integer:: it=1
    integer:: new=2
    integer:: nit=0
    real(8):: time=0d0, dtime=1d0, average
    real(8), allocatable:: t(:), dt(:)
  contains
    procedure:: init
    procedure:: update
    procedure:: slots
    procedure:: times
  end type
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self)
  class(task_t):: self
  !-----------------------------------------------------------------------------
  allocate (self%t (self%nt))
  allocate (self%dt(self%nt))
  self%t   = self%time
  self%dt  = self%dtime
  self%it  = 1
  self%new = 2
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE update (self)
  class(task_t):: self
  integer:: prv, new, nt
  real(8):: time, dtime
  !-----------------------------------------------------------------------------
  !$omp atomic read
  time = self%time
  !$omp atomic read
  dtime = self%dtime
  !$omp atomic read
  nt = self%nt
  !$omp atomic read
  new = self%new
  associate (stime=>self%time, st=>self%t(new), sdt=>self%dt(new), sit=>self%it, &
             sdtime=>self%dtime, snew=>self%new)
  !$omp atomic update
  time = time + dtime
  !$omp atomic write
  stime = time
  !$omp atomic write
  sdtime = dtime
  !$omp atomic write
  st = time
  !$omp atomic write
  sdt = dtime
  !$omp atomic write
  sit = new
  new = 1 + mod(new,nt)
  !$omp atomic write
  snew = new
  end associate
  !-----------------------------------------------------------------------------
  ! Compute the average, and update the update count
  !-----------------------------------------------------------------------------
  associate (average=>self%average)
  !$omp atomic write
  average = sum(self%t(self%slots()))/(self%nt-1)
  self%nit = self%nit+1
  end associate
END SUBROUTINE update

!===============================================================================
!> Slot indices, giving the most recently updated time slots. Result should be:
!>   it = 1  iit =  1  1  1  1  1
!>   it = 2  iit =  1  1  1  1  2
!>   it = 3  iit =  1  1  1  2  3
!>   it = 4  iit =  1  1  2  3  4
!>   it = 5  iit =  1  2  3  4  5
!>   it = 6  iit =  2  3  4  5  1
!>   it = 7  iit =  3  4  5  1  2
!> As this function shows, one can contstruct the time slot address array from
!> a knowledge of only the integer "it", which points to the most current time
!> slot.  This integer changes (atomically) once per time step for each task.
!> Each piece of time slot info, such as the time, task_t%t(:), or the actual
!> 3D data, %mem(), changes only once every nt time steps, so unless access is
!> very late (due to nbors with very different time steps), access should not
!> need to be protected (could be contingent on time slot info). Of course, it
!> is safer to lock access to %mem() during task update, and have nbor lock
!> their access during download of guard zone data.  It is natural to let that
!> lock be the one in the task_t data type.
!===============================================================================
FUNCTION slots (self)
  class(task_t):: self
  integer:: slots(self%nt-1)
  integer:: j, it, nt
  !-----------------------------------------------------------------------------
  !$omp atomic read
  nt = self%nt
  !$omp atomic read
  it = self%it
  do j=1,nt-1
    slots(j) = 1 + mod(j+it,nt)
    if (it < nt) then
      slots(1:nt-it) = 1
    end if
  end do
END FUNCTION slots

!===============================================================================
!> Publically avaiable times, produced by a function that does not require atomic
!> read access (except if some thread has time to update twice). The reason is
!> that only access to the time in slots (1:nt-1) is asked for, and an active
!> task changes the time in slot (nt).  If it has time to update once more after
!> this function has returned its values, it will be changing the time and mem
!> in slot number (1), which is only of interest if the accessing thread is much
!> delayed in time, relative to the nbor.  That can only happen if the accessing
!> thread is taking a very much larger time step than the nbor, and the grace
!> parameter allows it to export guard zone data extrapolated several nbor time
!> steps ahead in time -- this should probably be avoided in any case.
!===============================================================================
FUNCTION times (self)
  class(task_t):: self
  real(8):: times(1:self%nt-1)
  !.............................................................................
  times = self%t(self%slots())
END FUNCTION times

END MODULE task_mod

!===============================================================================
!===============================================================================
PROGRAM test
  USE task_mod
  USE omp_lock_mod
  USE omp_timer_mod
  implicit none
  type link_t
    type(link_t), pointer:: next => null()
    type(task_t), pointer:: task
    type(lock_t):: lock
  end type
  type(task_t):: task
  type(link_t), target:: head
  type(lock_t):: tail_lock
  type(link_t), pointer:: link, next, link2, tail
  integer:: ntask=100, it, nit=0, use_lock
  real(8):: total, thread_total, max_total=1d8, average, wait, start, now
  !-----------------------------------------------------------------------------
  ! Demonstrate retrieving memory slot indices from task%slots()
  !-----------------------------------------------------------------------------
  call task%init
  do it=1,10
    print *, it, task%slots()
    call task%update
  end do
  do use_lock=0,1
    !-----------------------------------------------------------------------------
    ! Create a task list with length ntask
    !-----------------------------------------------------------------------------
    link => head
    do it=1,ntask
      allocate (next)
      allocate (next%task)
      call next%task%init
      next%task%id = it
      link%next => next
      link => next
    end do
    tail => next
    !-----------------------------------------------------------------------------
    ! Loop over the list, and do whatever -- including deleting the tasks
    !-----------------------------------------------------------------------------
    start = wallclock()
    wait = 0d0
    total = 0d0
    nit = 0
    !$omp parallel private(thread_total, link, link2, next, average)
    do while (total < max_total)
      !---------------------------------------------------------------------------
      ! Lock the list and pop the first task
      !---------------------------------------------------------------------------
      now = wallclock()
      call head%lock%set
        wait = wait + (wallclock()-now)
        link => head%next
        head%next => head%next%next
        nit = nit+1
      call head%lock%unset
      !---------------------------------------------------------------------------
      ! Each thread computes the sum of the averages of the thread times. In the
      ! mean time, task are being popped from the list, and appended at the end.
      ! Hence %next should in principle be guaranteed to exist.  No link is
      ! removed in the parallel loop, so the only possible problem is that the
      ! tail%next value might be under change while being accessed.  However, this
      ! seems to never give a real problem, and using a lock here slows down the
      ! process very significantly
      !---------------------------------------------------------------------------
      call link%task%update
      thread_total = link%task%average
      link2 => head%next
      do while (associated(link2))
        if (use_lock > 0) then
          now = wallclock()
          call tail_lock%set
          wait = wait + (wallclock()-now)
        end if
        next => link2%next
        if (use_lock > 0) then
          call tail_lock%unset
        end if
        !$omp atomic read
        average = link%task%average
        thread_total = thread_total + average
        link2 => next
      end do
      !$omp atomic update
      total = total + thread_total
      !---------------------------------------------------------------------------
      ! Append the task at the end.  This needs to protected by a lock, to ensure
      ! that two threads don't modify tail at the same time, but the lock does not
      ! need to be the same that is used at the head of the list.
      !---------------------------------------------------------------------------
      now = wallclock()
      call tail_lock%set
        wait = wait + (wallclock()-now)
        nullify (link%next)                 ! make sure list always has an end
        tail%next => link                   ! this appends the link task
        tail => link                        ! redefine tail, for the next append
      call tail_lock%unset
      !---------------------------------------------------------------------------
    end do
    !$omp end parallel
    !-----------------------------------------------------------------------------
    ! Loop over the list, deleting the tasks
    !-----------------------------------------------------------------------------
    link => head%next
    do while (associated(link))
      !print *, link%task%id, link%task%nit, link%task%average
      next => link%next
      deallocate (link%task)
      !deallocate (link)
      link => next
    end do
    print *,nit, wait, wallclock()-start, &
      ' seconds '//merge('with','w/o ',use_lock>0)//' access lock'
  end do
END PROGRAM test
