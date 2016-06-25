!> *****************************************************************************
!> Module with list handling for generic class task_t objects
!> *****************************************************************************
MODULE list_mod
  USE omp_mod
  USE omp_lock_mod
  implicit none
  private
  !
  type, public:: task_t
    real(8):: time=0d0
    integer:: id=0
  end type
  !
  type, public:: link_t
    class(link_t), pointer:: next_time => null()    ! next link in time order
    class(task_t), pointer:: task => null()
    type(lock_t):: lock
  end type
  !
  type, public:: list_t
    integer:: n = 0
    class(link_t), pointer:: queue => null()
  contains
    procedure:: queue_by_time
    procedure:: make_link
    procedure:: count
    procedure:: print
  end type
CONTAINS

!> =============================================================================
!> Add a link node into the queue chain, in ascending time order.  When a link
!> to a task with time larger than this%task%time is found, or the end of the
!> list is fount, this link is inserted before.
!>
!> -> link0 -next_time-> link1 -next_time-> link2
!>     |                  |                  |
!>     t0                 t1                 t2
!>        tnew0   tnew1         tnew2
!>
!> In the sketch above the worst case is when two tasks are being queued at
!> nearly the same time, into the same interval.  This may be handled by having
!> a lock on each link, so when the two tasks have identified the interval
!> they should be in, which means they want to modify the next_time of the
!> previous link (link0), they need to compete about who gets the lock. The one
!> that succeeds modifies the next_time of link0 and its own next_time and then
!> unlocks the link.
!> =============================================================================
SUBROUTINE queue_by_time (self, this)
  class(list_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: next, prev
  !-----------------------------------------------------------------------------
  nullify (prev)
  !xomp atomic read
  !$omp flush
  next => self%queue
  do while (associated(next))
    if (associated(next%task, this%task)) then
      print*, omp_mythread(), ' WARNING: task', this%task%id, ' is already on ready queue'
      return
    !---------------------------------------------------------------------------
    ! This search is valid, even if we don't lock the links.  At some point it
    ! finds an interval it belongs to, in the sense that prev%task%time is
    ! smaller and next%task%time is larger.  Since we only lock after this
    ! happens there is no guarantee that prev%next still points to the same
    ! link, and that prev%next_time%task%time is still larger.  But we really
    ! should not care about being totally fine grained in the queue, so as long
    ! as we insert the get a lock on the prev link then we can insert the task
    ! there, ignoring the fact that another link may have been inserted in the
    ! mean time.  It is, however, important to use a lock here, since we must
    ! make sure that no other thread changes prev%next_time between the first
    ! and second statement.  It is important to flush prev, to make sure we are
    ! picking up the most recent prev%next time.
    !---------------------------------------------------------------------------
    else if (next%task%time > this%task%time) then
      if (associated(prev)) then
        call prev%lock%set ('queue_by_time')
        !$omp flush
        !-----------------------------------------------------------------------
        ! We now have a "prev" that for sure has a task%time smaller than this.
        ! We cannot be sure that the "next" we have actually is the same as
        ! prev%next_time, nor can we be sure that prev%next_time%task%time is
        ! larger than this%time.
        !-----------------------------------------------------------------------
        this%next_time => prev%next_time
        prev%next_time => this
        !$omp flush
        call prev%lock%unset ('queue_by_time')
      else
        call self%queue%lock%set ('queue_by_time')
        !$omp flush
        this%next_time => self%queue%next_time
        self%queue => this
        !$omp flush
        call self%queue%lock%unset ('queue_by_time')
      end if
      self%n = self%n+1
      return
    end if
    !---------------------------------------------------------------------------
    ! This could possibly skip some recently inserted link, but this is not
    ! serious.  What is more important is to make sure we do get a unique next,
    ! by performing an atomic read.
    !---------------------------------------------------------------------------
    !xomp atomic read
    !$omp flush
    prev => next
    next => prev%next_time
  end do
  !-----------------------------------------------------------------------------
  ! We get here if we fail to find a task whose task%time is larger than the
  ! this%time.  In that case, we should append this to the list, which means
  ! this%next_time should be null, and prev%next should point to this.  This is
  ! tricky, however, since prev%next time may have changed in the mean time,
  ! pointing to some new link added by another thread.
  !-----------------------------------------------------------------------------
  if (associated(prev)) then
    call prev%lock%set ('queue_by_time')
    !$omp flush
    next => prev%next_time
    prev%next_time => this
    this%next_time => next
    !$omp flush
    call prev%lock%unset ('queue_by_time')
  else
    call this%lock%set ('queue_by_time')
    !$omp flush
    nullify (next)
    if (associated(self%queue)) next => self%queue%next_time
    self%queue => this
    this%next_time => next
    !$omp flush
    call this%lock%unset ('queue_by_time')
  end if
  self%n = self%n+1
END SUBROUTINE queue_by_time

!> =============================================================================
!> =============================================================================
SUBROUTINE count (self)
  class(list_t):: self
  type(link_t), pointer:: next
  integer:: n
  !$omp flush
  next => self%queue
  n = 0
  do while (associated(next))
    n = n+1
    next => next%next_time
  end do
  !$omp flush
  if (n /= self%n) then
    print*, 'inconsistent counts:', n, self%n
  end if
END SUBROUTINE count

!> =============================================================================
!> =============================================================================
SUBROUTINE print (self)
  class(list_t):: self
  type(link_t), pointer:: next
  integer:: i
  next => self%queue
  i = 0
  do while (associated(next))
    i = i+1
    print*,i,next%task%id,next%task%time
    next => next%next_time
  end do
END SUBROUTINE print

!> =============================================================================
!> =============================================================================
FUNCTION make_link (self, time) RESULT (link)
  class(list_t):: self
  real(4):: time
  class(link_t), pointer:: link
  type(task_t), pointer:: task
  integer, save:: id=0
  !-----------------------------------------------------------------------------
  allocate (link, task)
  link%task => task
  task%time = time
  !$omp atomic
  id = id+1
  task%id = id
  call link%lock%init ('test')
END FUNCTION
END MODULE list_mod

!> =============================================================================
!> =============================================================================
PROGRAM test_queue_by_time
  USE list_mod
  USE random_mod
  implicit none
  type(list_t):: list
  class(link_t), pointer:: link
  type(task_t), pointer:: task
  integer:: i
  !
  !$omp parallel
  do i=1,50
    call list%queue_by_time (list%make_link(random%uniform_1(10.)))
    call list%count
  end do
  !$omp end parallel
  print *, list%n
END PROGRAM test_queue_by_time
