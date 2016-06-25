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
    procedure:: search
    procedure:: dealloc
    procedure:: print
  end type
  integer, public, parameter:: nit=10000, nadd=200, maxthread=20
  integer, public, save:: id=0
  real(8):: time(nit*nadd*maxthread)
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
!> 
!> With this version, most of the queue_by_time will involve (and lock) links
!> that are not being accessed by other threads, and when other threads are in
!> fact accessing the same link, the omp locks only affect those threads.
!>
!> This is a great time saving, relative to having a critical region for each
!> call to queue_by_time. On the Lenovo T430s adding up to 200 patches takes 
!> 4.3 mus/patch, which is to be divided by at least 2**12=4096 to get the cost
!> per cell, which is thus of the order 1 ns/cell.
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
    end if
    !---------------------------------------------------------------------------
    ! This search is valid, even if we don't lock the links.  At some point
    ! it finds an interval it belongs to, in the sense that prev%task%time is
    ! smaller and next%task%time is larger.  Since we only lock after this
    ! happens there is no guarantee that prev%next still points to the same
    ! link, and that prev%next_time%task%time is still larger.  But we really
    ! should not care about being totally fine grained in the queue, so as long
    ! as we get a lock on the prev link then we can insert the task there,
    ! ignoring the fact that another link may have been inserted in the mean
    ! time.  It is, however, important to use a lock here, since we must make
    ! sure that no other thread changes prev%next_time between the first and
    ! second statement. It is also important to flush before the lock and after
    ! the unlock, to make sure we are picking up the most recent prev%next time.
    !---------------------------------------------------------------------------
    !$omp flush
    if (next%task%time > this%task%time) then
      if (associated(prev)) then
        !$omp flush
        call prev%lock%set ('queue_by_time')
        !-----------------------------------------------------------------------
        ! We now have a "prev" that for sure has a task%time smaller than this.
        ! We cannot be sure that the "next" we have actually is the same as
        ! prev%next_time, nor can we be sure that prev%next_time%task%time is
        ! larger than this%time.
        !-----------------------------------------------------------------------
        this%next_time => prev%next_time
        prev%next_time => this
        !print*,omp%mythread(),' inserts',this%task%time,' before',next%task%time
        call prev%lock%unset ('queue_by_time')
        !$omp flush
      else
        !-----------------------------------------------------------------------
        ! When prev is not associated the link should be the new queue head.
        !-----------------------------------------------------------------------
        !$omp flush
        call self%queue%lock%set ('queue_by_time')
        this%next_time => self%queue
        self%queue => this
        !print*,omp%mythread(),' inserts',this%task%time,' before',next%task%time,' at head'
        call this%next_time%lock%unset ('queue_by_time')
        !$omp flush
      end if
      !$omp atomic
      self%n = self%n+1
      return
    end if
    !---------------------------------------------------------------------------
    ! Since atomic does not work in this context, we make a flush to ensure
    ! we have up-to-date pointers, and if prev is available we update next,
    ! in case prev%next_time has been changed in the mean time.
    !---------------------------------------------------------------------------
    !$omp flush
    if (associated(prev)) &
      next => prev%next_time
    prev => next
    next => next%next_time
  end do
  !-----------------------------------------------------------------------------
  ! We get here if we fail to find a task whose task%time is larger than the
  ! this%time.  In that case, we should append this to the list, which in 
  ! principle means this%next_time should be null, and prev%next should point
  ! to this.  This is tricky, however, since prev%next time may have changed
  ! in the mean time, pointing to some new link added by another thread. This
  ! link is then added in place, even if that is not the tail.
  !-----------------------------------------------------------------------------
  if (associated(prev)) then
    !$omp flush
    call prev%lock%set ('queue_by_time')
    next => prev%next_time
    prev%next_time => this
    this%next_time => next
    if (associated(next)) then
      !print*,omp%mythread(),' inserts',this%task%time,' before',next%task%time,' near tail'
    else
      !print*,omp%mythread(),' inserts',this%task%time,' at tail'
    end if
    call prev%lock%unset ('queue_by_time')
    !$omp flush
  else
    !$omp flush
    if (associated(self%queue)) then
      !$omp flush
      call self%queue%lock%set ('queue_by_time')
      next => self%queue%next_time
      this%next_time => next
      self%queue => this
      if (associated(next)) then
        !print*,omp%mythread(),' inserts',this%task%time,' before',next%task%time,' at head'
      else
        !print*,omp%mythread(),' inserts',this%task%time,' as the first one, at head(2)'
      end if
      call self%queue%lock%unset ('queue_by_time')
      !$omp flush
    else
      nullify (next, this%next_time)
      self%queue => this
      !print*,omp%mythread(),' inserts',this%task%time,' as the first one, at head'
      !$omp flush
    end if
  end if
  !$omp atomic
  self%n = self%n+1
END SUBROUTINE queue_by_time

!> =============================================================================
!> =============================================================================
SUBROUTINE count (self)
  class(list_t):: self
  type(link_t), pointer:: next
  integer:: i, n
  !$omp flush
  next => self%queue
  n = 0
  do while (associated(next))
    n = n+1
    next => next%next_time
  end do
  !$omp flush
  if (n == self%n) then
    !print*, n, self%n
  else
    print*, n, self%n, ' inconsistent'
    do i=1,id
      call self%search (time(i))
    end do
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
FUNCTION make_link (self, t, i) RESULT (link)
  class(list_t):: self
  real(8):: t
  integer:: i
  class(link_t), pointer:: link
  type(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  allocate (link, task)
  link%task => task
  task%time = t
  !$omp atomic
  id = id+1
  task%id = id
  time(id) = t
  call link%lock%init ('test')
  !print*,id,omp%mythread(),t,i
END FUNCTION

!> =============================================================================
!> =============================================================================
SUBROUTINE search (self, t)
  class(list_t):: self
  class(link_t), pointer:: next
  real(8):: t
  next => self%queue
  do while (associated(next))
    if (next%task%time == t) then
      return
    end if
    next => next%next_time
  end do
  print*, 'missing:', t
END SUBROUTINE search

!> =============================================================================
!> =============================================================================
SUBROUTINE dealloc (self)
  class(list_t):: self
  class(link_t), pointer:: next, old
  real(8):: t
  next => self%queue
  do while (associated(next))
    old => next
    next => next%next_time
    deallocate (old%task)
    deallocate (old)
  end do
  nullify (self%queue)
  self%n = 0
END SUBROUTINE dealloc

END MODULE list_mod

!> =============================================================================
!> =============================================================================
PROGRAM test_queue_by_time
  USE list_mod
  USE random_mod
  USE omp_timer_mod
  USE omp_mod
  implicit none
  type(list_t):: list
  class(link_t), pointer:: link
  type(task_t), pointer:: task
  integer:: i, j, n, nth, nthreads, add_per_thread
  real(8):: t, wc
  !
  call omp%init
  nthreads = omp%nthreads
  add_per_thread = nadd/nthreads
  print*,nthreads,nadd,add_per_thread
  wc = wallclock()
  do j=1,nit
    n = 0
    id = 0
    nth = 0
    !$omp parallel private(i,t)
    do i=1,add_per_thread
      t = random%uniform_1(10.)
      call list%queue_by_time (list%make_link(t,i))
      !$omp atomic
      n = n+1
    end do
    !!$omp critical
    !nth = nth+1
    !print*,omp%mythread(),nth,list%n,n
    !!$omp end critical
    !$omp end parallel
    call list%count
    call list%dealloc
  end do
  print*,n, 1e6*(wallclock()-wc)/(nit*nadd)
END PROGRAM test_queue_by_time
