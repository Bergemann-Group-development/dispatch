!===============================================================================
!> Doubly linked list (DLL), carrying anything, as simply as possible
!===============================================================================
MODULE dll_mod
  USE omp_lib
  USE trace_mod
  USE io_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE random_mod
  implicit none
  private
  type, public:: dll_node_t
    class(dll_node_t), pointer:: prev => null()
    class(dll_node_t), pointer:: next => null()
    type(lock_t):: lock
    class(*), pointer:: car => null()
    integer:: id = -1
  contains
    procedure:: init => init_node
  end type
  type, public:: dll_t
    class(dll_node_t), pointer:: head => null()
    class(dll_node_t), pointer:: tail => null()
    type(lock_t):: lock
    integer:: n = 0
    character(len=32):: name='dll'
    logical:: initialized = .false.
  contains
    procedure:: init
    procedure:: append
    procedure:: prepend
    procedure:: insert_before
    procedure:: insert_after
    procedure:: remove
    procedure:: purge
    procedure:: pop
    procedure:: test1
    procedure:: test2
    procedure:: test3
    procedure, private:: check_initialized
  end type
  logical:: local = .false.
  integer:: verbose = 0
CONTAINS

!===============================================================================
!> Initialize a list node
!===============================================================================
SUBROUTINE init_node (self)
  class(dll_node_t):: self
  integer, save:: id=-1
  !-----------------------------------------------------------------------------
  if (self%id==-1) then
    self%id = id
    !$omp atomic update
     id = id+1
    call self%lock%init
  end if
END SUBROUTINE init_node

!===============================================================================
!> Print node info
!===============================================================================
SUBROUTINE print_node (self, label)
  class(dll_node_t), pointer:: self
  character(len=*):: label
  !-----------------------------------------------------------------------------
  if (associated(self)) then
    if (local) call self%lock%set ('dll_node_t%print_node')
      print *,trim(label)//'%id:',self%id
      if (associated(self%prev)) then
        print *,trim(label)//'%prev%id:',self%prev%id
      else
        print *,trim(label)//'%prev%id: null'
      end if
      if (associated(self%next)) then
        print *,trim(label)//'%next%id:',self%next%id
      else
        print *,trim(label)//'%next%id: null'
      end if
    if (local) call self%lock%unset ('dll_node_t%print_node')
  end if
END SUBROUTINE print_node

!===============================================================================
!> Initialize a doubly-linked-list (DLL), with anonymous carry.
!>
!> Note that, iIn order to be able to use the lock immidiately, it must be a
!> static part of the data type.
!>
!> To simplify handling, head and tail nodes are permanently present.
!===============================================================================
SUBROUTINE init (self, name)
  class(dll_t):: self
  character(len=*), optional:: name
  !.............................................................................
  call trace%begin ('dll_t%init')
  call self%lock%init
  call self%lock%set ('dll_t%init')
  if (present(name)) then
    self%name = name
  else
    self%name = 'dll'
  end if
  allocate (self%head, self%tail)
  call self%head%init
  call self%tail%init
  self%head%next => self%tail
  self%tail%prev => self%head
  self%n = 0
  self%initialized = .true.
  call self%lock%unset ('dll_t%init')
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Add an item at the tail.  The locking of two nodes is potentially dangerous
!> since it could cause a deadlock if the two locks were engaged in a different
!> order by another thread.
!===============================================================================
SUBROUTINE append (self, new)
  class(dll_t):: self
  class(dll_node_t), pointer:: new
  !.............................................................................
  call trace%begin ('dll_t%append')
  if (local) then
    call self%check_initialized
    call self%tail%lock%set ('dll_t%append')
    call       new%lock%set ('dll_t%append')
      new%next => self%tail
      new%prev => self%tail%prev
      new%next%prev => new
      new%prev%next => new
    call       new%lock%unset ('dll_t%append')
    call self%tail%lock%unset ('dll_t%append')
    associate (self_n => self%n)
    !$omp atomic update
    self_n = self_n + 1
    end associate
  else
    call self%lock%set ('dll_t%append')
      call self%check_initialized
      new%next => self%tail
      new%prev => self%tail%prev
      new%next%prev => new
      new%prev%next => new
      self%n = self%n+1
    call self%lock%unset ('dll_t%append')
  end if
  call trace%end ()
END SUBROUTINE append

!===============================================================================
!> Add an item at the head
!===============================================================================
SUBROUTINE prepend (self, new)
  class(dll_t):: self
  class(dll_node_t), pointer:: new
  !.............................................................................
  call trace%begin ('dll_t%prepend')
  if (local) then
    call self%check_initialized
    call self%head%lock%set ('dll_t%prepend')
    call       new%lock%set ('dll_t%prepend')
      new%prev => self%head
      new%next => self%head%next
      new%next%prev => new
      new%prev%next => new
    call       new%lock%unset ('dll_t%prepend')
    call self%head%lock%unset ('dll_t%prepend')
    associate (self_n => self%n)
    !$omp atomic update
    self_n = self_n + 1
    end associate
  else
    call self%lock%set ('dll_t%prepend')
      call self%check_initialized
      new%prev => self%head
      new%next => self%head%next
      new%next%prev => new
      new%prev%next => new
      self%n = self%n+1
    call self%lock%unset ('dll_t%prepend')
  end if
  call trace%end ()
END SUBROUTINE prepend

!===============================================================================
!> Insert a new item before an old one
!> before:
!>   H--next->O--next->T
!>   H<-prev--O<-prev--T
!> after:
!>   H--next->N--next->O--next->T
!>   H<-prev--N<-prev--O<-prev--T
!===============================================================================
SUBROUTINE insert_before (self, old, new)
  class(dll_t):: self
  class(dll_node_t), pointer:: old, new
  !.............................................................................
  call trace%begin ('dll_t%insert_before')
  if (local) then
    call self%check_initialized
    call old%lock%set ('dll_t%insert_before')
    call new%lock%set ('dll_t%insert_before')
    new%next => old
    new%prev => old%prev
    new%prev%next => new
    new%next%prev => new
    call new%lock%unset ('dll_t%insert_before')
    call old%lock%unset ('dll_t%insert_before')
    associate (self_n => self%n)
    !$omp atomic update
    self_n = self_n + 1
    end associate
  else
    call self%check_initialized
    new%next => old
    new%prev => old%prev
    new%prev%next => new
    new%next%prev => new
    self%n = self%n + 1
  end if
  call trace%end ()
END SUBROUTINE insert_before

!===============================================================================
!> Insert a new item after an old one
!> before:
!>   H--next->O--next---------->T
!>   H<-prev--O<----------prev--T
!> after:
!>   H--next->O--next->N--next->T
!>   H<-prev--O<-prev--N<-prev--T
!===============================================================================
SUBROUTINE insert_after (self, old, new)
  class(dll_t):: self
  class(dll_node_t), pointer:: old, new
  !.............................................................................
  call trace%begin ('dll_t%insert_after')
  if (local) then
    call self%check_initialized
    call old%lock%set ('dll_t%insert_after')
    call new%lock%set ('dll_t%insert_before')
    new%prev => old
    new%next => old%next
    new%prev%next => new
    new%next%prev => new
    call new%lock%unset ('dll_t%insert_before')
    call old%lock%unset ('dll_t%insert_before')
    associate (self_n => self%n)
    !$omp atomic update
    self_n = self_n + 1
    end associate
  else
    call self%check_initialized
    new%prev => old
    new%next => old%next
    new%prev%next => new
    new%next%prev => new
    self%n = self%n + 1
  end if
  call trace%end ()
END SUBROUTINE insert_after

!===============================================================================
!> Remove an item from the DLL; note that the item is not deallocated!
!> before:
!>   H--next->P--next->O--next->N--next->T
!>   H<-prev--P<-prev--O<-prev--N<-prev--T
!> after:
!>   H--next->P--next---------->N--next->T
!>   H<-prev--P<----------prev--N<-prev--T
!===============================================================================
SUBROUTINE remove (self, old)
  class(dll_t):: self
  class(dll_node_t), pointer:: old
  !.............................................................................
  class(dll_node_t), pointer:: rem, prev, next
  integer:: n
  !-----------------------------------------------------------------------------
  call trace%begin ('dll_t%remove')
  if (local) then
    call self%check_initialized
    call self%lock%set ('dll_t%remove')
    n = self%n
    self%n = self%n-1
    call self%lock%unset ('dll_t%remove')
    if (n > 0) then
      !-------------------------------------------------------------------------
      ! Note that it is necessary to cache the pointer to old, since it is the
      ! same as (not only pointing to the same as) next%prev and prev%next !!
      !-------------------------------------------------------------------------
      rem => old
      call rem%lock%set ('dll_t%remove')
        prev => rem%prev
        next => rem%next
        if (associated(prev).and.associated(next)) then
          if (verbose > 1) then
            call print_node (rem, 'rem(before)')
            call print_node (prev, 'prev(before)')
            call print_node (next, 'next(before)')
          end if
          call next%lock%set ('dll_t%remove')
            next%prev => prev
          call next%lock%unset ('dll_t%remove')
          call prev%lock%set ('dll_t%remove')
            prev%next => next
          call prev%lock%unset ('dll_t%remove')
          nullify (rem%next, rem%prev)
          if (verbose > 1) then
            call print_node (rem, 'rem(after)')
            call print_node (prev, 'prev(after)')
            call print_node (next, 'next(after)')
          end if
        else if (verbose >= 0) then
          print *, rem%id, ' appears to have already been removed'
        end if
      call rem%lock%unset ('dll_t%remove')
    end if
  else
    call self%lock%set ('dll_t%remove')
    call self%check_initialized
    if (self%n > 0) then
      !---------------------------------------------------------------------------
      ! Note that it is necessary to cache the pointer to old, since it is the
      ! same as (not only pointing to the same as) next%prev and prev%next !!
      !---------------------------------------------------------------------------
      rem => old
      prev => rem%prev
      next => rem%next
      if (associated(prev).and.associated(next)) then
        if (verbose > 1) then
          call print_node (rem, 'remo(before)')
          call print_node (prev, 'prev(before)')
          call print_node (next, 'next(before)')
        end if
        next%prev => prev
        prev%next => next
        nullify (rem%next, rem%prev)
        if (verbose > 1) then
          call print_node (rem, 'remo(after)')
          call print_node (prev, 'prev(after)')
          call print_node (next, 'next(after)')
        end if
      else if (verbose >= 0) then
        print *, rem%id, ' appears to have already been removed'
      end if
      self%n = self%n - 1
    end if
    call self%lock%unset ('dll_t%remove')
  end if
  call trace%end ()
END SUBROUTINE remove

!===============================================================================
!> Pop the last appended item off the list.  The list node and its carry should
!> be deallocated after use.
!===============================================================================
FUNCTION pop (self)
  class(dll_t):: self
  class(dll_node_t), pointer:: pop
  !.............................................................................
  call trace%begin ('dll_t%pop')
  if (local) then
    call self%check_initialized
    call self%tail%lock%set ('dll_t%pop')
      pop => self%tail%prev
      call self%remove (pop)
    call self%tail%lock%unset ('dll_t%pop')
  else
    call self%lock%set ('dll_t%pop')
    call self%check_initialized
    pop => self%tail%prev
    call self%remove (pop)
    call self%lock%unset ('dll_t%pop')
  end if
  call trace%end ()
END FUNCTION pop

!===============================================================================
!> Delete a double-linked list, including carried content
!===============================================================================
SUBROUTINE purge (self)
  class(dll_t):: self
  !.............................................................................
  class (dll_node_t), pointer:: item, next
  !-----------------------------------------------------------------------------
  call trace%begin ('dll_t%purge')
  call self%check_initialized
  call self%lock%set ('dll_t%purge')
  item => self%head
  do while (associated(item))
    next => item%next
    if (associated(item%car)) deallocate (item%car)
    deallocate (item)
    item => next
  end do
  call self%init (name=self%name)
  call self%lock%unset ('dll_t%purge')
  call trace%end ()
END SUBROUTINE purge

!===============================================================================
!> Check that the dll has been initialized
!===============================================================================
SUBROUTINE check_initialized (self)
  class(dll_t):: self
  !-----------------------------------------------------------------------------
  if (.not. self%initialized) then
    !all io%abort ('dll_t instance not initialized')
    stop 'not initialized'
  end if
END SUBROUTINE check_initialized

!===============================================================================
!> Run a self test
!===============================================================================
SUBROUTINE test1  (self)
  class(dll_t):: self
  !.............................................................................
  class(dll_node_t), pointer:: next, p, pp(:)
  real(8):: start, now
  real, pointer:: rp
  type(random_t):: random
  integer:: i, n=1000
  real(8):: used
  !-----------------------------------------------------------------------------
  print *,'==================== Timing ========================================'
  call self%init
  io%do_trace = .true.
  !-----------------------------------------------------------------------------
  ! Allocate
  !-----------------------------------------------------------------------------
  used = wallclock()
  allocate (pp(n))
  used = wallclock()-used
  print '(6p,f7.3,a)',used/n,' mus/alloc'
  !-----------------------------------------------------------------------------
  ! Initialize
  !-----------------------------------------------------------------------------
  used = wallclock()
  do i=1,n
    call pp(i)%init
  end do
  used = wallclock()-used
  print '(6p,f7.3,a)',used/n,' mus/init'
  !-----------------------------------------------------------------------------
  ! Append
  !-----------------------------------------------------------------------------
  used = wallclock()
  do i=1,n
    next => pp(i)
    call self%append (next)
  end do
  used = wallclock()-used
  print '(6p,f7.3,a)',used/n,' mus/append'
  !-----------------------------------------------------------------------------
  ! Count
  !-----------------------------------------------------------------------------
  i = 0
  used = wallclock()
  next => self%head%next
  do while (.not.associated(next,self%tail))
    i = i+1
    next => next%next
  end do
  used = wallclock()-used
  print '(6p,f7.3,a)',used/n,' mus/count'
  print *,i,' counted'
  !-----------------------------------------------------------------------------
  print *,'------------ remove head%next 1st time -------------'
  if (verbose > 0) call print_node (self%head%next, 'head%next')
  call self%remove (self%head%next)
  if (verbose > 0) call print_node (self%head%next, 'head%next')
  print *,' n =',self%n
  print *,'------------ remove head%next 2nd time -------------'
  print *,'remove head%next 1st time'
  call self%remove (self%head%next)
  if (verbose > 0) call print_node (self%head%next, 'head%next')
  print *,' n =',self%n
  print *,'------------ remove head%next 3rd time -------------'
  call self%remove (self%head%next)
  if (verbose > 0) call print_node (self%head%next, 'head%next')
  print *,' n =',self%n
  print *,'------------ remove pp(5) --------------------------'
  p => pp(5)
    if (verbose > 0) call print_node (p, 'pp(5)')
  call self%remove (p)
    if (verbose > 0) call print_node (p, 'pp(5)')
  p => self%head%next
  if (verbose > 0) call print_node (self%head%next, 'head%next')
  next => pp(6); call self%remove (next)
  next => pp(n); call self%remove (next)
  print *,' n =',self%n
  next => pp(1); call self%prepend (next)
  next => pp(2); call self%prepend (next)
  next => pp(3); call self%prepend (next)
  next => pp(6); call self%prepend (next)
  next => pp(n); call self%prepend (next)
  print *,' n =',self%n
  !-----------------------------------------------------------------------------
  ! Time remove -- pp(5) already missing, on purpose
  !-----------------------------------------------------------------------------
  used = wallclock()
  do i=1,n
    next => pp(i)
    call self%remove (next)
  end do
  used = wallclock()-used
  print '(6p,f7.3,a)',used/n,' mus/remove'
  print *,' n =',self%n
END SUBROUTINE test1

!===============================================================================
!> Run a self test
!===============================================================================
SUBROUTINE test2 (self)
  class(dll_t):: self
  !.............................................................................
  class(dll_node_t), pointer:: next
  real(8):: start, now
  real:: p
  real, pointer:: rp
  type(random_t):: random
  integer:: i
  !-----------------------------------------------------------------------------
  ! Append five numbers
  !-----------------------------------------------------------------------------
  print *,'==================== Doubly linked list test ======================='
  call self%init
  call random%init
  start = wallclock()
  now = start
  do while (now < start+2)
    p = random%ran1()
    if (p < 0.5) then
      allocate (rp)
      rp = p
      allocate (next)
      next%car => rp
      call self%prepend (next)
    else if (self%n > 0) then
      call self%remove (self%tail%prev)
      if (self%n > 0) then
        call self%remove (self%head%next)
      end if
    end if
    now = wallclock()
    !print *, now, self%n
  end do
  !-----------------------------------------------------------------------------
  ! List them
  !-----------------------------------------------------------------------------
  print *, 'list:', self%n
  i = 0
  next => self%head%next
  do while (.not.associated(next,self%tail))
    i = i+1
    next => next%next
  end do
  print *,i,' counted'
END SUBROUTINE test2

!===============================================================================
!> Run a self test
!===============================================================================
SUBROUTINE test3 (self)
  class(dll_t):: self
  !.............................................................................
  class(dll_node_t), pointer:: next
  real(8):: start, now, time
  real:: p
  real, pointer:: rp
  type(random_t):: random
  integer:: i
  integer, save:: n_threads = 2
  !-----------------------------------------------------------------------------
  ! Append five numbers
  !-----------------------------------------------------------------------------
  print *,'==================== Multi-threaded DLL test ======================='
  call self%init
  call random%init
  start = wallclock()
  now = start
  call omp_set_num_threads (n_threads)
  print *, 'n_threads  =', n_threads
  print *, 'local lock =', local
  !$omp parallel &
  !$omp private(p,rp,next,time) shared(now,start) firstprivate(random)
  do while (now < start+2)
    !---------------------------------------------------------------------------
    ! Each thread looks up a private random number
    !---------------------------------------------------------------------------
    p = random%ran1()
    if (p < 0.5) then
      !-------------------------------------------------------------------------
      ! If less than 0.5 it allocates a node, and adds it to the list, using
      ! either local or global locking
      !-------------------------------------------------------------------------
      allocate (rp)
      rp = p
      allocate (next)
      next%car => rp
      if (local) then
        call self%prepend (next)
      else
        call self%lock%set
        !print *, omp_get_thread_num(), self%n
        call self%append (next)
        !print *, omp_get_thread_num(), self%n
        call self%lock%unset
      end if
    else if (local) then
      call self%tail%prev%lock%set
      !if (self%n > 0) call self%remove (self%tail%prev)
      call self%tail%prev%lock%unset
      call self%head%next%lock%set
      !if (self%n > 0) call self%remove (self%head%next)
      call self%head%next%lock%unset
    else
      call self%lock%set ('test3')
      !if (self%n > 0) call self%remove (self%tail%prev)
      !if (self%n > 0) call self%remove (self%head%next)
      call self%lock%unset ('test3')
    end if
    time = wallclock()
    !$omp atomic write
    now = time
    !print *, now, self%n
  end do
  !$omp end parallel
  !-----------------------------------------------------------------------------
  ! List them
  !-----------------------------------------------------------------------------
  print *, 'list:', self%n
  i = 0
  next => self%head%next
  do while (.not.associated(next,self%tail))
    i = i+1
    next => next%next
  end do
  print *,i,' counted'
END SUBROUTINE test3

END MODULE dll_mod
