! $Id$
!===============================================================================
!> The lock module uses nested locks, to allow versatile use of locks, where
!> a procedure may want to make sure to lock a data type, even though it may
!> (or may not) already have been locked by another procedure in the calling
!> hierarchy.
!>
!> NOTE: Enabling trace_mod here, via the TRACE macro, causes a dependency 
!> loop, and can only be compiled if/when the code is already successfully 
!> compiled.  It should only be used to trace problems with locks.
!===============================================================================
MODULE omp_lock_mod
  USE io_mod
  USE io_unit_mod
  USE omp_mod
  USE omp_timer_mod
  USE omp_lib
!define TRACE
#ifdef TRACE 
  USE trace_mod
#endif
  implicit none
  private
  integer, parameter:: mkinds=20, maxproc=100
  type, public:: lock_t
    integer:: id=-1
    integer:: verbose=0
    integer(kind=omp_nest_lock_kind):: lock
    integer:: thread=-1
    character(len=4):: kind = 'orig'
    integer:: level=0
    logical:: on
    logical:: initialized=.false.
    real(8), pointer::  wait_time(:) => null()    ! time waiting for lock
    real(8), pointer::  hold_time(:) => null()    ! time holding lock
    real(8), pointer::  start_time(:) => null()   ! start time
    integer(8), pointer:: n_hold(:) => null()     ! number of times holding lock
    real(8):: time(maxproc) = 0d0                 ! lock waiting times
    logical:: links=.false.             ! lock nbor lists
    logical:: tasks=.false.             ! lock task memory updates
    class(lock_t), pointer:: next => null()
  contains
    procedure:: init
    procedure, nopass:: check_pointer
    procedure, private:: prepend
    procedure:: set
    procedure:: test
    procedure:: check
    procedure:: unset
    procedure, nopass:: destroy
    procedure:: info
  end type
  real, save:: verbose_time=0.0
  integer, save:: id_save=-1, n_lock=0
  character(len=6):: log='lock'
  character(len=6):: ext, ext1, ext2
  integer:: unit
  real(8):: alternate_time=0d0, alternate_next=0d0
  type(lock_t), save, public:: omp_lock
CONTAINS

!===============================================================================
!> Initialize a lock, taking care not to do it more than once
!===============================================================================
SUBROUTINE init (self, kind, id)
  class(lock_t), target:: self
  integer, optional:: id
  character(len=*), optional:: kind
  character(len=64):: filename
  logical omp_in_parallel
  logical, save:: tasks=.true., links=.false.
  integer, save:: verbose=0
  integer:: iostat, i
  logical:: set_id
  logical, save:: on = .true.
  logical, save:: first_time=.true.
  namelist /lock_params/ on, verbose, verbose_time, log, tasks, links, &
    alternate_time
  !-----------------------------------------------------------------------------
#ifdef TRACE 
  call trace%begin('lock_t%init')
#endif
  if (.not.self%initialized) then
    call io%lock%set
    if (.not.self%initialized) then
      if (first_time) then
        rewind (io_unit%input)
        read (io_unit%input, lock_params, iostat=iostat)
        write (io_unit%output, lock_params)
        unit = io_unit%lock1
        if (alternate_time > 0d0) then
          filename = trim(io_unit%rankbase)//'.lock1'
          open (io_unit%lock1, file=trim(filename), form='formatted', &
                status='unknown')
          filename = trim(io_unit%rankbase)//'.lock2'
          open (io_unit%lock2, file=trim(filename), form='formatted', &
                status='unknown')
          alternate_next = wallclock() + alternate_time
        else
          if (trim(log)=='queue') then
            unit = io_unit%queue
          else if (trim(log)=='log') then
            unit = io_unit%log
          else
            unit = io_unit%lock1
          end if
          filename = trim(io_unit%rankbase)//'.lock'
          open (io_unit%lock1, file=trim(filename), form='formatted', &
                status='unknown')
        end if
        omp_lock%verbose = verbose
        omp_lock%tasks = tasks
        omp_lock%links = links
        if (omp%nthreads == 1) then
          on = .false.
        end if
        omp_lock%on = on
        first_time = .false.
      end if
      call omp_init_nest_lock (self%lock)
      if (present(id)) then
        self%id = id
        id_save = max(id_save,id)
        set_id = .true.
      else
        id_save = id_save+1
        self%id = id_save
        set_id = .false.
      end if
      if (present(kind)) then
        self%kind = kind(1:4)
      else
        self%kind = 'void'
      end if
      self%on  = on
      allocate (self%start_time(0:omp%nthreads-1))
      allocate (self%wait_time (0:omp%nthreads-1))
      allocate (self%hold_time (0:omp%nthreads-1))
      allocate (self%n_hold    (0:omp%nthreads-1))
      self%start_time = 0d0
      self%wait_time = 0d0
      self%hold_time = 0d0
      self%n_hold = 0
      if (n_lock > 1) then
        call self%prepend
      end if
      n_lock = n_lock+1
      if (wallclock() < verbose_time .or. omp_lock%verbose > 1) &
        write(unit,*) omp_get_thread_num(), 'initialized lock kind, id ', &
          self%kind, self%id, set_id
      self%initialized = .true.
    end if
    call io%lock%unset
  end if
#ifdef TRACE 
  call trace%end()
#endif
END SUBROUTINE init

!===============================================================================
!> Prepend a lock to a simple linked list.  The first lock added will inherit
!> the null pointer that ends the list, from the original omp_lock%next.
!===============================================================================
SUBROUTINE prepend (self)
  class(lock_t), target:: self
  !-----------------------------------------------------------------------------
  if (.not.omp_lock%on) return
  self%next => omp_lock%next
  omp_lock%next => self
END SUBROUTINE prepend

!===============================================================================
!> Check a pointer to a lock
!===============================================================================
SUBROUTINE check_pointer (lock, label, id)
  type(lock_t), pointer:: lock
  character(len=*):: label
  integer:: id
  !-----------------------------------------------------------------------------
  if (.not. associated(lock)) then
    write (io_unit%log,*) label//' ERROR: lock not assoc', id
    allocate (lock)
    call lock%init ('err1')
  end if
  if (.not. lock%initialized) then
    write (io_unit%log,*) label//' ERROR: lock not initialized', id
    allocate (lock)
    call lock%init ('err2')
  end if
  if (.not. associated(lock%start_time)) then
    write (io_unit%log,*) label//' ERROR: start_time not associated', id
    lock%initialized = .false.
    call lock%init ('err3')
  end if
END SUBROUTINE check_pointer

!===============================================================================
!> Return with the lock, waiting if necessary, and creating the lock if needed
!===============================================================================
SUBROUTINE set (self, label, itimer)
  class(lock_t):: self
  character(len=*), optional:: label
  integer, optional:: itimer
  character(len=32):: filename
  real(8):: start, end, used
  integer:: thread
  logical:: log
  !-----------------------------------------------------------------------------
  if (.not.omp_lock%on) return
  if (.not.self%initialized) call self%init
#ifdef TRACE 
  call trace%begin('lock_t%set')
#endif
  !-----------------------------------------------------------------------------
  ! Must print log info before the lock, to reveal dead-locks
  !-----------------------------------------------------------------------------
  call alternate_log (self)
  start = wallclock()
  log = start < verbose_time .or. omp_lock%verbose > 1
  if (log) then
    if (present(label)) then
      write (unit,'(f10.6,a,a4,i6,i3,i4,2x,a)') start, &
        ' lock kind, id, level, thread = ', self%kind, self%id, &
        self%level, omp%thread, '   get at '//label
    else
      write (unit,'(f10.6,a,a4,i6,i3,i4,2x,a)') start, &
        ' lock kind, id, level, thread = ', self%kind, self%id, &
        self%level, omp%thread, '   get'
    end if
    !flush (unit)
  end if
  start = wallclock()
  !$omp flush
  call omp_set_nest_lock (self%lock)
  !$omp flush
  end = wallclock()
  used = end-start
  thread = omp_get_thread_num()
  self%start_time(thread) = start
  self%wait_time(thread) = &
  self%wait_time(thread) + used
  self%n_hold(thread) = self%n_hold(thread) + 1
  if (present(itimer)) then
    !$omp atomic update
    omp_lock%time(itimer) = &
    omp_lock%time(itimer) + used
  end if
  self%level = self%level+1
  !$omp flush
  if (log) then
    if (present(label)) then
      write (unit,'(f10.6,a,a4,i6,i3,i4,6p,f10.1,2x,a)') end, &
        ' lock kind, id, level, thread = ', self%kind, self%id, &
        self%level, omp%thread, used, '   set at '//label
    else
      write (unit,'(f10.6,a,a4,i6,i3,i4,6p,f10.1,2x,a)') end, &
        ' lock kind, id, level, thread = ', self%kind, self%id, &
        self%level, omp%thread, used, '   set'
    end if
    !flush (unit)
  end if
#ifdef TRACE 
  call trace%end()
#endif
END SUBROUTINE set

!===============================================================================
!> Alternate between two log files
!===============================================================================
SUBROUTINE alternate_log (self)
  class(lock_t):: self
  real(8):: wc
  !-----------------------------------------------------------------------------
  if (.not.omp_lock%on) return
  if (alternate_time > 0d0) then
    wc = wallclock()
    if (wc > alternate_next) then
      write (io_unit%log,*) omp%thread,'waiting at lock_cr'
      flush (io_unit%log)
      !$omp critical (lock_cr)
      write (io_unit%log,*) omp%thread,'in lock_cr'
      if (wc > alternate_next) then
        if (unit == io_unit%lock1) then
          rewind (io_unit%lock2)
          unit = io_unit%lock2
        else
          rewind (io_unit%lock1)
          unit = io_unit%lock1
        end if
        alternate_next = wc + alternate_time
      end if
      write (io_unit%log,*) omp%thread,'exit lock_cr'
      !$omp end critical (lock_cr)
    end if
  end if
END SUBROUTINE alternate_log

!===============================================================================
!> Return true if we already have the lock, or if we get it now
!===============================================================================
FUNCTION test (self) RESULT (out)
  class(lock_t):: self
  logical:: out
  real(8):: wc
  !-----------------------------------------------------------------------------
  if (.not.omp_lock%on) return
#ifdef TRACE 
  call trace%begin('lock_t%test')
#endif
  if (omp%thread == self%thread) then
    out = .true.
    !print *, omp%thread, ' already has lock id', self%id
  else
    !print *, omp%thread, ' omp_test_lock for id', self%id
    out =  omp_test_nest_lock (self%lock)
    if (out) then
      self%thread = omp%thread
      !$omp atomic update
      self%level = self%level+1
      wc = wallclock()
      if (wc < verbose_time .or. omp_lock%verbose > 1 .or. alternate_time > 0d0) then
        write (io_unit%log,'(f10.6,i4,a,i6,i3,a)') wc, omp%thread, ' lock id', &
          self%id, self%level, '  test'
      end if
    else
      !print *, omp%thread, ' failed to acquire lock id', self%id, self%thread
    end if
  end if
  !$omp flush
#ifdef TRACE 
  call trace%end()
#endif
END FUNCTION test

!===============================================================================
!> Return true only if we have rightful access to the locked item
!===============================================================================
FUNCTION check (self) RESULT (out)
  class(lock_t):: self
  logical:: out
  !-----------------------------------------------------------------------------
  if (.not.omp_lock%on) return
  out = (omp%thread == self%thread)
END FUNCTION check

!===============================================================================
!> Unset the associated thread number and release the lock (in that order!)
!===============================================================================
SUBROUTINE unset (self, label)
  class(lock_t):: self
  character(len=*), optional:: label
  integer:: level, thread
  real(8):: wc, used
  !-----------------------------------------------------------------------------
  if (.not.omp_lock%on) return
#ifdef TRACE 
  call trace%begin('lock_t%unset')
#endif
  self%level = self%level-1
  level = self%level
  !$omp flush
  call omp_unset_nest_lock (self%lock)
  !$omp flush
  wc = wallclock()
  thread = omp_get_thread_num()
  used = wc - self%start_time(thread)
  self%hold_time(thread) = &
  self%hold_time(thread) + used
  !call alternate_log (self)
  if (wc < verbose_time .or. omp_lock%verbose > 1) then
    if (present(label)) then
      write (unit,'(f10.6,a,a4,i6,i3,i4,6p,f10.1,2x,a)') wc, &
        ' lock kind, id, level, thread = ', self%kind, self%id, level, &
        thread, used, ' unset at '//label
    else
      write (unit,'(f10.6,a,a4,i6,i3,i4,6p,f10.1,2x,a)') wc, &
        ' lock kind, id, level, thread = ', self%kind, self%id, level, &
        thread, used, ' unset'
    end if
  end if
#ifdef TRACE 
  call trace%end()
#endif
END SUBROUTINE unset

!===============================================================================
SUBROUTINE destroy (this)
  type(lock_t), pointer:: this
  !-----------------------------------------------------------------------------
  if (omp_lock%on .and. associated(this)) then
    if (this%initialized) then
      call omp_destroy_nest_lock (this%lock)
      call remove (this)
      deallocate (this%start_time)
      deallocate (this%wait_time)
      deallocate (this%hold_time)
      deallocate (this%n_hold)
    end if
    if (omp_lock%verbose > 0) then
      write (unit,'(f10.6,a,a4,i6,i3,i4,6p,f10.1,2x,a)') wallclock(), &
        ' lock kind, id, level, thread = ', this%kind, this%id, this%level, &
        omp_get_thread_num(), 0.0, ' destroy'
    end if
    deallocate (this)
    nullify (this)
  end if
END SUBROUTINE destroy

!===============================================================================
SUBROUTINE remove (this)
  type(lock_t), pointer:: this, lock, prev
  !-----------------------------------------------------------------------------
  if (.not.omp_lock%on) return
  !$omp critical (lock_cr)
  nullify (prev)
  lock => omp_lock%next
  do while (associated(lock))
    if (associated(this,lock)) then
      if (associated(prev)) then
        prev%next => lock%next
      else
        omp_lock%next => lock%next
      end if
    end if
    prev => lock
    lock => lock%next
  end do
  n_lock = n_lock-1
  !$omp end critical (lock_cr)
END SUBROUTINE remove

!===============================================================================
SUBROUTINE info (self, unit, p_verbose)
  class(lock_t):: self
  integer:: unit, p_verbose
  real(8):: time
  character(len=4):: kinds(mkinds)
  integer(8), dimension(mkinds):: n_hold
  real(8), dimension(mkinds):: wait_time, hold_time
  integer:: nkinds, i
  class(lock_t), pointer:: lock
  !-----------------------------------------------------------------------------
  if (.not.omp_lock%on) return
  !!omp critical (lock_info_cr)
  time = 0d0
  if (omp_lock%verbose > -1) then
    !---------------------------------------------------------------------------
    ! Find the kind of locks that exist
    !---------------------------------------------------------------------------
    do i=1,mkinds
      kinds(i) = '    '
    end do
    n_hold = 0
    wait_time = 0d0
    hold_time = 0d0
    lock => omp_lock%next
    do while (associated(lock))
      nkinds = 1
      kinds(1) = 'void'
      do i=1,mkinds-1
        if (lock%kind(1:4) == kinds(i)) then
          exit
        else if (kinds(i) == '    ') then
          kinds(i) = lock%kind(1:4)
          nkinds = i
          exit
        end if
      end do
      if (omp_lock%verbose > 2) &
        write (stdout,*) lock%id, lock%kind, lock%initialized, associated(lock%wait_time)
      if (.not.associated(lock%wait_time)) then
        print *,'omp_lock_t%init: ERROR',omp%thread,lock%kind
        stop
      end if
      wait_time(i) = wait_time(i) + sum(lock%wait_time)
      hold_time(i) = hold_time(i) + sum(lock%hold_time)
      n_hold(i)    = n_hold(i)    + sum(lock%n_hold)
      lock%start_time = 0.0_8
      lock%wait_time = 0.0_8
      lock%hold_time = 0.0_8
      lock%n_hold = 0_8
      lock => lock%next
    end do
    time = sum(wait_time)
    if (omp_lock%verbose >= p_verbose) then
      call write_table (stdout)
    end if
    write (stdout,'(a,1p,g11.2,a,10f8.1)') ' total lock waiting time =', &
      time, ' s'
    wait_time = 0d0
    hold_time = 0d0
    n_hold = 0
  end if
  !!omp end critical (lock_info_cr)
  return
contains
subroutine write_table (unit)
  integer:: unit
  real(8):: wait_sum, hold_sum
  write (unit,'(11x,a,a)') '  i  lock       wait      wait(%)     mus/n', &
                                     '       hold      hold(%)     mus/n           n'
  wait_sum = max(sum(wait_time(1:nkinds)),tiny(1.0))
  hold_sum = max(sum(hold_time(1:nkinds)),tiny(1.0))
  do i=1,nkinds
    if (wait_time(i) > 0d0) then
      write (unit,'(11x,i3,2x,a4,2(2x,1p,e12.3,2p,f10.1,6p,f10.1,0p),i12)') i, kinds(i), &
      wait_time(i), wait_time(i)/wait_sum, wait_time(i)/max(n_hold(i),1_8), &
      hold_time(i), hold_time(i)/hold_sum, hold_time(i)/max(n_hold(i),1_8), &
      n_hold(i)
    end if
  end do
end subroutine write_table
END SUBROUTINE info

!===============================================================================
END MODULE omp_lock_mod
