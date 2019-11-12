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
  type, public:: lock_t
    integer:: id=-1
    integer(kind=omp_nest_lock_kind):: lock
    integer:: thread=-1
    character(len=4):: kind = 'void'
    integer:: level=0
    logical:: on
    !real:: verbose=0.0
    real(8):: wait_time=0.0_8           ! time waiting for lock
    real(8):: hold_time=0.0_8           ! time holding lock
    real(8):: start_time=0.0_8          ! start time
    real(8):: used_time=0.0_8           ! used time
    integer(8):: n_hold=0_8             ! number of times holding lock
    logical:: links=.false.             ! lock nbor lists
    logical:: tasks=.false.             ! lock task memory updates
    type(lock_t), pointer:: next => null()
  contains
    procedure:: init
    procedure:: append
    procedure:: set
    procedure:: test
    procedure:: check
    procedure:: unset
    procedure:: destroy
    procedure, nopass:: info
  end type
  integer, save:: verbose=0
  real, save:: verbose_time=0.0
  integer, save:: id_save=-1, n_lock=0
  logical, save:: on = .true.
  character(len=6):: log='queue'
  character(len=6):: ext, ext1, ext2
  integer:: unit
  real(8):: wait_time=0d0, alternate_time=0d0, alternate_next=0d0
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
  integer:: iostat
  logical:: set_id
  logical:: first_time=.true.
  namelist /lock_params/ on, verbose, verbose_time, log, tasks, links, &
    alternate_time
  !-----------------------------------------------------------------------------
  if (.not.on) return
#ifdef TRACE 
  call trace%begin('lock_t%init')
#endif
  !$omp critical (lock_cr)
  if (self%id == -1) then
    if (first_time) then
      rewind (io_unit%input)
      read (io_unit%input, lock_params, iostat=iostat)
      write (io_unit%output, lock_params)
      if (trim(log)=='queue') then
        unit = io_unit%queue
      else if (trim(log)=='log') then
        unit = io_unit%log
      else
        unit = io_unit%output
      end if
      if (alternate_time > 0d0) then
        close (unit)
        ext = '.lock1'
        ext1 = '.lock1'
        ext2 = '.lock2'
        filename = trim(io_unit%rankbase)//ext
        open (unit, file=trim(filename), form='formatted', &
              status='unknown')
        alternate_next = wallclock() + alternate_time
      end if
      first_time = .false.
    end if
    call omp_init_nest_lock (self%lock)
    omp_lock%tasks = tasks
    omp_lock%links = links
    if (present(id)) then
      self%id = id
      id_save = max(id_save,id)
      set_id = .true.
    else
      id_save = id_save+1
      self%id = id_save
      set_id = .false.
    end if
    if (present(kind)) self%kind = kind(1:4)
    if (wallclock() < verbose_time .or. verbose > 1) &
      write(unit,*) omp_get_thread_num(), 'initialized lock kind, id ', &
        self%kind, self%id, set_id
    omp_lock%on = on
  end if
  !$omp end critical (lock_cr)
  if (omp%nthreads == 1) &
    on = .false.
#ifdef TRACE 
  call trace%end()
#endif
END SUBROUTINE init

!===============================================================================
!> Append a lock to a simple linked list
!===============================================================================
SUBROUTINE append (self)
  class(lock_t), target:: self
  !-----------------------------------------------------------------------------
  if (.not.on) return
  !$omp critical (lock_cr)
  self%next => omp_lock%next
  omp_lock%next => self
  n_lock = n_lock+1
  !$omp end critical (lock_cr)
END SUBROUTINE append

!===============================================================================
!> Return with the lock, waiting if necessary, and creating the lock if needed
!===============================================================================
SUBROUTINE set (self, label)
  class(lock_t):: self
  character(len=*), optional:: label
  character(len=32):: filename
  integer:: thread, othread, level
  real(8):: wc
  !-----------------------------------------------------------------------------
  if (.not.on) return
  if (self%id == -1) call self%init
#ifdef TRACE 
  call trace%begin('lock_t%set')
#endif
  thread = omp_get_thread_num()
  !-----------------------------------------------------------------------------
  ! Must print log info before the lock, to reveal dead-locks
  !-----------------------------------------------------------------------------
  call alternate_log (self)
  wc = wallclock()
  if (wc < verbose_time .or. verbose > 1 .or. alternate_time > 0d0) then
    if (present(label)) then
      write (unit,'(f10.6,a,a4,i6,i3,i4,2x,a)') wc, &
        ' lock kind, id, level, thread = ', self%kind, self%id, &
        self%level, thread, '   get at '//label
    else
      write (unit,'(f10.6,a,a4,i6,i3,i4,2x,a)') wc, &
        ' lock kind, id, level, thread = ', self%kind, self%id, &
        self%level, thread, '   get'
    end if
  end if
  wc = wallclock()
  self%start_time = wc
  call omp_set_nest_lock (self%lock)
  wc = wallclock()
  self%used_time = wc - self%start_time
  self%wait_time = self%wait_time + self%used_time
  self%n_hold = self%n_hold + 1
  !$omp atomic update
  wait_time = wait_time + self%used_time
  self%thread = thread
  !$omp atomic update
  self%level = self%level+1
  !$omp flush
  if (wc < verbose_time .or. verbose > 1 .or. alternate_time > 0d0) then
    write (unit,'(f10.6,a,a4,i6,i3,i4,2x,a)') wc, &
      ' lock kind, id, level, thread = ', self%kind, self%id, &
      self%level, thread, '   set'
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
  if (alternate_time > 0d0) then
    wc = wallclock()
    if (wc > alternate_next) then
      ext = ext2
      ext2 = ext1
      ext1 = ext
      close (unit)
      open (unit, file=trim(io_unit%rankbase)//ext, form='formatted', &
            status='unknown')
      alternate_next = wc + alternate_time
    end if
  end if
END SUBROUTINE alternate_log

!===============================================================================
!> Return true if we already have the lock, or if we get it now
!===============================================================================
FUNCTION test (self) RESULT (out)
  class(lock_t):: self
  logical:: out
  integer:: thread
  real(8):: wc
  !-----------------------------------------------------------------------------
  if (.not.on) return
#ifdef TRACE 
  call trace%begin('lock_t%test')
#endif
  thread = omp_get_thread_num()
  !print *, thread, ' set lock id, thread', self%id, self%thread
  !$omp flush
  !print *, thread, ' testing lock id', self%id, self%thread
  if (thread == self%thread) then
    out = .true.
    !print *, thread, ' already has lock id', self%id
  else
    !print *, thread, ' omp_test_lock for id', self%id
    out =  omp_test_nest_lock (self%lock)
    if (out) then
      self%thread = thread
      !$omp atomic update
      self%level = self%level+1
      wc = wallclock()
      if (wc < verbose_time .or. verbose > 1 .or. alternate_time > 0d0) then
        write (io_unit%log,'(f10.6,i4,a,i6,i3,a)') wc, thread, ' lock id', &
          self%id, self%level, '  test'
      end if
    else
      !print *, thread, ' failed to acquire lock id', self%id, self%thread
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
  integer:: thread
  !-----------------------------------------------------------------------------
  thread = omp_get_thread_num()
  out = (thread == self%thread)
END FUNCTION check

!===============================================================================
!> Unset the associated thread number and release the lock (in that order!)
!===============================================================================
SUBROUTINE unset (self, label)
  class(lock_t):: self
  character(len=*), optional:: label
  integer:: thread, level
  real(8):: wc
  !-----------------------------------------------------------------------------
  if (.not.on) return
#ifdef TRACE 
  call trace%begin('lock_t%unset')
#endif
  thread = omp_get_thread_num()
  !$omp atomic update
  self%level = self%level-1
  level = self%level
  if (self%level==0) then
    self%thread = -1
  end if
  !$omp flush
  call omp_unset_nest_lock (self%lock)
  wc = wallclock()
  self%used_time = wc - self%start_time
  self%hold_time = self%hold_time + self%used_time
  !write(unit,*) &
  !  'lock%unset: kind, id, thread, level = ', self%kind, self%id, self%thread, self%level
  call alternate_log (self)
  if (wc < verbose_time .or. verbose > 1 .or. alternate_time > 0d0) then
    if (present(label)) then
      write (unit,'(f10.6,a,a4,i6,i3,i4,2x,a)') wc, &
        ' lock kind, id, level, thread = ', self%kind, self%id, level, thread, ' unset at '//label
    else
      write (unit,'(f10.6,a,a4,i6,i3,i4,2x,a)') wc, &
        ' lock kind, id, level, thread = ', self%kind, self%id, level, thread, ' unset'
    end if
  end if
#ifdef TRACE 
  call trace%end()
#endif
END SUBROUTINE unset

!===============================================================================
SUBROUTINE destroy (self)
  class(lock_t):: self
  !-----------------------------------------------------------------------------
  if (.not.on) return
  call omp_destroy_nest_lock (self%lock)
  self%thread = -1
END SUBROUTINE destroy

!===============================================================================
SUBROUTINE info (unit)
  integer:: unit
  real(8):: time
  type(lock_t), pointer:: lock
  !-----------------------------------------------------------------------------
  if (.not.on) return
  if (verbose > 0) then
    lock => omp_lock%next
    !$omp critical (lock_cr)
    write (unit,'(1x,a)') '    id type    wait        wait/n       hold        hold/n       n'
    do while (associated(lock))
      write (unit,'(i6,2x,a4,4f12.6,i12)') lock%id, lock%kind, &
        lock%wait_time, 1d3*lock%wait_time/max(lock%n_hold,1_8), &
        lock%hold_time, 1d3*lock%hold_time/max(lock%n_hold,1_8), &
        lock%n_hold
      lock%wait_time = 0.0_8
      lock%hold_time = 0.0_8
      lock%n_hold = 0_8
      lock => lock%next
    end do
    !$omp end critical (lock_cr)
  end if
  !$omp atomic read
  time = wait_time
  write (unit,'(a,f10.3)') ' total lock waiting time =', time
  !$omp atomic write
  wait_time = 0d0
END SUBROUTINE info

!===============================================================================
END MODULE omp_lock_mod
