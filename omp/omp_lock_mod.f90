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
    real(8):: time(maxproc) = 0d0                 ! lock waiting times
    real(8):: start_time = 0d0
    logical:: links=.false.             ! lock nbor lists
    logical:: tasks=.false.             ! lock task memory updates
  contains
    procedure:: init
    procedure, nopass:: check_pointer
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
  real(8):: wait_time = 0d0
  type(lock_t), save, public:: omp_lock
CONTAINS

!===============================================================================
!> Initialize a lock, taking care not to do it more than once
!===============================================================================
SUBROUTINE init (self, kind, id)
  class(lock_t), target:: self
  character(len=*), optional:: kind
  integer, optional:: id
  !.............................................................................
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
  !$omp critical (input_cr)
  if (.not.self%initialized) then
    !---------------------------------------------------------------------------
    ! This part is executed only once, reading parameters and opening log files
    !---------------------------------------------------------------------------
    if (first_time) then
      rewind (io_unit%input)
      read (io_unit%input, lock_params, iostat=iostat)
      write (io_unit%nml, lock_params)
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
    !---------------------------------------------------------------------------
    ! The part below is executed once for all locks
    !---------------------------------------------------------------------------
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
      self%kind = kind(1:min(len(kind),4))
    else
      self%kind = 'void'
    end if
    self%on  = on
    n_lock = n_lock+1
    if (wallclock() < verbose_time .or. omp_lock%verbose > 1) &
      write(unit,*) omp_get_thread_num(), 'initialized lock kind, id ', &
        self%kind, self%id, set_id
    self%initialized = .true.
  end if
  if (verbose > 0) then
    write (io_unit%log,'(a,i8,l4)') &
      'omp_lock_t%init: id, initialized =', self%id, self%initialized
  end if
  !$omp end critical (input_cr)
#ifdef TRACE 
  call trace%end()
#endif
END SUBROUTINE init

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
#ifdef TRACE 
  call trace%begin('lock_t%set')
#endif
  if (.not. self%initialized) then
    write (stderr,'(a,i4)') &
    'omp_lock_t%set WARNING: lock not initialized, thread =', omp%thread
    write (io_unit%log,'(f12.6,2x,a,i4)') wallclock(), &
    'omp_lock_t%set WARNING: lock not initialized, thread =', omp%thread
  call self%init
  end if
  !-----------------------------------------------------------------------------
  ! Must print log info before the lock, to reveal dead-locks
  !-----------------------------------------------------------------------------
  call alternate_log (self)
  start = wallclock()
  log = start < verbose_time .or. omp_lock%verbose > 1
  if (log) then
    if (present(label)) then
      write (unit,'(f10.6,a,a4,i6,i3,i4,12x,a)') start, &
        ' lock kind, id, level, thread = ', self%kind, self%id, &
        self%level, omp_get_thread_num(), '   get at '//label
    else
      write (unit,'(f10.6,a,a4,i6,i3,i4,12x,a)') start, &
        ' lock kind, id, level, thread = ', self%kind, self%id, &
        self%level,   omp_get_thread_num(), '   get'
    end if
    flush (unit)
  end if
  start = wallclock()
  !$omp flush
  call omp_set_nest_lock (self%lock)
  !$omp flush
  end = wallclock()
  used = end-start
  self%start_time = end
  !$omp atomic update
  wait_time = wait_time + used
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
    flush (unit)
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
  integer:: level
  real(8):: wc, used
  !-----------------------------------------------------------------------------
  if (.not.omp_lock%on) return
#ifdef TRACE 
  call trace%begin('lock_t%unset')
#endif
  !-----------------------------------------------------------------------------
  ! It should be impossible to come here with an uninitialized lock
  !-----------------------------------------------------------------------------
  if (.not. self%initialized) then
    write (stderr,'(a,i4)') &
      'omp_lock_t%unset ERROR: lock not initialized, thread =', omp%thread
    stop 'omp_lock_t%unset'
  end if
  self%level = self%level-1
  level = self%level
  !$omp flush
  call omp_unset_nest_lock (self%lock)
  !$omp flush
  wc = wallclock()
  used = wc - self%start_time
  !call alternate_log (self)
  if (wc < verbose_time .or. omp_lock%verbose > 1) then
    if (present(label)) then
      write (unit,'(f10.6,a,a4,i6,i3,i4,6p,f10.1,2x,a)') wc, &
        ' lock kind, id, level, thread = ', self%kind, self%id, level, &
        omp%thread, used, ' unset at '//label
    else
      write (unit,'(f10.6,a,a4,i6,i3,i4,6p,f10.1,2x,a)') wc, &
        ' lock kind, id, level, thread = ', self%kind, self%id, level, &
        omp%thread, used, ' unset'
    end if
    flush (unit)
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
SUBROUTINE info (self, unit, p_verbose)
  class(lock_t):: self
  integer:: unit, p_verbose
  class(lock_t), pointer:: lock
  real(8):: time
  !-----------------------------------------------------------------------------
  if (.not.omp_lock%on) return
  if (omp_lock%verbose > -1) then
    !$omp atomic read
    time = wait_time
    write (stdout,'(a,1p,g11.2,a,10f8.1)') ' total lock waiting time =', &
      time, ' core-s'
    !$omp atomic write
    wait_time = 0d0
  end if
  return
END SUBROUTINE info

!===============================================================================
END MODULE omp_lock_mod
