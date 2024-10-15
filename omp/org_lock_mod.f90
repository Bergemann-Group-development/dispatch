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
MODULE org_lock_mod
    USE io_unit_mod
    USE omp_mod
    USE omp_timer_mod
    USE omp_lib
    USE shared_mod
!define TRACE
#ifdef TRACE 
    USE trace_mod
#endif
    implicit none
    private
    integer, parameter:: mkinds=20, maxproc=100
    type, public:: org_lock_t
      integer:: id=-1
      integer:: verbose=0, unit=0
      integer(kind=omp_nest_lock_kind):: lock
      integer:: thread=-1
      character(len=5):: kind = 'orig'
      integer:: level=0, readers=0
      logical:: on = .true.
      logical:: log = .false.
      logical:: initialized = .false.
      real(8):: time(maxproc) = 0d0                 ! lock waiting times
      real(8):: start_time = 0d0
      real(8):: wait_time = 0d0
      logical:: links = .false.             ! lock nbor lists
      logical:: tasks = .true.              ! lock task memory updates
      logical:: fall_through = .true.
    contains
      procedure:: init
      procedure:: init_unit
      procedure:: set
      procedure:: test
      procedure:: check
      procedure:: unset
      procedure:: destroy
      procedure:: info
      procedure:: unset_log
    end type
    real, save:: verbose_time=0.0
    integer, save:: id_save=-1, n_lock=0
    integer, save:: verbose=0
    integer, save:: unit=0
    integer, save:: single_thread=-1
    !$omp threadprivate (unit)
    logical, save:: fall_through=.true.
    character(len=6):: log='lock'
    character(len=6):: ext, ext1, ext2
    real(8):: log_interval=0d0, alternate_next=0d0
    real(8):: wait_list = 0d0
    real(8):: wait_link = 0d0
    real(8):: wait_mem = 0d0
    real(8):: log_time = 0d0
  CONTAINS
  
  !===============================================================================
  !> Initialize a lock, taking care not to do it more than once
  !===============================================================================
  SUBROUTINE init (self, kind, id)
    class(org_lock_t), target:: self
    character(len=*), optional:: kind
    integer, optional:: id
    !.............................................................................
    logical omp_in_parallel
    logical, save:: tasks=.true., links=.false.
    integer:: iostat, i, u
    logical:: set_id
    logical, save:: on = .true.
    logical, save:: first_time=.true.
    namelist /lock_params/ on, verbose, verbose_time, log, tasks, links, &
      log_interval, fall_through, single_thread
    character(len=120):: ids= &
    '$Id$'
    !-----------------------------------------------------------------------------
    call io_unit%print_id (ids)
    !-----------------------------------------------------------------------------
#ifdef TRACE 
    call trace%begin('org_lock_t%init')
#endif
    !-----------------------------------------------------------------------------
    ! This part is executed only once -- called from critical region
    !-----------------------------------------------------------------------------
    if (first_time) then
      !---------------------------------------------------------------------------
      ! By default don't use locks when running single-threaded, but allow
      ! override, for recording lock sequence etc.
      !---------------------------------------------------------------------------
      if (omp%nthreads == 1) then
        on = .false.
      end if
      rewind (io_unit%input)
      read (io_unit%input, lock_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io_unit%input)
        read (io_unit%input, lock_params)
      end if
      write (io_unit%nml, lock_params)
      first_time = .false.
    end if
    self%on = on
    self%verbose = verbose
    self%tasks = tasks
    self%links = links
    self%fall_through = fall_through
    !-----------------------------------------------------------------------------
    ! This part is executed once per lock, initializing it
    !-----------------------------------------------------------------------------
    if (.not.self%initialized) then
      !$omp critical (lock_cr)
      if (.not.self%initialized) then
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
          self%kind = kind(1:min(len(kind),5))
        else
          self%kind = 'void'
        end if
        self%on  = on
        n_lock = n_lock+1
        self%initialized = .true.
      end if
      !$omp end critical (lock_cr)
    end if
    self%unit = unit
#ifdef TRACE 
    call trace%end()
#endif
  END SUBROUTINE init
  
  !===============================================================================
  !> Choose output files and units -- this must be called AFTER log files have
  !> been opened!
  !===============================================================================
  SUBROUTINE init_unit (self)
    class(org_lock_t):: self
    character(len=256):: filename
    !-----------------------------------------------------------------------------
    ! This part is executed once per thread, choosing unit, depending on params
    !-----------------------------------------------------------------------------
    if (unit == 0) then
      !$omp critical (unit_cr)
      if (unit == 0) then
        if (log_interval > 0d0) then
          filename = trim(io_unit%rankbase)//'.lock1'
          open (io_unit%lock1, file=trim(filename), form='formatted', &
                status='unknown')
          filename = trim(io_unit%rankbase)//'.lock2'
          open (io_unit%lock2, file=trim(filename), form='formatted', &
                status='unknown')
          alternate_next = wallclock() + log_interval
        else
          filename = trim(io_unit%rankbase)//'.lock'
          open (io_unit%lock1, file=trim(filename), form='formatted', &
                status='unknown')
        end if
        unit = io_unit%lock1
      end if
      !$omp end critical (unit_cr)
    end if
  END SUBROUTINE init_unit
  
  !===============================================================================
  !> Local clock function, returning the fine-grained OMP time when verbose > 0
  !===============================================================================
  FUNCTION wclock()
    real(8):: wclock
    !.............................................................................
    if (verbose > 0) then
      wclock = wallclock()
    else
      wclock = 0d0
    end if
  END FUNCTION
  
  !===============================================================================
  !> Return with the lock, waiting if necessary, and creating the lock if needed
  !===============================================================================
  SUBROUTINE set (self, label, itimer)
    class(org_lock_t):: self
    character(len=*), optional:: label
    integer, optional:: itimer
    character(len=32):: filename
    real(8):: start, end, used
    logical:: self_log
    integer:: thread
    !-----------------------------------------------------------------------------
    if (.not. self%on) return
#ifdef TRACE 
    call trace%begin('org_lock_t%set')
#endif
    if (.not. self%initialized) then
      write (stderr,'(a,i4)') &
      'org_lock_t%set WARNING: lock not initialized, thread =', omp%thread
      write (io_unit%log,'(f14.6,2x,a,i4)') wallclock(), &
      'org_lock_t%set WARNING: lock not initialized, thread =', omp%thread
      call self%init
    end if
    !-----------------------------------------------------------------------------
    ! Optionally, redirect lock-logging, to avoid slow-down
    !-----------------------------------------------------------------------------
    if (log == 'stdlog') then
      if (unit /= io_unit%log) then
        !$omp atomic write
        unit = io_unit%log
      end if
    end if
    !-----------------------------------------------------------------------------
    ! Must print log info before the lock, to reveal dead-locks
    !-----------------------------------------------------------------------------
    call alternate_log (self)
    start = wclock()
    self_log = (start < verbose_time) .or. (verbose > 1)
    !$omp atomic write
    self%log = self_log
    if ((self%log .and. verbose > 2) .or. omp%thread == single_thread) then
      if (present(label)) then
        write (unit,'(f14.6,a,a5,i10,2i3,i4,12x,a)') start, &
          ' lock: kind,id,level,thread= ', self%kind, self%id, &
          self%level, self%readers, omp_get_thread_num(), '   get at '//label
      else
        write (unit,'(f14.6,a,a5,i10,2i3,i4,12x,a)') start, &
          ' lock: kind,id,level,thread= ', self%kind, self%id, &
          self%level, self%readers, omp_get_thread_num(), '   get'
      end if
      !flush (unit)
      used = wclock()-start
      !$omp atomic update
      log_time = log_time + used
    end if
    !$omp flush
    call omp_set_nest_lock (self%lock)
    !$omp flush
    end = wclock()
    used = end-start
    self%start_time = end
    !$omp atomic update
    self%wait_time = self%wait_time + used
    if (self%kind == 'list') then
      !$omp atomic update
      wait_list = wait_list + used
    else if (self%kind == 'link') then
      !$omp atomic update
      wait_link = wait_link + used
    else if (self%kind == 'meml') then
      !$omp atomic update
      wait_mem = wait_mem + used
    end if
    self%level = self%level+1
    !$omp flush
    if (self%log .or. (omp%thread == single_thread)) then
      if (present(label)) then
        write (unit,'(f14.6,a,a5,i10,2i3,i4,6p,f10.1,2x,a)') end, &
          ' lock: kind,id,level,thread= ', self%kind, self%id, &
          self%level, self%readers, omp%thread, used, '   set at '//label
      else
        write (unit,'(f14.6,a,a5,i10,2i3,i4,6p,f10.1,2x,a)') end, &
          ' lock: kind,id,level,thread= ', self%kind, self%id, &
          self%level, self%readers, omp%thread, used, '   set'
      end if
      !flush (unit)
      used = wclock()-start
      !$omp atomic update
      log_time = log_time + used
    end if
#ifdef TRACE 
    call trace%end()
#endif
  END SUBROUTINE set
  
  !===============================================================================
  !> Alternate between two log files
  !===============================================================================
  SUBROUTINE alternate_log (self)
    class(org_lock_t):: self
    real(8):: wc
    !-----------------------------------------------------------------------------
    if (.not. self%on) return
    if (log_interval > 0d0) then
      wc = wallclock()
      if (wc > alternate_next) then
        write (io_unit%log,*) omp%thread,'waiting at lock_cr'
        if (io_unit%flush_logs) &
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
          alternate_next = wc + log_interval
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
    class(org_lock_t):: self
    logical:: out
    real(8):: wc
    !-----------------------------------------------------------------------------
    if (.not. self%on) return
#ifdef TRACE 
    call trace%begin('org_lock_t%test')
#endif
    if (omp%thread == self%thread) then
      out = .true.
      !print *, omp%thread, ' already has lock id', self%id
    else
      !print *, omp%thread, ' omp_test_lock for id', self%id
      out =  ( omp_test_nest_lock (self%lock) .gt. 0 )
      if (out) then
        self%thread = omp%thread
        !$omp atomic update
        self%level = self%level+1
        wc = wallclock()
        if (wc < verbose_time .or. verbose > 1 .or. log_interval > 0d0) then
          write (io_unit%log,'(f14.6,i4,a,i10,i3,a)') wc, omp%thread, ' lock id', &
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
    class(org_lock_t):: self
    logical:: out
    !-----------------------------------------------------------------------------
    if (.not. self%on) return
    out = (omp%thread == self%thread)
  END FUNCTION check
  
  !===============================================================================
  !> Unset the associated thread number and release the lock (in that order!)
  !===============================================================================
  SUBROUTINE unset (self, label)
    class(org_lock_t):: self
    character(len=*), optional:: label
    integer:: level
    real(8):: wc, used
    !-----------------------------------------------------------------------------
    if (.not. self%on) return
#ifdef TRACE 
    call trace%begin('org_lock_t%unset')
#endif
    !-----------------------------------------------------------------------------
    ! It should be impossible to come here with an uninitialized lock
    !-----------------------------------------------------------------------------
    if (.not. self%initialized) then
      write (stderr,'(a,i4)') &
        'org_lock_t%unset ERROR: lock not initialized, thread =', omp%thread
      stop 'org_lock_t%unset'
    end if
    self%level = self%level-1
    level = self%level
    !$omp flush
    call omp_unset_nest_lock (self%lock)
    !$omp flush
    wc = wclock()
    used = wc - self%start_time
    !call alternate_log (self)
    if ((self%log .and. verbose > 2) .or. omp%thread == single_thread) then
      if (present(label)) then
        write (unit,'(f14.6,a,a5,i10,2i3,i4,6p,f10.1,2x,a)') wc, &
          ' lock: kind,id,level,thread= ', self%kind, self%id, level, &
          self%readers, omp%thread, used, ' unset at '//label
      else
        write (unit,'(f14.6,a,a5,i10,2i3,i4,6p,f10.1,2x,a)') wc, &
          ' lock: kind,id,level,thread= ', self%kind, self%id, level, &
          self%readers, omp%thread, used, ' unset'
      end if
      if (io_unit%flush_logs) &
        flush (unit)
    end if
#ifdef TRACE 
    call trace%end()
#endif
  END SUBROUTINE unset
  
  !===============================================================================
  SUBROUTINE destroy (self)
    class(org_lock_t):: self
    !-----------------------------------------------------------------------------
    if (self%initialized) then
      call omp_destroy_nest_lock (self%lock)
      self%initialized = .false.
    end if
  END SUBROUTINE destroy
  
  !===============================================================================
  SUBROUTINE info (self, unit, p_verbose)
    class(org_lock_t):: self
    integer:: unit, p_verbose, i
    class(org_lock_t), pointer:: lock
    real(8):: time, times(3), total
    real(8), save:: prev=0d0, prevs(3)=0d0
    !-----------------------------------------------------------------------------
    if (.not. self%on) return
    if (verbose > 0) then
      !$omp atomic read
      time = self%wait_time
      !$omp atomic read
      times(1) = wait_list
      !$omp atomic read
      times(2) = wait_link
      !$omp atomic read
      times(3) = wait_mem
      times = times-prevs
      total = sum(times)
      write (stdout,'(4(a,f6.1,3x),a,2p,f6.1,1x,a)') &
        ' lock waiting time, lists:', times(1), &
        ' links:', times(2), &
        ' mem:', times(3), &
        ' total:', total, &
        'core-s =', total/(shared%sec_per_report*omp%nthreads), '%'
      !$omp atomic write
      prev = time
      do i=1,3
        !$omp atomic write
        prevs(i) = times(i)
      end do
      !---------------------------------------------------------------------------
      !$omp atomic read
      time = log_time
      write (stdout,'(a,f6.2,a)') ' lock     log time =', &
        time, ' core-s'
      !$omp atomic write
      log_time = 0d0
    end if
    return
  END SUBROUTINE info
  
  !===============================================================================
  SUBROUTINE unset_log (self, label)
    class(org_lock_t):: self
    character(len=*), optional:: label
    real(8):: wc, used
    !-----------------------------------------------------------------------------
    if (self%log .or. (omp%thread == single_thread)) then
      wc = wclock()
      used = wc - self%start_time
      if (present(label)) then
        write (unit,'(f14.6,a,a5,i10,2i3,i4,6p,f10.1,2x,a)') wc, &
          ' lock: kind,id,level,thread= ', self%kind, self%id, self%level, &
          self%readers, omp%thread, used, 'unread at R'//label
      else
        write (unit,'(f14.6,a,a5,i10,2i3,i4,6p,f10.1,2x,a)') wc, &
          ' lock: kind,id,level,thread= ', self%kind, self%id, self%level, &
          self%readers, omp%thread, used, 'unread'
      end if
      if (io_unit%flush_logs) &
        flush (unit)
    end if
  END SUBROUTINE unset_log
  
  !===============================================================================
  END MODULE org_lock_mod
  