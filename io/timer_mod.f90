!===============================================================================
!> Each thread uses a private timer data type, with arrays for start time and
!> total time for each registered procedure.  At begin, it stores the start time,
!> and at end it increments the call counter and sums the time used into the
!> total counter.
!>
!> We suspend counting of the active procedure if a new begin occurs, and resume
!> it again on the next end call.  Each subsequent begin call suspends the
!> ongoing call, so we need a list of active counters.
!===============================================================================
MODULE timer_mod
  USE mpi_mod
  USE io_unit_mod
  USE omp_mod
  USE omp_timer_mod
  USE omp_lock_mod
  USE omp_lib
  implicit none
  private
  integer, parameter:: maxproc=100, maxdepth=20
  integer, save:: ntimer=0
  type:: latency_t
    real:: max=0.0
    real(8):: aver=0d0
    integer:: n=0
  end type
  type:: timer_t
    integer(8):: calls(maxproc)=0
    real(8), dimension(maxproc):: start=0d0, total=0d0
    integer, dimension(maxdepth):: stack
    integer:: istack=0
    real(8):: bytes_recv=0.0_8
    real(8):: sec_per_report=10d0
    integer(8):: n_master(4)=0
    integer(8):: n_update=0
    integer(8):: n_recv=0_8
    integer(8):: mpi_test=0_8, mpi_hit=0_8
    integer(8):: mem_test=0_8, mem_hit=0_8
    real(8):: busy_time=0d0, spin_time=0d0
    integer:: nq_send_max=0
    integer:: nq_recv_max=0
    integer:: n_lines=0
    integer:: n_mhd=0, n_solve=0
    type(latency_t):: latency
    real(8), allocatable:: levelcost(:)
    integer(8), allocatable:: levelpop(:)
    integer:: levelmin, levelmax
    real:: dead_mans_hand=0.0
    logical:: detailed=.false.
  contains
    procedure:: init
    procedure:: begin
    procedure:: end
    procedure:: print
    procedure, nopass:: tic
    procedure, nopass:: toc_1
    procedure, nopass:: toc_i4
    procedure, nopass:: toc_i8
    procedure, nopass:: toc_r
    generic:: toc => toc_1, toc_i4, toc_i8, toc_r
  end type
  character(len=64), dimension(maxproc):: names
  type(timer_t), pointer:: timers(:)
  type(timer_t), public:: timer
  integer, save:: nprint=20
  integer, save:: verbose=0
  real(8), save:: next_report_sec=0d0
  real(8), save:: prev2=0d0
  real, save:: min_fraction=0.0005
  interface toc
    module procedure toc_i4, toc_i8, toc_1
  end interface
  real(8), save:: wt=0d0
  !$omp threadprivate (wt)
  public tic, toc
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self)
  class(timer_t):: self
  logical, save:: first_time=.true., detailed=.false.
  integer:: ignore, i
  real(8):: sec_per_report=10.
  namelist /timer_params/ sec_per_report, min_fraction, verbose, detailed
  !.............................................................................
  !$omp critical (timer_cr)
  if (first_time) then
    first_time = .false.
    rewind (io_unit%input)
    read (io_unit%input, timer_params, iostat=ignore)
    if (io_unit%master) write (*, timer_params)
    self%sec_per_report = sec_per_report
    next_report_sec = wallclock() + sec_per_report
    omp%nthreads = max(1,omp%nthreads)
    timer%detailed = detailed
    allocate (timers(0:omp%nthreads-1))
  end if
  !$omp end critical (timer_cr)
  ntimer = 0
  do i=0,omp%nthreads-1
    timers(i)%calls = 0
    timers(i)%start = 0d0
    timers(i)%total = 0d0
  end do
END SUBROUTINE init

!===============================================================================
!> If itimer==0, find a suitable value for it, and store the name with it.
!> The name argument is not needed when itimer has already been assigned a value.
!===============================================================================
SUBROUTINE begin (self, name, itimer)
  class(timer_t):: self
  character(len=*):: name
  integer:: itimer
  !.............................................................................
  integer:: thread, i, otimer, istack
  logical:: ok
  real(8):: wc
  !-----------------------------------------------------------------------------
  if (io_unit%do_validate) return
  thread = omp_get_thread_num()
  if (itimer==0) then
    !$omp critical (timer_cr)
    if (itimer==0) then
      ok = .true.
      do i=1,ntimer
        ok = ok .and. trim(names(i)) /= trim(name)
      end do
      if (ok) then
        ntimer = ntimer+1
        if (ntimer > maxproc) call mpi%abort ('timer_mod: increase maxproc')
        itimer = ntimer
        names(itimer) = name
      else
        write (stdout,*) 'timer_t%begin: WARNING, duplicate name: '//trim(name)
      end if
    end if
    !$omp end critical (timer_cr)
    !---------------------------------------------------------------------------
    ! If the name already exists in names then itimer remains =0
    !---------------------------------------------------------------------------
    if (itimer <= 0) return
  end if
  !-----------------------------------------------------------------------------
  ! If a procedure is being timed, accumulate to its total
  !-----------------------------------------------------------------------------
  istack = timers(thread)%istack
  wc = wallclock()
  if (istack>0) then
    otimer = timers(thread)%stack(istack)
    timers(thread)%total(otimer) = timers(thread)%total(otimer) &
                                 + (wc - timers(thread)%start(otimer))
  end if
  !-----------------------------------------------------------------------------
  ! Start the new timer, and record it on the stack
  !-----------------------------------------------------------------------------
  timers(thread)%start(itimer) = wc
  istack = istack+1
  if (istack>maxdepth) call mpi%abort ('too many timer levels;'//trim(name))
  timers(thread)%stack(istack) = itimer
  timers(thread)%istack = istack
  if (verbose > 0) then
    print '(a,2i5,3x,a)', 'timer%begin: istack, itimer, name =', istack, itimer, trim(name)
  end if
END SUBROUTINE begin

!===============================================================================
!> At end call, increment number of calls and total time for the active timer
!> and resume any suspended timer
!===============================================================================
SUBROUTINE end (self, itimer)
  class(timer_t):: self
  integer, optional:: itimer
  integer::thread, istack, otimer
  real(8):: wc
  !.............................................................................
  if (io_unit%do_validate) return
  !-----------------------------------------------------------------------------
  ! Find the active timer for this thread
  !-----------------------------------------------------------------------------
  thread = omp_get_thread_num()
  istack = timers(thread)%istack
  if (istack <= 0) then
    print *, 'WARNING, timer%end out-of-sync: thread, istack =', omp%thread, istack
    timers(thread)%istack = 0
    return
  end if
  otimer = timers(thread)%stack(istack)
  if (verbose > 0) then
    print '(a,2i5,3x,a)', 'timer%end:   istack, otimer, name =', istack, otimer, trim(names(otimer))
  end if
  !-----------------------------------------------------------------------------
  ! Increment call counter and total time used
  !-----------------------------------------------------------------------------
  wc = wallclock()
  if (present(itimer)) then
    if (itimer /= otimer) then
      write (stdout,'(i6,2(3x,a,i4))') &
        omp%thread, 'WARNING: itimer =',itimer, & 
        ' not equal to expected value =', otimer
    end if
  end if
  timers(thread)%calls(otimer) = timers(thread)%calls(otimer) + 1
  timers(thread)%total(otimer) = timers(thread)%total(otimer) + &
                                 wc - timers(thread)%start(otimer)
  !-----------------------------------------------------------------------------
  ! Check if there is a suspended timer present, and if so restart it
  !-----------------------------------------------------------------------------
  istack = istack-1
  if (istack>0) then
    otimer = timers(thread)%stack(istack)
    timers(thread)%start(otimer) = wc
  else if (istack<0) then
    print *,'WARNNING: too many timer%end calls'
    istack = 0
  end if
  timers(thread)%istack = istack
  !-----------------------------------------------------------------------------
  ! Check if it is time to do a printout
  !-----------------------------------------------------------------------------
  if (wc > next_report_sec) then
    !$omp critical (timer_print_cr)
    if (wc > next_report_sec) then
      call real_print (self)
      next_report_sec = (nint(wc/self%sec_per_report)+1)*self%sec_per_report
    end if
    !$omp end critical (timer_print_cr)
  end if
END SUBROUTINE end

!===============================================================================
!> Printout times used in routines using more than 0.05% of the time, and reset
!> counters.
!===============================================================================
SUBROUTINE print (self)
  class(timer_t):: self
  if (omp%master) call real_print (self)
END SUBROUTINE print
SUBROUTINE real_print (self)
  class(timer_t):: self
  real(8):: total, proc(maxproc), ncalls(maxproc)
  real(8):: updates, musppt, calls
  real(8), save:: n_previous=0
  real(8):: wc, sec
  integer:: i, n_recv, lntimer
  !.............................................................................
  if (io_unit%do_validate) return
  !-----------------------------------------------------------------------------
  ! Prevent several treads from printing at the same time
  !-----------------------------------------------------------------------------
  wc = wallclock()
  sec = wc - prev2
  total = 0d0
  calls = 0d0
  !$omp atomic read
  lntimer = ntimer
  do i=1,lntimer
    proc(i) = sum(timers(:)%total(i))
    ncalls(i) = sum(timers(:)%calls(i))
    total = total + proc(i)
    calls = calls + ncalls(i)
  end do
  !-----------------------------------------------------------------------------
  ! 'sec' measure wall clock time on a single thread, while 'total' is a sum over
  ! all procedures being measured, over all threads.
  ! 'updates' is the total number of active-cell updates (i.e. not including
  ! ghost cells) across *all* threads.
  ! `musppt` is defined to be the number of microseconds per update per thread.
  !-----------------------------------------------------------------------------
  updates = max(1.0,real(timer%n_update-n_previous))
  musppt = sec*1d6*min(omp%nthreads,omp%ncores)/updates
  n_recv = max(timer%n_recv,1)
  if (mpi%master) &
    call write_to (io_unit%output)
  call write_to (io_unit%mpi)
  !$omp atomic write
  timer%n_lines = 1
  timer%bytes_recv = 0.0_8
  timer%n_recv = 0_8
  timer%nq_send_max = 0
  timer%nq_recv_max = 0
  n_previous = timer%n_update
  prev2 = wc
  timer%latency%max  = 0.0
  timer%latency%aver = 0d0
  timer%latency%n    = 0
  timer%mpi_test = 0_8
  timer%mpi_hit = 0_8
  do i=1,ntimer
    timers(:)%total(i) = 0d0
    timers(:)%calls(i) = 0
  end do
  return
contains
  subroutine write_to (unit)
    integer:: unit, level, popu
    real:: cost, per_call
    logical:: warn, mesg
    write (unit,'(23x,a,7x,a,10x,a,7x,a,5x,a)') 'procedure', 'calls', 'time', 'percent', ' s/call'
    if (unit==io_unit%output) then
      !$omp atomic write
      timer%n_lines = 1
    end if
    if (total > 0.0) then
      mesg = .false.
      do i=1,ntimer
        if (proc(i)/total > min_fraction) then
          per_call = proc(i)/max(ncalls(i),1d0)
          warn = per_call < 0.5e-6
          mesg = mesg .or. warn
          write (unit,'(a32,1p,2g15.3,0p,f10.1,0p,f12.6,1x,a1)') &
            trim(names(i)), ncalls(i), proc(i), proc(i)/total*100., &
            per_call, merge('W',' ',warn)
        end if
      end do
      if (mesg) &
        write (unit,*) 'W: much of the item time may be due to timer calls'
    end if
    !---------------------------------------------------------------------------
    ! Output AMR level costs
    !---------------------------------------------------------------------------
    if (allocated(timer%levelcost) .and. timer%levelmax > timer%levelmin) then
      do level=timer%levelmin,timer%levelmax
        !$omp atomic read
        cost = timer%levelcost(level)
        !$omp atomic read
        popu = timer%levelpop(level)
        write(unit,'(22x,a,i3,g12.3,i7)') &
          'level, cost, population =', level, cost, popu
      end do
      do level=timer%levelmin,timer%levelmax
        !$omp atomic write
        timer%levelcost(level) = 0d0
      end do
    end if
    !---------------------------------------------------------------------------
    ! MPI and OMP statistics
    !---------------------------------------------------------------------------
    write (unit,'(a32,1p,2g15.3,g14.4,2(g12.3,2x,a,5x))') &
      "TOTAL thread time, calls", calls, total, 100., musppt, 'core-mus/cell-upd', wc, 'wall sec'
    if (mpi%size>1) then
      timer%latency%aver = timer%latency%aver/max(1,timer%latency%n)
      write (unit,'(1x,"MPI recv:",f10.1," MB/s",f11.3," MB/mesg",i8," nq_send_max",i8, &
        " nq_recv_max",2f8.3," max, aver latency",f7.2, " f_unpk",f7.2," f_mem",f8.3," f_q")') &
        timer%bytes_recv/1024.**2/sec, timer%bytes_recv/1024.**2/n_recv, &
        timer%nq_send_max, timer%nq_recv_max, timer%latency%max, timer%latency%aver, &
        timer%mpi_hit/real(max(1,timer%mpi_test)), &
        timer%mem_hit/real(max(1,timer%mem_test)), &
        timer%busy_time/max(timer%spin_time + timer%busy_time,1d-30)
      timer%mpi_hit = 0
      timer%mpi_test = 0
      timer%busy_time = 0d0
      timer%spin_time = 0d0
    end if
    call omp_lock%info (unit)
    flush (unit)
    self%n_mhd = 0; self%n_solve=0
  end subroutine write_to
END SUBROUTINE real_print

!===============================================================================
!===============================================================================
SUBROUTINE tic (time)
  implicit none
  real(8), optional:: time
  !.............................................................................
  if (present(time)) then
    time = wallclock()
  else
    wt = wallclock()
  end if
END SUBROUTINE tic

!===============================================================================
!===============================================================================
SUBROUTINE toc_1 (label, time)
  implicit none
  character(len=*):: label
  real(8), optional:: time
  integer:: n
  !.............................................................................
  call toc_r (label, 1.0, time)
END SUBROUTINE toc_1

!===============================================================================
!===============================================================================
SUBROUTINE toc_i4 (label, n, time)
  implicit none
  character(len=*):: label
  real(8), optional:: time
  integer:: n
  !.............................................................................
  call toc_r (label, real(n), time)
END SUBROUTINE toc_i4

!===============================================================================
!===============================================================================
SUBROUTINE toc_i8 (label, n, time)
  implicit none
  character(len=*):: label
  real(8), optional:: time
  integer(8):: n
  !.............................................................................
  call toc_r (label, real(n), time)
END SUBROUTINE toc_i8

!===============================================================================
FUNCTION str_f4 (t) RESULT (out)
  real(8):: t
  character(len=8):: out
  !-----------------------------------------------------------------------------
  write (out,'(f8.3)') t
  out = adjustl(out)
  out = trim(out(1:4))
END FUNCTION

!===============================================================================
FUNCTION time_str (t) RESULT (out)
  real(8):: t
  character(len=8):: out
  real, parameter:: ns=1e-9, mu=1e-6, ms=1e-3, s=1e0, mi=6e1, hr=36e2, dy=24*hr
  real, parameter:: wk=7*dy, yr=365.*dy, mo=yr/12.
  !-----------------------------------------------------------------------------
  if      (t/mu < 1d0) then
    out=trim(str_f4(t/ns))//' ns'
    1 format(f7.3,1x,a2)
  else if (t/ms < 1d0) then
    out=trim(str_f4(t/mu))//' mus'
  else if (   t < 1d0) then
    out=trim(str_f4(t/ms))//' ms'
  else if (   t < 100) then
    out=trim(str_f4(   t))//' s '
  else if (t/mi < 100) then
    out=trim(str_f4(t/mi))//' mn'
  else if (t/dy < 1d0) then
    out=trim(str_f4(t/hr))//' hr'
  else if (t/wk < 1d0) then
    out=trim(str_f4(t/dy))//' dy'
  else if (t/mo < 1d0) then
    out=trim(str_f4(t/wk))//' wk'
  else if (t/yr < 1d0) then
    out=trim(str_f4(t/mo))//' mo'
  else
    out=trim(str_f4(t/yr))//' yr'
  end if
END FUNCTION

!===============================================================================
!> Print the time used per thread.  Overuse is compensated for by the factor f =
!> the number of threads per core.
!===============================================================================
SUBROUTINE toc_r (label, n, time)
  implicit none
  character(len=*):: label
  real:: n
  real(8), optional:: time
  real(8):: dtime, now
  integer:: n_cores, n_threads
  !.............................................................................
  if (io_unit%do_validate) return
  now = wallclock()
  if (present(time)) then
    dtime = (now - time)
    time = now
  else
    dtime = (now - wt)
    wt = now
  end if
  n_cores = max(omp%ncores,1)
  n_threads = max(omp%nthreads,1)
  if (n_cores > n_threads) n_cores = min(n_cores,n_threads)

  1 format (1x, a, ':  ', a, ',  ', a, a, 1p, e10.2, 1x, a, 3(i4, 1x, a))
  2 format (1x, a, ':  ', a, ',  ', a, a, 1p, e10.2, 1x, a, 4(i4, 1x, a))
  if (n_cores == n_threads) then
    write (io_unit%log,1) label, trim(time_str(dtime)), &
      trim(time_str(dtime*n_cores/max(n,1.))), '/upd,', n, 'updates/process,', &
      mpi%size, merge('processes,','process,  ',mpi%size>1), &
      n_cores, merge('cores/process','core/process ',n_cores>1)
  else
    write (io_unit%log,1) label, trim(time_str(dtime)), &
      trim(time_str(dtime*n_cores/max(n,1.))), '/upd,', n, 'updates/process,', &
      mpi%size, merge('processes,','process,  ',mpi%size>1), &
      n_cores, merge('cores/process','core/process ',n_cores>1), &
      n_threads, merge('threads/process','thread/process ',n_threads>1)
  end if
  flush (io_unit%log)
  if (mpi%rank==0.and.omp%thread==0) then
    if (n>1.0) then
      if (n_cores == n_threads) then
        print 1, label, trim(time_str(dtime)), &
          trim(time_str(dtime*n_cores/max(n,1.))), '/upd,', n, 'updates/process,', &
          mpi%size, merge('processes,','process,  ',mpi%size>1), &
          n_cores, merge('cores/process','core/process ',n_cores>1)
      else
        print 1, label, trim(time_str(dtime)), &
          trim(time_str(dtime*n_cores/max(n,1.))), '/upd,', n, 'updates/process,', &
          mpi%size, merge('processes,','process,  ',mpi%size>1), &
          n_cores, merge('cores/process','core/process ',n_cores>1), &
          n_threads, merge('threads/process','thread/process ',n_threads>1)
      end if
    else
      print 1, label, trim(time_str(dtime))
    end if
    flush (6)
  end if
  if (present(time)) then
    time = wallclock()
  else
    wt = wallclock()
  end if
END SUBROUTINE toc_r

END MODULE timer_mod
