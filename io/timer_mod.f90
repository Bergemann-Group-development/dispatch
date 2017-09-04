!> =============================================================================
!> Each thread uses a private timer data type, with arrays for start time and
!> total time for each registered procedure.  At begin, it stores the start time,
!> and at end it increments the call counter and sums the time used into the
!> total counter.
!>
!> We suspend counting of the active procedure if a new begin occurs, and resume
!> it again on the next end call.  Each subsequent begin call suspends the
!> ongoing call, so we need a list of active counters.
!> =============================================================================
MODULE timer_mod
  USE mpi_mod
  USE io_unit_mod
  USE omp_mod
  USE omp_timer_mod
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
    integer(8):: n_recv=0_8, n_test=0_8
    integer:: nq_send_max=0
    integer:: nq_recv_max=0
    real:: gb=0.0
    type(latency_t):: latency
  contains
    procedure:: init
    procedure:: begin
    procedure:: end
    procedure, nopass:: print
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
  real(8), save:: sec=0d0, next_report_sec=0d0
  real(8), save:: prev2=0d0
  integer(8), save:: n_update=0
  interface toc
    module procedure toc_i4, toc_i8, toc_1
  end interface
  real(8), save:: wt=0d0
  !$omp threadprivate (wt)
  public tic, toc, n_update
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self)
  class(timer_t):: self
  logical, save:: first_time=.true.
  integer:: ignore, i
  real(8):: sec_per_report=10.
  namelist /timer_params/ sec_per_report
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

!> =============================================================================
!> If itimer==0, find a suitable value for it, and store the name with it.
!> The name argument is not needed when itimer has already been assigned a value.
!> =============================================================================
SUBROUTINE begin (self, name, itimer)
  class(timer_t):: self
  character(len=*), optional:: name
  integer:: itimer, thread, i, otimer, istack
  logical:: ok
  real(8):: wc
  !.............................................................................
  if (io_unit%do_validate) return
  thread = omp_get_thread_num()
  if (itimer==0) then
    !$omp critical (timer_cr)
    ok = .true.
    do i=1,ntimer
      ok = ok .and. trim(names(i)) /= trim(name)
    end do
    if (ok) then
      ntimer = ntimer+1
      if (ntimer > maxproc) call mpi%abort ('timer_mod: increase maxproc')
      itimer = ntimer
      names(itimer) = name
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
  if (istack>0) then
    otimer = timers(thread)%stack(istack)
    !print *,'suspending timer for '//trim(names(otimer))
    wc = wallclock()
    timers(thread)%total(otimer) = timers(thread)%total(otimer) &
                                 + (wc - timers(thread)%start(otimer))
  end if
  !-----------------------------------------------------------------------------
  ! Start the new timer, and record it on the stack
  !-----------------------------------------------------------------------------
  timers(thread)%start(itimer) = wallclock()
  istack = istack+1
  if (istack>maxdepth) call mpi%abort ('too many timer levels')
  timers(thread)%stack(istack) = itimer
  timers(thread)%istack = istack
END SUBROUTINE begin

!> =============================================================================
!> At end call, increment number of calls and total time for the active timer
!> and resume any suspended timer
!> =============================================================================
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
  otimer = timers(thread)%stack(istack)
  !-----------------------------------------------------------------------------
  ! Increment call counter and total time used
  !-----------------------------------------------------------------------------
  wc = wallclock()
  if (present(itimer)) then
    timers(thread)%calls(otimer) = timers(thread)%calls(otimer) + 1
    if (itimer == 0) then
      print *, 'WARNING: itimer should not be zero!'
    else
      timers(thread)%total(otimer) = timers(thread)%total(otimer) + &
                                     wc - timers(thread)%start(itimer)
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Check if there is a suspended timer present, and if so restart it
  !-----------------------------------------------------------------------------
  istack = istack-1
  if (istack>0) then
    otimer = timers(thread)%stack(istack)
    timers(thread)%start(otimer) = wc
    !print *,'restarting timer for '//trim(names(otimer))
  else if (istack<0) then
    print *,'WARNNING: too many timer%end calls'
    istack = 0
  end if
  timers(thread)%istack = istack
  !-----------------------------------------------------------------------------
  ! Check if it is time to do a printout
  !-----------------------------------------------------------------------------
  if (io_unit%master) then
    if (wc > next_report_sec) then
      call print()
      next_report_sec = wc + self%sec_per_report
    end if
  end if
END SUBROUTINE end

!> =============================================================================
!> =============================================================================
SUBROUTINE print ()
  real(8):: total, proc(maxproc), ncalls(maxproc)
  real(8):: updates, musppt, calls
  real(8), save:: n_previous=0
  real(8):: wc
  integer:: i, n_recv
  !.............................................................................
  if (io_unit%do_validate) return
  !-----------------------------------------------------------------------------
  ! Prevent several treads from printing at the same time
  !-----------------------------------------------------------------------------
  wc = wallclock()
  sec = wc - prev2
  if (sec < 0.01d0) return
  total = 0d0
  calls = 0d0
  do i=1,ntimer
    proc(i) = sum(timers(:)%total(i))
    ncalls(i) = sum(timers(:)%calls(i))
    total = total + proc(i)
    calls = calls + ncalls(i)
  end do
  if (io_unit%master) then
    print 4, 'procedure', 'calls', 'time', 'percent', 'time/call'
    if (omp%master .and. io_unit%verbose > 0) &
    write (io_unit%mpi,4) 'procedure', 'calls', 'time', 'percent', 'time/call'
    4 format(23x,a,7x,a,10x,a,7x,a,5x,a)
    do i=1,ntimer
      if (proc(i)/total > 5e-4) then
        print 1, trim(names(i)), ncalls(i), proc(i), proc(i)/total*100., proc(i)/max(ncalls(i),1d0)
        if (omp%master .and. io_unit%verbose > 0) &
        write (io_unit%mpi,1) trim(names(i)), ncalls(i), proc(i), proc(i)/total*100., proc(i)/max(ncalls(i),1d0)
        1 format(a32,1p,2g15.3,0p,f10.1,f14.6)
      end if
    end do
  end if
  updates = max(1.0,real(n_update-n_previous))
  !-----------------------------------------------------------------------------
  ! 'sec' measure wall clock time on a single thread, while 'total' is a sum over
  ! all procedures being measured, over all threads.
  ! 'updates' is the total number of active-cell updates (i.e. not including
  ! ghost cells) across *all* threads.
  ! `musppt` is defined to be the number of microseconds per update per thread.
  !-----------------------------------------------------------------------------
  musppt = sec*1d6*min(omp%nthreads,omp%ncores)/updates
  n_recv = max(timer%n_recv,1)
  if (io_unit%master) then
    print 2, "TOTAL thread time, calls", calls, total, 100., musppt, 'core-mus/cell-upd', &
      sec, 'wall sec'
    2 format(a32,1p,2g15.3,g12.3,2(g12.3,2x,a,5x))
    if (mpi%size>1) then
      timer%latency%aver = timer%latency%aver/max(1,timer%latency%n)
      print 3, timer%bytes_recv/1024.**2/sec, timer%bytes_recv/1024.**2/n_recv, &
        timer%nq_send_max, timer%nq_recv_max, timer%latency%max, timer%latency%aver, timer%n_test
      3 format(8x,"MPI recv:",f10.1," MB/s",f11.3," MB/mesg",i8," nq_send_max",i8, &
        " nq_recv_max",2f8.3," max, aver latency",i8, " n_test")
      timer%latency%max  = 0.0
      timer%latency%aver = 0d0
      timer%latency%n    = 0
      timer%n_test = 0_8
    end if
  end if
  if (omp%master .and. io_unit%verbose>0) then
    write (io_unit%mpi,2) "TOTAL thread time, calls", calls, total, 100., musppt, 'core-mus/cell-upd', &
      sec, 'wall sec'
    write (io_unit%mpi,3) timer%bytes_recv/1024.**2/sec, timer%bytes_recv/1024.**2/n_recv, &
      timer%nq_send_max, timer%nq_recv_max
  end if
  timer%bytes_recv = 0.0_8
  timer%n_recv = 0_8
  timer%nq_send_max = 0
  timer%nq_recv_max = 0
  if (omp%master .and. io_unit%verbose>0) then
    write (io_unit%mpi,'(a,1p,e9.2,3x,a,g10.3)') 'timer_mod::print, updates =', &
      updates,'mus/p =', musppt
    flush (io_unit%mpi)
  end if
  n_previous = n_update
  prev2 = wc
  do i=1,ntimer
    timers(:)%total(i) = 0d0
    timers(:)%calls(i) = 0
  end do
END SUBROUTINE print

!> =============================================================================
!> =============================================================================
SUBROUTINE tic (time)
  implicit none
  real(8), optional:: time
  !.............................................................................
  real(8):: omp_get_wtime
  if (present(time)) then
    time = omp_get_wtime()
  else
    wt = omp_get_wtime()
  end if
END SUBROUTINE tic

!> =============================================================================
!> =============================================================================
SUBROUTINE toc_1 (label, time)
  implicit none
  character(len=*):: label
  real(8), optional:: time
  integer:: n
  !.............................................................................
  call toc_r (label, 1.0, time)
END SUBROUTINE toc_1

!> =============================================================================
!> =============================================================================
SUBROUTINE toc_i4 (label, n, time)
  implicit none
  character(len=*):: label
  real(8), optional:: time
  integer:: n
  !.............................................................................
  call toc_r (label, real(n), time)
END SUBROUTINE toc_i4

!> =============================================================================
!> =============================================================================
SUBROUTINE toc_i8 (label, n, time)
  implicit none
  character(len=*):: label
  real(8), optional:: time
  integer(8):: n
  !.............................................................................
  call toc_r (label, real(n), time)
END SUBROUTINE toc_i8

!> =============================================================================
FUNCTION str_f4 (t) RESULT (out)
  real(8):: t
  character(len=8):: out
  !-----------------------------------------------------------------------------
  write (out,'(f8.3)') t
  out = adjustl(out)
  out = trim(out(1:4))
END FUNCTION

!> =============================================================================
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

!> =============================================================================
!> Print the time used per thread.  Overuse is compensated for by the factor f =
!> the number of threads per core.
!> =============================================================================
SUBROUTINE toc_r (label, n, time)
  implicit none
  character(len=*):: label
  real:: n
  real(8), optional:: time
  real(8):: omp_get_wtime, dtime
  integer:: n_cores, n_threads
  !.............................................................................
  if (io_unit%do_validate) return
  if (present(time)) then
    dtime = (omp_get_wtime() - time)
  else
    dtime = (omp_get_wtime() - wt)
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
    time = omp_get_wtime()
  else
    wt = omp_get_wtime()
  end if
END SUBROUTINE toc_r

END MODULE timer_mod
