!> =============================================================================
!> =============================================================================
MODULE timer_mod
  USE mpi_mod
  USE io_unit_mod
  USE omp_mod
  USE omp_timer_mod
  USE omp_lib
  implicit none
  private
  integer, parameter:: maxproc=100
  integer, save:: ntimer=0
  type:: timer_t
    integer(8):: calls(maxproc)=0
    real(8), dimension(maxproc):: start=0d0, total=0d0
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
  real(8), save:: sec_per_report=10d0, sec=0d0, next_report_sec
  real(8), save:: prev1=0d0
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
  integer, save:: ncores=0
  namelist /timer_params/ sec_per_report, ncores
  !.............................................................................
  !$omp critical (timer_cr)
  if (first_time) then
    first_time = .false.
    ncores = omp%ncores
    rewind (io_unit%input)
    read (io_unit%input, timer_params, iostat=ignore)
    if (io_unit%master) write (*, timer_params)
    next_report_sec = sec_per_report
    omp%nthreads = max(1,omp%nthreads)
    allocate (timers(0:omp%nthreads-1))
  end if
  !$omp end critical (timer_cr)
  omp%ncores = ncores
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
  integer:: itimer, thread, i
  logical:: ok
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
  end if
  timers(thread)%start(itimer) = wallclock()
END SUBROUTINE begin

!> =============================================================================
!> =============================================================================
SUBROUTINE end (self, itimer)
  class(timer_t):: self
  integer:: itimer, thread, i
  optional:: itimer
  real(8):: wc
  !.............................................................................
  if (io_unit%do_validate) return
  thread = omp_get_thread_num()
  if (itimer > 0) then
    wc = wallclock()
    timers(thread)%calls(itimer) = timers(thread)%calls(itimer) + 1
    timers(thread)%total(itimer) = timers(thread)%total(itimer) + &
      wc - timers(thread)%start(itimer)
    timers(thread)%start(itimer) = wc
    if (io_unit%master) then
      sec = wc - prev1
      if (sec > sec_per_report) then
        call print()
        prev1 = wc
      end if
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
  integer:: i
  !.............................................................................
  !if (.not. io_unit%master) return
  if (io_unit%do_validate) return
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
    print '(23x,a,7x,a,10x,a,7x,a,5x,a)', &
      'procedure', 'calls', 'time', 'percent', 'time/call'
    do i=1,ntimer
      if (proc(i)/total > 1e-4) then
        print 1, trim(names(i)), ncalls(i), proc(i), proc(i)/total*100., proc(i)/max(ncalls(i),1d0)
        1 format(a32,1p,2g15.3,2g12.3)
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
  if (io_unit%master) &
    print 2, "TOTAL thread time, calls", calls, total, 100., musppt, 'core-mus/cell-upd', &
      sec, 'wall sec'
    2 format(a32,1p,2g15.3,g12.3,2(g12.3,2x,a,5x))
  if (io_unit%verbose>1) &
    write (io_unit%log,'(a,1p,e9.2,3x,a,g10.3)') 'timer_mod::print, updates =', &
      updates,'mus/p =', musppt
  n_previous = n_update
  prev1 = wc
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
  integer:: n_cores
  !.............................................................................
  if (io_unit%do_validate) return
  if (present(time)) then
    dtime = (omp_get_wtime() - time)
  else
    dtime = (omp_get_wtime() - wt)
  end if
  n_cores = max(min(omp%nthreads,omp%ncores),1)
  write (io_unit%log,1) label, trim(time_str(dtime)), &
    trim(time_str(dtime*n_cores/max(n,1.))), '/upd,', n, 'updates/process,', &
    mpi%size, merge('processes,','process,  ',mpi%size>1), &
    n_cores, merge('cores/process','core/process ',n_cores>1)
  1 format (1x, a, ':  ', a, ',  ', a, a, 1p, e10.2, 1x, a, 3(i4, 1x, a))
  flush (io_unit%log)
  if (mpi%rank==0.and.omp%thread==0) then
    if (n>1.0) then
      print 1, label, trim(time_str(dtime)), &
        trim(time_str(dtime*n_cores/max(n,1.))), '/upd,', n, 'updates/process,', &
        mpi%size, merge('processes,','process,  ',mpi%size>1), &
        n_cores, merge('cores/process','core/process ',n_cores>1)
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
