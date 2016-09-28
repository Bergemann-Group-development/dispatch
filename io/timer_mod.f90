!> =============================================================================
!> =============================================================================
MODULE timer_mod
  USE mpi_mod
  USE io_unit_mod
  USE omp_mod
  USE omp_timer_mod
  implicit none
  private
  integer, parameter:: maxproc=100
  integer, save:: ntimer=0
  type:: timer_t
    integer(8):: calls(maxproc)=0
    real(8), dimension(maxproc):: start, total
  contains
    procedure:: init
    procedure:: begin
    procedure:: end
    procedure, nopass:: print
  end type
  character(len=64), dimension(maxproc):: names
  type(timer_t), pointer:: timers(:)
  type(timer_t), public:: timer
  integer, save:: nprint=20
  real(8), save:: sec_per_report=10d0, sec=0d0, next_report_sec
  real(8), save:: prev1=0d0
  real(8), save:: prev2=0d0
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self)
  class(timer_t):: self
  logical, save:: first_time=.true.
  integer:: ignore
  namelist /timer_params/ sec_per_report
  !.............................................................................
  !$omp critical (timer_cr)
  if (first_time) then
    first_time = .false.
    rewind (io_unit%input)
    read (io_unit%input, timer_params, iostat=ignore)
    if (io_unit%master) write (*, timer_params)
    next_report_sec = sec_per_report
    allocate (timers(0:omp%nthreads-1))
  end if
  !$omp end critical (timer_cr)
END SUBROUTINE init

!> =============================================================================
!> =============================================================================
SUBROUTINE begin (self, name, itimer)
  class(timer_t):: self
  character(len=*):: name
  integer:: itimer, thread, i
  logical:: ok
  !.............................................................................
  if (io_unit%do_validate) return
  thread = omp%mythread()
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
  thread = omp%mythread()
  if (itimer > 0) then
    wc = wallclock()
    timers(thread)%calls(itimer) = timers(thread)%calls(itimer) + 1
    timers(thread)%total(itimer) = timers(thread)%total(itimer) + &
      wc - timers(thread)%start(itimer)
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
      print 1, trim(names(i)), ncalls(i), proc(i), proc(i)/total*100., proc(i)/max(ncalls(i),1d0)
      1 format(a32,1p,2g15.3,2g12.3)
    end do
  end if
  updates = max(1.0,real(n_update-n_previous))
  musppt = min(omp%ncores,omp%nthreads)*sec*1d6/updates
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

END MODULE timer_mod
