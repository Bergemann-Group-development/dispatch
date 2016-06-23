!> =============================================================================
!> =============================================================================
MODULE timer_mod
  USE io_unit_mod
  USE omp_mod
  USE omp_timer_mod
  USE omp_lock_mod
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
  character(len=32), dimension(maxproc):: names
  type(lock_t):: lock
  type(timer_t), pointer:: timers(:)
  type(timer_t), public:: timer
  integer, save:: nprint=20
  integer, save:: calls_per_report=10000
  integer(8):: calls, next_report
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self)
  class(timer_t):: self
  logical, save:: first_time=.true.
  integer:: ignore
  namelist /timer_params/ calls_per_report
  !.............................................................................
  !$omp critical (timer_cr)
  if (first_time) then
    first_time = .false.
    read (io_unit%input, timer_params, iostat=ignore)
    if (io_unit%master) write (*, timer_params)
    next_report = calls_per_report
    call lock%init('timer')
    allocate (timers(0:omp%nthreads-1))
    print*,'timer::init', omp%nthreads, associated(timers), size(timers)
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
  thread = omp%mythread()
  if (itimer==0) then
    !$omp critical (timer_cr)
    ok = .true.
    do i=1,ntimer
      ok = ok .and. trim(names(i)) /= trim(name)
    end do
    if (ok) then
      ntimer = ntimer+1
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
  real(8):: total, proc(maxproc)
  !.............................................................................
  thread = omp%mythread()
  if (itimer > 0) then
    timers(thread)%total(itimer) = timers(thread)%total(itimer) &
      + wallclock() - timers(thread)%start(itimer)
    timers(thread)%calls(itimer) = timers(thread)%calls(itimer) + 1
    !$omp atomic
    calls = calls + 1
    if (calls > next_report) then
      next_report = next_report + calls_per_report
      call print()
    end if
  end if
END SUBROUTINE end

!> =============================================================================
!> =============================================================================
SUBROUTINE print ()
  real(8):: total, proc(maxproc),ncalls(maxproc)
  integer:: i
  !.............................................................................
  total = 0d0
  do i=1,ntimer
    proc(i) = sum(timers(:)%total(i))
    ncalls(i) = sum(timers(:)%calls(i))
    total = total + proc(i)
  end do
  do i=1,ntimer
    print 1, trim(names(i)), proc(i), proc(i)/total*100., proc(i)/ncalls(i)
    1 format(a32,1p,g15.3,2g10.2)
  end do
  print 1, "TOTAL time, calls", total, real(calls)
END SUBROUTINE print

END MODULE timer_mod
