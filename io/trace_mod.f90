!> =============================================================================
!> $Id$
!> =============================================================================
MODULE trace_mod
  USE omp_mod
  USE omp_timer_mod
  USE io_mod
  USE io_unit_mod
  USE mpi_mod
  USE timer_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! The public instance of the trace data type allows trace_begin calls to be
  ! replaced by trace%begin calls, for consistency
  !-----------------------------------------------------------------------------
  type trace_t
  contains
    procedure, nopass:: begin => trace_begin
    procedure, nopass:: tag   => trace_tag
    procedure, nopass:: end   => trace_end
  end type
  type(trace_t), public:: trace
  !-----------------------------------------------------------------------------
  ! The data manipulated by the trace routines are thread private, and do not
  ! need to be visible from the calling code.
  !-----------------------------------------------------------------------------
  integer, parameter:: maxlev=30                                ! max number of levels
  character(len=80), dimension(maxlev), save:: tracing          ! active routine
  integer, save:: verbosity(maxlev)                             ! verbosity level
  integer, private, save:: indent=4                             ! trace indent
  integer, private, save:: level=1                              ! trace level
  integer, private, save:: otimer=0                             ! previous value
  !$omp threadprivate(indent,level,tracing,verbosity,otimer)
PUBLIC trace_begin, trace_tag, trace_end
CONTAINS

!> =============================================================================
!> Start a new trace level, and a new timer epoch.
!> =============================================================================
SUBROUTINE trace_begin (id, set_verbose, itimer)
  implicit none
  integer, optional:: set_verbose, itimer
  character(len=*) id
  character(len=mch):: fmt
  integer:: verbose
  real(8):: wc
  !-----------------------------------------------------------------------------
  ! If a different timer epoch is still open, as signalled by the previous itimer
  ! value, stored in otimer, then close that timer epoch by calling timer_end
  !-----------------------------------------------------------------------------
  if (present(itimer)) then
    if (itimer/=otimer .and. otimer>0) call timer%end (otimer)
    call timer%begin (id, itimer)
    otimer = itimer
  !-----------------------------------------------------------------------------
  ! If no new itimer value is given, but an old timer epoch is still open, then
  ! close it
  !-----------------------------------------------------------------------------
  else if (otimer > 0) then
    call timer%end (otimer)
  end if
  if (.not.io%do_trace) return
  tracing(level) = trim(id)
  if (present(set_verbose)) then
    verbosity(level) = set_verbose
  else
    verbosity(level) = 0
  end if
  verbose = verbosity(level)
  if (level<maxlev) then
    level = level+1
    indent = indent + 3
  end if
  write(fmt,'(a,i3.3,a)') '("trace:",i3,f12.3,',indent,'x,a,i5)'
  if (io%verbose>=verbose) then
    wc = wallclock()
    if (io_unit%do_validate) wc = 0d0
    if (io%omp_trace) then
      write(io_unit%log,fmt) omp_mythread, wc, trim(id)//' begin'
      flush(io_unit%log)
    else if (io%master) then
      write(io%output,fmt) omp_mythread, wc, trim(id)//' begin'
      flush(io%output)
    end if
  end if
END SUBROUTINE trace_begin

!> =============================================================================
!> =============================================================================
SUBROUTINE trace_tag (id)
  implicit none
  character(len=*) id
  character(len=mch):: fmt
  integer:: verbose
  real(8):: wc
  !.............................................................................
  if (.not.io%do_trace) return
  verbose = verbosity(level)
  if (io%verbose>=verbose) then
    write(fmt,'(a,i3.3,a)') '("trace:",i3,f12.3,',indent,'x,a)'
    wc = wallclock()
    level  = max(level-1,1)
    if (io%omp_trace) then
      write(io_unit%log,fmt) omp_mythread, wc, trim(tracing(level))//' '//trim(id)
      flush(io_unit%log)
    else if (io%master) then
      write(io%output,fmt) omp_mythread, wc, trim(tracing(level))//' '//trim(id)
      flush(io%output)
    end if
    level  = min(level+1,maxlev)
  end if
END SUBROUTINE trace_tag

!> =============================================================================
!> End a trace level.  If the previous trace%begin call had an itimer argument,
!> the value is remembered in otimer, and that timer epoch can thus be closed,
!> even if an itimer argument is not given.
!> =============================================================================
SUBROUTINE trace_end (itimer)
  implicit none
  integer, optional:: itimer
  integer i
  character(len=mch):: fmt
  integer:: verbose
  real(8):: wc
  !-----------------------------------------------------------------------------
  if (present(itimer)) then
    call timer%end(itimer)
  else if (otimer > 0) then
    call timer%end(otimer)
  end if
  otimer = 0
  if (.not.io%do_trace) return
  level  = max(level-1,1)
  indent = max(indent,4)
  verbose = verbosity(level)
  if (io%verbose>=verbose) then
    write(fmt,'(a,i3.3,a)') '("trace:",i3,f12.3,',indent,'x,a,i5)'
    wc = wallclock()
    if (io_unit%do_validate) wc = 0d0
    if (io%omp_trace) then
      write(io_unit%log,fmt) omp_mythread, wc, trim(tracing(level))//' end'
      flush(io_unit%log)
    else if (io%master) then
      write(io%output,fmt) omp_mythread, wc, trim(tracing(level))//' end'
      flush(io%output)
    end if
  end if
  indent = max(indent-3,4)
END SUBROUTINE trace_end

END MODULE trace_mod
