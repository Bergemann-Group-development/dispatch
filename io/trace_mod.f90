!===============================================================================
!> $Id$
!===============================================================================
MODULE trace_mod
  USE omp_mod
  USE omp_timer_mod
  USE io_unit_mod
  !USE mpi_mod
  USE timer_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! The public instance of the trace data type allows trace_begin calls to be
  ! replaced by trace%begin calls, for consistency
  !-----------------------------------------------------------------------------
  type trace_t
    integer:: indent=4
  contains
    procedure, nopass:: begin => trace_begin
    procedure, nopass:: tag   => trace_tag
    procedure, nopass:: end   => trace_end
    procedure, nopass:: print_id
    !procedure, nopass:: print_hl
    procedure, nopass:: back
  end type
  type(trace_t), public:: trace
  character(len=120), public:: pad_trace=''
  !-----------------------------------------------------------------------------
  ! The data manipulated by the trace routines are thread private, and do not
  ! need to be visible from the calling code.
  !-----------------------------------------------------------------------------
  integer, parameter:: maxlev=30                          ! max number of levels
  character(len=80), dimension(maxlev), save:: tracing    ! active routine
  integer, save:: verbosity(maxlev)                       ! verbosity level
  integer, public, save:: indent=4                        ! trace indent
  integer, private, save:: level=1                        ! trace level
  !$omp threadprivate(indent,level,tracing,verbosity)
PUBLIC trace_begin, trace_tag, trace_end
CONTAINS

!===============================================================================
!> Start a new trace level, and a new timer epoch.  A call that has a
!> detailed_timer logical present will call the timer level only if that
!> logical is true, while call that doesn't have that optional argument
!> will call the timer if and only if the itimer argument is present.
!===============================================================================
SUBROUTINE trace_begin (id, set_verbose, itimer, detailed_timer)
  implicit none
  integer, optional:: set_verbose, itimer
  logical, optional:: detailed_timer
  character(len=*) id
  character(len=mch):: fmt
  integer:: verbose
  real(8):: wc
  !-----------------------------------------------------------------------------
  tracing(level) = id
  if (present(set_verbose)) then
    verbosity(level) = set_verbose
  else
    verbosity(level) = -1
  end if
  verbose = verbosity(level)
  if (level<maxlev) then
    level = level+1
    indent = indent+2
  end if
  if (io_unit%do_trace) then
    write(fmt,'(a,i3.3,a)') '("trace:",i3,f12.6,',indent,'x,a,i5)'
    if (io_unit%verbose >= verbose) then
      wc = wallclock()
      if (io_unit%do_validate) wc = 0d0
      if (io_unit%omp_trace) then
        call print_out (io_unit%log)
      else
        !$omp critical (trace_cr)
        call print_out (io_unit%output)
        !$omp end critical (trace_cr)
      end if
    end if
  end if
  !-----------------------------------------------------------------------------
  if (present(itimer)) then
    if (present(detailed_timer)) then
      if (detailed_timer) then
        call timer%begin (id, itimer)
      end if
    else
      call timer%begin (id, itimer)
    end if
  end if
  !-----------------------------------------------------------------------------
contains
!===============================================================================
subroutine print_out (unit)
  integer:: unit
  if (present(itimer)) then
    write(unit,fmt) &
      omp%thread, wc, trim(id)//' begin, itimer =', itimer
  else
    write(unit,fmt) &
      omp%thread, wc, trim(id)//' begin'
  end if
  flush(unit)
end subroutine
END SUBROUTINE trace_begin

!===============================================================================
!===============================================================================
SUBROUTINE trace_tag (id)
  implicit none
  character(len=*) id
  character(len=mch):: fmt
  integer:: verbose
  real(8):: wc
  !.............................................................................
  if (.not.io_unit%do_trace) return
  verbose = verbosity(level)
  if (io_unit%verbose>=verbose) then
    write(fmt,'(a,i3.3,a)') '("trace:",i3,f12.6,',indent,'x,a)'
    wc = wallclock()
    level  = max(level-1,1)
    if (io_unit%omp_trace) then
      write(io_unit%log,fmt) omp%thread, wc, trim(tracing(level))//' '//trim(id)
      flush(io_unit%log)
    else
      !$omp critical (trace_cr)
      write(io_unit%output,fmt) omp%thread, wc, trim(tracing(level))//' '//trim(id)
      flush(io_unit%output)
      !$omp end critical (trace_cr)
    end if
    level  = min(level+1,maxlev)
  end if
END SUBROUTINE trace_tag

!===============================================================================
!> End a trace level.
!===============================================================================
SUBROUTINE trace_end (itimer, detailed_timer)
  implicit none
  integer, optional:: itimer
  logical, optional:: detailed_timer
  integer:: i, ind
  character(len=mch):: fmt
  integer:: verbose
  real(8):: wc
  !-----------------------------------------------------------------------------
  if (present(itimer)) then
    if (present(detailed_timer)) then
      if (detailed_timer) then
        call timer%end(itimer)
      end if
    else
      call timer%end(itimer)
    end if
  end if
  level  = max(level-1,1)
  ind    = max(indent,4)
  indent = max(indent-2,4)
  if (.not.io_unit%do_trace) return
  verbose = verbosity(level)
  if (io_unit%verbose>=verbose) then
    write(fmt,'(a,i3.3,a)') '("trace:",i3,f12.6,',ind,'x,a,i5)'
    wc = wallclock()
    if (io_unit%do_validate) wc = 0d0
    if (io_unit%omp_trace) then
      call print_out (io_unit%log)
    else
      !$omp critical (trace_cr)
      call print_out (io_unit%output)
      !$omp end critical (trace_cr)
    end if
  end if
contains
!===============================================================================
subroutine print_out (unit)
  integer:: unit
  if (present(itimer)) then
    write (unit,fmt) &
      omp%thread, wc, trim(tracing(level))//' end, itimer =', itimer
  else
    write (unit,fmt) &
      omp%thread, wc, trim(tracing(level))//' end'
  end if
  flush (unit)
end subroutine
END SUBROUTINE trace_end

!===============================================================================
!> Write a trace-back to stdout
!===============================================================================
SUBROUTINE back
  integer:: l
  !-----------------------------------------------------------------------------
  write (stdout,'(/," thread",i4,2x,"traceback:")') omp%thread
  do l=1,level-1
    write (stdout,'(i3,2x,a)') l,trim(tracing(l))
  end do
END SUBROUTINE back

!===============================================================================
SUBROUTINE print_hl
  character(len=120), save:: hl= &
    '--------------------------------------------------------------------------------'
  !..............................................................................
  if (io_unit%master) then
    write (io_unit%output,'(a)') trim(hl)
  end if
END SUBROUTINE print_hl

!===============================================================================
SUBROUTINE print_id (id)
  character(len=120) id
  character(len=120), save:: hl= &
    '--------------------------------------------------------------------------------'
  !..............................................................................
  if (id /= '') then
    !$omp critical (print_id_cr)
    if (id /= '' .and. len_trim(id) > 8) then
      write (io_unit%hash,'(a)') id(6:len_trim(id)-2)
      flush (io_unit%hash)
      id = ''
    end if
    !$omp end critical (print_id_cr)
  end if
END SUBROUTINE print_id

END MODULE trace_mod
