!> =============================================================================
!> $Id$
!> =============================================================================
MODULE trace_mod
  USE omp_mod
  USE omp_timer_mod
  USE io_mod
  USE mpi_mod
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
  !$omp threadprivate(indent,level,tracing,verbosity)
PUBLIC trace_begin, trace_tag, trace_end
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE trace_begin (id, set_verbose)
  implicit none
  integer, optional:: set_verbose
  character(len=*) id
  character(len=mch):: fmt
  integer:: verbose
  !.............................................................................
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
  if (io%master.and.io%verbose>=verbose) then
    write(*,fmt) omp_mythread(), wallclock(), trim(id)//' begin'
    call flush(6)
  end if
  write (io%log_unit,fmt) omp_mythread(), wallclock(), trim(id)//' begin'
  call flush(io%log_unit)
END SUBROUTINE trace_begin

!> =============================================================================
!> =============================================================================
SUBROUTINE trace_tag (id)
  implicit none
  character(len=*) id
  character(len=mch):: fmt
  integer:: verbose, level
  if (.not.io%do_trace) return
  level = max(1,level-1)
  verbose = verbosity(level)
  write(fmt,'(a,i3.3,a)') '("trace:",i3,f12.3,',indent,'x,a)'
  if (io%master.and.io%verbose>=verbose) then
    write(*,fmt) omp_mythread(), wallclock(), trim(tracing(level))//' '//trim(id)
    call flush(6)
  end if
END SUBROUTINE trace_tag

!> =============================================================================
!> =============================================================================
SUBROUTINE trace_end
  implicit none
  integer i
  character(len=mch):: fmt
  integer:: verbose
  if (.not.io%do_trace) return
  level  = max(level-1,1)
  indent = max(indent,4)
  verbose = verbosity(level)
  write(fmt,'(a,i3.3,a)') '("trace:",i3,f12.3,',indent,'x,a,i5)'
  if (io%master.and.io%verbose>=verbose) then
    write(*,fmt) omp_mythread(), wallclock(), trim(tracing(level))//' end'
  end if
  call flush(6)
  write (io%log_unit,fmt) omp_mythread(), wallclock(), trim(tracing(level))//' end'
  call flush(io%log_unit)
  indent = max(indent-3,4)
END SUBROUTINE trace_end

END MODULE trace_mod
