!> =============================================================================
!> $Id$
!> =============================================================================
MODULE trace_mod
  USE omp_mod
  USE omp_timer_mod
  USE io_mod
  implicit none
  private
  character(len=80), dimension(30), save:: tracing                              ! active routine
  integer, save:: verbosity(30)
  integer, private, save:: indent=4                                             ! trace indent
  integer, private, save:: level=1                                              ! trace level
  !$omp threadprivate(indent,level,tracing,verbosity)
PUBLIC trace_begin, trace, trace_end
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
  tracing(level) = trim(id)
  if (present(set_verbose)) then
    verbosity(level) = set_verbose
  else
    verbosity(level) = 0
  end if
  verbose = verbosity(level)
  indent = indent + 3
  level = level+1
  if (io%do_trace) then
    write(fmt,'(a,i3.3,a)') '("trace:",i3,f12.3,',indent,'x,a,i5)'
    if (io%master.and.io%verbose>=verbose) then
      write(*,fmt) omp_mythread(), wallclock(), trim(id)//' begin'
    end if
  end if
END SUBROUTINE trace_begin

!> =============================================================================
!> =============================================================================
SUBROUTINE trace (id)
  implicit none
  character(len=*) id
  character(len=mch):: fmt
  integer:: verbose, level
  level = max(1,level-1)
  verbose = verbosity(level)
  if (io%do_trace) then
    write(fmt,'(a,i3.3,a)') '("trace:",i3,f12.3,',indent,'x,a)'
    if (io%master.and.io%verbose>=verbose) then
      write(*,fmt) omp_mythread(), wallclock(), trim(tracing(level))//' '//trim(id)
    end if
  end if
END SUBROUTINE trace

!> =============================================================================
!> =============================================================================
SUBROUTINE trace_end
  implicit none
  integer i
  character(len=mch):: fmt
  integer:: verbose
  level  = max(level-1,1)
  indent = max(indent,4)
  verbose = verbosity(level)
  if (io%do_trace) then
    write(fmt,'(a,i3.3,a)') '("trace:",i3,f12.3,',indent,'x,a,i5)'
    if (io%master.and.io%verbose>=verbose) then
      write(*,fmt) omp_mythread(), wallclock(), trim(tracing(level))//' end'
    end if
    call flush(6)
  end if
  indent = max(indent-3,4)
END SUBROUTINE trace_end

END MODULE trace_mod
