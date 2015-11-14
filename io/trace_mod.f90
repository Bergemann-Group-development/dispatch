!> =============================================================================
!> $Id$
!> =============================================================================
MODULE trace_mod
  USE io_mod
  USE mpi_timer
  USE mpi_basic
  implicit none
  character(len=64), dimension(30), save:: tracing                              ! active routine
  integer, private, save:: indent=4                                             ! trace indent
  integer, private, save:: level=1                                              ! trace level
  !$omp threadprivate(indent,level,tracing)
PRIVATE
PUBLIC trace_begin, trace, trace_end, verbose
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE trace_begin (id)
  implicit none
  character(len=*) id
  character(len=mch):: fmt
  !.............................................................................
  tracing(level) = trim(id)
  indent = indent + 3
  level = level+1
  !write(output,*) 'indent+:', indent, level, omp_mythread()
  if (do_trace .or. do_debug) then
    write(fmt,'(a,i3.3,a)') '("trace:",i3,f12.3,',indent,'x,a)'
    if (master.and.verbose>0 .or. verbose>1) then
      write(*,fmt) omp_mythread(), wallclock(), trim(id)//' begin'
    end if
    write(output,fmt) omp_mythread(), wallclock(), trim(id)//' begin'
    call flush (output)
  end if
END SUBROUTINE trace_begin

!> =============================================================================
!> =============================================================================
SUBROUTINE trace (id)
  implicit none
  character(len=*) id
  character(len=mch):: fmt
  if (do_trace .or. do_debug) then
    write(fmt,'(a,i3.3,a)') '("trace:",i3,f12.3,',indent,'x,a)'
    if (master.and.verbose>0 .or. verbose>1) then
      write(*,fmt) omp_mythread(), wallclock(), trim(tracing(level-1))//' '//trim(id)
    end if
    write(output,fmt) omp_mythread(), wallclock(), trim(tracing(level-1))//' '//trim(id)
    call flush (output)
  end if
END SUBROUTINE trace

!> =============================================================================
!> =============================================================================
SUBROUTINE trace_end
  implicit none
  integer i
  character(len=mch):: fmt
  level  = max(level-1,1)
  indent = max(indent,4)
  if (do_trace .or. do_debug) then
    write(fmt,'(a,i3.3,a)') '("trace:",i3,f12.3,',indent,'x,a)'
    if (master.and.verbose>0 .or. verbose>1) then
      write(*,fmt) omp_mythread(), wallclock(), trim(tracing(level))//' end'
    end if
    write(output,fmt) omp_mythread(), wallclock(), trim(tracing(level))//' end'
    call flush (output)
  end if
  indent = max(indent-3,4)
  !write(output,*) 'indent-:', indent, level, omp_mythread()
END SUBROUTINE trace_end

!> =============================================================================
!> =============================================================================
SUBROUTINE trace_mpi (message)
  character(len=*):: message
!..............................................................................
  if (mpi_trace) print*,mpi_rank,' trace: '//trim(message)
END SUBROUTINE trace_mpi

END MODULE trace_mod
