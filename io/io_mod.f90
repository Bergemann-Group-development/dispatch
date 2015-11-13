!> =============================================================================
!> $Id$
!> =============================================================================
MODULE io_mod
  USE mpi_basic
  USE mpi_timer
  USE mpi_mod
  implicit none
PRIVATE
PUBLIC init_io, print_id, do_debug, do_trace, verbose, input, output, master, &
       mpi_size
CONTAINS

!> =============================================================================
SUBROUTINE init_io
  character(len=mch), save:: id = &
  '$Id$'
  namelist /io_params/ verbose, do_debug, do_trace
  !.............................................................................
  call print_id (id)
  read (*,io_params)
END SUBROUTINE init_io

!> =============================================================================
SUBROUTINE print_id (id)
  character(len=mch) id
  character(len=mch), save:: hl= &
    '--------------------------------------------------------------------------------'
  !..............................................................................
  if (id .ne. ' ' .or. mpi_trace) then
    if (master) then
      print'(a)',trim(hl)
      print'(1x,a,f12.2)', trim(id), wallclock()
    end if
    if (.not.mpi_trace) id = ' '
  end if
END SUBROUTINE print_id

END MODULE io_mod
