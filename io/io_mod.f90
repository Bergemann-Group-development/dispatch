!> =============================================================================
!> $Id$
!> =============================================================================
MODULE io_mod
  USE mpi_basic
  USE mpi_timer
  implicit none
  type io_t
    integer:: verbose, input, output
    logical:: master, do_debug, do_trace
  contains
    procedure:: init  => init_io
    procedure:: print => print_id
  end type
PRIVATE
PUBLIC io_t
CONTAINS

!> =============================================================================
SUBROUTINE init_io (self)
  class(io_t):: self
  character(len=mch), save:: id = &
  '$Id$'
  namelist /io_params/ verbose, do_debug, do_trace
  !.............................................................................
  call print_id (self, id)
  read (*,io_params)
END SUBROUTINE init_io

!> =============================================================================
SUBROUTINE print_id (self, id)
  class(io_t):: self
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
