!> =============================================================================
!> $Id$
!> =============================================================================
MODULE io_mod
  USE mpi_mod
  USE mpi_timer
  implicit none
  logical:: mpi_trace=.false.
  type io_t
    integer:: verbose, input, output
    logical:: master, do_debug, do_trace
  contains
    procedure:: init  => init_io
    procedure:: print => print_id_
  end type
  integer, parameter:: mch=120
  type(io_t):: io
PRIVATE
PUBLIC io_t, io, print_id, mch
CONTAINS

!> =============================================================================
SUBROUTINE init_io (self)
  class(io_t):: self
  character(len=120), save:: id = &
  '$Id$'
  integer:: verbose
  logical:: do_debug, do_trace
  namelist /io_params/ verbose, do_debug, do_trace
  !.............................................................................
  call print_id (id)
  read (*,io_params)
  self%do_debug = do_debug
  self%do_trace = do_trace
  self%verbose = verbose
  io = self
END SUBROUTINE init_io

!> =============================================================================
SUBROUTINE print_id_ (self, id)
  class(io_t):: self
  character(len=120) id
  !..............................................................................
  call print_id (id)
END SUBROUTINE

!> =============================================================================
SUBROUTINE print_id (id)
  character(len=120) id
  character(len=120), save:: hl= &
    '--------------------------------------------------------------------------------'
  !..............................................................................
  if (id .ne. ' ' .or. mpi_trace) then
    if (mpi%master) then
      print'(a)',trim(hl)
      print'(1x,a,f12.2)', trim(id), wallclock()
    end if
    if (.not.mpi_trace) id = ' '
  end if
END SUBROUTINE print_id

END MODULE io_mod
