!> =============================================================================
!> $Id$
!> =============================================================================
MODULE io_mod
  USE mpi_mod
  USE mpi_timer
  implicit none
  logical:: mpi_trace=.false.
  type io_t
    integer:: verbose=0, input=1, output=2, data_unit=3, log_unit=4
    logical:: master, do_debug=.false., do_trace=.false.
  contains
    procedure:: init  => init_io
    procedure:: print => print_id_
  end type
  integer, parameter:: mch=120
  logical:: first_time=.true.
  type(io_t):: io
PRIVATE
PUBLIC io_t, io, print_id, mch
CONTAINS

!> =============================================================================
SUBROUTINE init_io (self)
  class(io_t):: self
  character(len=120), save:: id = &
  'io_mod.f90 $Id$'
  integer:: verbose
  character(len=mch):: inputfile
  logical:: do_debug, do_trace
  namelist /io_params/ verbose, do_debug, do_trace
  !.............................................................................
  call print_id (id)
  call getarg(1,inputfile)                                                      ! 1st command line argument
  if (inputfile==' ') inputfile = 'input.nml'                                   ! default = input.nml
  do_trace = io%do_trace
  do_debug = io%do_debug
  verbose  = io%verbose
  open (io%input, file=inputfile, form='formatted', status='old')                  ! open for reading
  rewind (io%input)
  read (io%input, io_params)
  io%master = mpi%rank==0
  io%do_debug = do_debug
  io%do_trace = do_trace
  io%verbose = verbose
  if (io%master.and.first_time) write (*,io_params)
  self%master = mpi%rank==0
  self%do_debug = do_debug
  self%do_trace = do_trace
  self%verbose = verbose
  first_time = .false.
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
