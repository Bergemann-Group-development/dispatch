!> =============================================================================
!> $Id$
!> =============================================================================
MODULE io_mod
  USE omp_timer_mod
  implicit none
  logical:: mpi_trace=.false.
  type io_t
    integer:: verbose=0, input=1, output=2, data_unit=3, log_unit=4
    logical:: master=.true., do_debug=.false., do_trace=.false.
    integer:: levelmax=10
  contains
    procedure:: init
    procedure:: print_id => print_id_
    procedure, nopass:: print_hl
  end type
  integer, parameter:: mch=120
  type(io_t):: io
  character(len=64), target:: experimentname
PRIVATE
PUBLIC io_t, io, print_id, mch, experimentname
CONTAINS

!> =============================================================================
SUBROUTINE init (self, name)
  class(io_t):: self
  character(len=120), save:: id = &
  'io_mod.f90 $Id$'
  integer:: verbose=0, levelmax=10
  character(len=*), optional:: name
  character(len=mch):: inputfile, argfile
  logical:: first_time=.true.
  logical:: do_debug=.false., do_trace=.false., exist
  namelist /io_params/ verbose, do_debug, do_trace, levelmax
  !.............................................................................
  call print_id (id)
  if (first_time) then
    first_time=.false.
    if (present(name)) inputfile = trim(name)//'.nml'                             ! experiment_name given?
    inquire (file=inputfile, exist=exist)                                         ! experiment.nml ?
    if (.not.exist) inputfile = 'input.nml'                                       ! default
    call getarg(1,argfile)                                                        ! command line argument
    if (argfile /= ' ') inputfile = argfile                                       ! use, if present
    do_trace = io%do_trace
    do_debug = io%do_debug
    verbose  = io%verbose
    open (io%input, file=inputfile, form='formatted', status='old')               ! open for reading
    experimentname = inputfile(1:index(inputfile,'.')-1)
    if (io%master) then
      print*,'======================================================================='
      print*,'NOTE: Reading parameters from '//trim(inputfile)
      print*,'======================================================================='
    end if
    rewind (io%input)
    read (io%input, io_params)
  end if
  io%do_debug = do_debug
  io%do_trace = do_trace
  io%verbose  = verbose
  io%levelmax = levelmax
  if (io%master.and.first_time) write (*,io_params)
  self%do_debug = do_debug
  self%do_trace = do_trace
  self%verbose  = verbose
END SUBROUTINE init

!> =============================================================================
SUBROUTINE print_id_ (self, id)
  class(io_t):: self
  character(len=120) id
  !..............................................................................
  call print_id (id)
END SUBROUTINE print_id_

!> =============================================================================
SUBROUTINE print_hl
  character(len=120), save:: hl= &
    '--------------------------------------------------------------------------------'
  !..............................................................................
  if (io%master) then
    print'(a)',trim(hl)
  end if
END SUBROUTINE print_hl

!> =============================================================================
SUBROUTINE print_id (id)
  character(len=120) id
  character(len=120), save:: hl= &
    '--------------------------------------------------------------------------------'
  !..............................................................................
  if (id .ne. ' ' .or. mpi_trace) then
    if (io%master) then
      call print_hl
      print'(1x,a,f12.2)', trim(id), wallclock()
    end if
    if (.not.mpi_trace) id = ' '
  end if
END SUBROUTINE print_id

END MODULE io_mod
