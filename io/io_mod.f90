!> =============================================================================
!> $Id$
!> =============================================================================
MODULE io_mod
  USE omp_mod
  USE omp_timer_mod
  USE mpi_mod
  USE io_unit_mod
  implicit none
  private
  logical:: mpi_trace=.false.
  type io_t
    integer:: verbose=0, input, output, data_unit, log_unit
    logical:: master=.true., do_debug=.false., do_trace=.false., do_output=.true.
    logical:: do_flags, do_legacy
    integer:: levelmax
    integer:: word_size=4
    character(len=64):: inputname, outputname, datadir='data'
    character(len=1):: sep = '/'
  contains
    procedure:: init
    procedure:: print_id => print_id_
    procedure, nopass:: print_hl
  end type
  character(len=64), save:: inputname, outputname
  type(io_t), public:: io
  public print_id, mch,io_unit
CONTAINS

!> =============================================================================
!> Initialize I/O related parameters.  The first thread that enters here blocks
!> access in a critical region, reads namelist parameters, and stores them in a
!> global, public instance, which we normally refer to.  Other calls pick up the
!> global values and store them in other instances (which may be the same).
!> =============================================================================
SUBROUTINE init (self, name)
  class(io_t):: self
  character(len=120), save:: id = &
  'io_mod.f90 $Id$'
  integer:: verbose=0, levelmax=10, iostat
  character(len=*), optional:: name
  character(len=64):: filename
  logical:: first_time=.true.
  logical:: do_debug=.false., do_trace=.false., do_output=.false., exist, &
    do_legacy=.true., do_flags=.false.
  namelist /io_params/ verbose, do_debug, do_trace, do_output, do_flags, &
    do_legacy, levelmax
  real:: test
  !.............................................................................
  call init_omp
  self%master = self%master .and. (mpi%rank==0)
  io%master = self%master
  call print_id (id)
  !-----------------------------------------------------------------------------
  ! I/O unit numbers, for convenience and backwards compatibility.  All other
  ! unit numbers should be taken directly from io_unit_mod
  !-----------------------------------------------------------------------------
  self%input     = io_unit%input
  self%output    = io_unit%output
  self%log_unit  = io_unit%log
  self%data_unit = io_unit%data
  !-----------------------------------------------------------------------------
  ! Use critical region to avoid threads using the input file simultaneously
  !-----------------------------------------------------------------------------
  !$omp critical (input_cr)
  if (first_time) then
    !---------------------------------------------------------------------------
    ! The input filename is input.nml, or as given on the command line
    !---------------------------------------------------------------------------
    call getarg(1,filename)                                     ! command line
    if (filename /= ' ') then
      self%inputname = filename                                 ! run_name.nml
      self%outputname = trim(self%datadir)//self%sep// &        ! data/run_name/
        trim(filename(1:index(filename,'.')-1))//self%sep
    else
      self%inputname = 'input.nml'
      self%outputname = trim(self%datadir)//self%sep            ! data/
    end if
    call system ('mkdir -p '//trim(self%outputname))
    !---------------------------------------------------------------------------
    ! Open an MPI rank specific log file
    !---------------------------------------------------------------------------
    write (filename,'(a,i4.4,".log")') trim(self%outputname), mpi%rank
    if (self%master) print*, 'logfile:', filename
    open (self%log_unit, file=filename, form='formatted', status='unknown')
    open (self%input, file=self%inputname, form='formatted', status='old')
    if (self%master) then
      print*,'======================================================================='
      print*,'NOTE: Reading parameters from '//trim(self%inputname)
      print'(a,i1)', ' This version was compiled with default real KIND=', kind(test)
      print*,'======================================================================='
    end if
    rewind (io_unit%input)
    read (io_unit%input, io_params,iostat=iostat)
    if (self%master.and.first_time) write (*,io_params)
    first_time=.false.
  end if
  !$omp end critical (input_cr)
  self%verbose   = verbose
  self%do_debug  = do_debug
  self%do_trace  = do_trace
  self%do_output = do_output
  self%do_flags  = do_flags
  self%do_legacy = do_legacy
  self%levelmax  = levelmax
  self%inputname = inputname
  self%master = self%master .and. (mpi%rank==0)
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
