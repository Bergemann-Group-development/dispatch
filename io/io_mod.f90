!> =============================================================================
!> $Id$
!> =============================================================================
MODULE io_mod
  USE omp_mod
  USE omp_timer_mod
  USE mpi_mod
  USE io_unit_mod
  USE timer_mod
  implicit none
  private
  logical:: mpi_trace=.false.
  type io_t
    integer:: verbose=0, input, output, data_unit, log_unit
    logical:: master=.true.
    logical:: do_legacy, do_direct, omp_trace
    character(len=16):: method='legacy'
    integer:: levelmax
    integer:: word_size=4
    integer:: iodir = -1
    integer:: format=0
    integer:: nwrite=0, ntask=0, ntotal=0
    character(len=64):: inputname, outputname, datadir='data', top, inputdir
    character(len=1):: sep = '/'
    character(len=72):: hl='-------------------------------------------------------------------------'
    character(len=72):: hs='*************************************************************************'
    ! flag values
    logical:: do_trace=.false.
    logical:: do_output=.true.
    logical:: do_debug=.false.
    logical:: do_flags=.false.
    logical:: do_stop=.false.
    logical:: guard_zones=.false.
    logical:: time_derivs=.false.
    integer:: id_debug=0
    integer:: id_track=0
    real(8):: end_time=1d30
    real(8):: out_time=1d0
    real(8):: print_time=0d0
    real:: flag_time=1.
    real:: flag_max=10.
    real:: ampl=0.0
    real(8):: gb=0.0
  contains
    procedure:: init
    procedure:: debug
    procedure:: check_flags
    procedure:: print_id => print_id_
    procedure, nopass:: mkdir
    procedure, nopass:: mkdir_no
    procedure, nopass:: print_hl
    procedure, nopass:: namelist_warning
    procedure, nopass:: abort
  end type
  character(len=64), save:: inputname, outputname
  type(io_t), public:: io
  public print_id, mch, io_unit, log_unit
CONTAINS

!> =============================================================================
!> Function suitable for testing in debug print statements
!> =============================================================================
LOGICAL FUNCTION debug (self, verbose)
  class(io_t):: self
  integer:: verbose
  debug = self%master .and. (self%verbose >= verbose)
END FUNCTION debug

!> =============================================================================
!> Utility to create a new directory
!> =============================================================================
SUBROUTINE mkdir (dir)
#ifdef __INTEL_COMPILER
  USE ifport
#else
  use iso_c_binding
  interface
    function fmkdir(path,mode) bind(c,name="mkdir")
      use iso_c_binding
      integer(c_int) :: fmkdir
      character(kind=c_char,len=1) :: path(*)
      integer(c_int16_t), value :: mode
    end function fmkdir
  end interface
  integer :: i, iter
#endif
  character(len=*):: dir
  logical:: exists
  real, parameter:: ms=100.
  !.............................................................................
#ifdef __INTEL_COMPILER
  inquire (directory=dir, exist=exists)
  if (exists) return
  if (.not.io_unit%do_validate) then
    exists = MakeDirQQ (dir)
    if (exists) then
      print *,omp%thread,' Intel created directory '//dir
    else
      print *,omp%thread,' WARNING: Intel failed to create directory '//dir
    end if
  end if
#else
  !$omp critical (system_cr)
  inquire (file=dir, exist=exists)
  if (.not.exists) then
    do iter=1,10
      i = fmkdir(dir, int(o'772',c_int16_t))
      call mpi%delay (ms=ms)
      inquire (file=dir, exist=exists)
      if (exists) then
        if (.not.io_unit%do_validate) &
          print *, mpi%rank, omp%thread,' C-binding call created directory '//dir
        exit
      else
        if (.not.io_unit%do_validate) &
          print *, mpi%rank, omp%thread,' WARNING: C-binding call failed to create directory '//dir
      end if
      call system ('mkdir -p '//dir)
      call mpi%delay (ms=ms)
      inquire (file=dir, exist=exists)
      if (exists) then
        if (.not.io_unit%do_validate) &
          print *, mpi%rank, omp%thread,' system call created directory '//dir
        exit
      else
        if (.not.io_unit%do_validate) &
          print *, mpi%rank, omp%thread,' WARNING: system call failed to create directory '//dir
      end if
    end do
  end if
  !$omp end critical (system_cr)
#endif
END SUBROUTINE mkdir

!> =============================================================================
!> Make sure the next directory is created ahead of time
!> =============================================================================
SUBROUTINE mkdir_no (iout)
  integer:: iout
  character(len=120):: filename
  !-----------------------------------------------------------------------------
  !$omp master
  write (filename,'(a,i5.5,"/")') trim(io%outputname), iout
  if (io%iodir/=iout) call io%mkdir (trim(filename))
  io%iodir = iout
  !$omp end master
END SUBROUTINE mkdir_no 

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
  integer, save:: verbose=0, levelmax=10
  character(len=*), optional:: name
  character(len=64):: filename, datadir='data'
  character(len=64), save:: top='../../'
  character(len=8), save:: method='legacy'
  integer, save:: id_debug=-1
  logical, save:: first_time=.true., omp_trace=.false., do_validate=.false.
  logical, save:: do_debug=.false., do_trace=.false., do_output=.false., exist, &
    do_legacy=.false., do_flags=.false., do_direct=.false.
  namelist /io_params/ verbose, do_debug, do_trace, do_output, do_flags, &
    do_validate, do_legacy, do_direct, levelmax, omp_trace, id_debug, top, &
    datadir, method
  real:: test
  integer:: iostat
  !.............................................................................
  self%master = self%master .and. (mpi%rank==0)
  io%master = self%master
  !$omp parallel
  io_unit%master = io_unit%master.and.self%master
  !$omp end parallel
  !-----------------------------------------------------------------------------
  ! I/O unit numbers, for convenience and backwards compatibility.  All other
  ! unit numbers should be taken directly from io_unit_mod
  !-----------------------------------------------------------------------------
  self%input     = io_unit%input
  self%output    = io_unit%output
  self%log_unit  = io_unit%log
  self%data_unit = io_unit%data
  !-----------------------------------------------------------------------------
  ! Use serial region to avoid threads using the input file simultaneously
  !-----------------------------------------------------------------------------
  if (first_time) then
    !---------------------------------------------------------------------------
    ! The input filename is input.nml, or as given on the command line
    !---------------------------------------------------------------------------
    call getarg(1,filename)                                     ! command line
    if (filename/=' '.and.trim(filename)/='input.nml') then
      self%inputname = filename                                 ! run_name.nml
      self%outputname = trim(datadir)//self%sep// &             ! data/run_name/
        trim(filename(1:index(filename,'.')-1))//self%sep
    else
      self%inputname = 'input.nml'
      self%outputname = trim(datadir)//self%sep                 ! data/
    end if
    open (self%input, file=self%inputname, form='formatted', status='old')
    rewind (io_unit%input)
    read (io_unit%input, io_params, iostat=iostat)
    if (self%master.and.first_time) write (*,io_params)
    if (mpi%master .and. .not. do_validate) then
      print'(a,i4)', ' n_socket =', omp%nsockets
      print'(a,i4)', '   n_core =', omp%ncores
      print'(a,i4)', ' n_thread =', omp_nthreads
    end if
    !---------------------------------------------------------------------------
    ! OMP specific setup
    !---------------------------------------------------------------------------
    !$omp parallel
    !$omp critical (open_cr)
    if (omp%master) then
      io_unit%verbose = verbose
    else
      io_unit%verbose = -2
    end if
    !---------------------------------------------------------------------------
    ! Open MPI rank and OMP thread-specific log files
    !---------------------------------------------------------------------------
    io_unit%do_validate = do_validate
    call self%mkdir (trim(self%outputname))
    write (filename,'(a,i5.5,"/")') trim(io%outputname), 0
    call self%mkdir (trim(filename))
    if (omp_trace .and. omp%nthreads>1) then
      io_unit%log = 110+omp%mythread()
      write (filename,'(a,"thread_",i5.5,"_",i3.3,".log")') trim(self%outputname), mpi%rank, omp%thread
      open (io_unit%log, file=filename, form='formatted', status='unknown')
    else
      io_unit%log = io_unit%mpi
    end if
    write (filename,'(a,"rank_",i5.5,".nq")') trim(self%outputname), mpi%rank
    open (io_unit%queue, file=filename, form='formatted', status='unknown')
    write (filename,'(a,"rank_",i5.5,".log")') trim(self%outputname), mpi%rank
    open (io_unit%mpi, file=filename, form='formatted', status='unknown')
    write (filename,'(a,"rank_",i5.5,".disp")') trim(self%outputname), mpi%rank
    open (io_unit%dispatcher, file=filename, form='formatted', status='unknown')
    log_unit = io_unit%log
    write (io_unit%log,'(a,i4,2l4)') &
      'io_mod::init io_unit%log, io_unit%master, io_unit%verbose:', &
      io_unit%log, io_unit%master, io_unit%verbose
    !$omp end critical (open_cr)
    !$omp end parallel
    if (verbose>1.or.self%master) &
      print*, 'logfile:', filename
    self%log_unit = io_unit%log
    if (self%master) then
      print*,'======================================================================='
      print*,'NOTE: Reading parameters from '//trim(self%inputname)
      print'(a,i1)', ' This version was compiled with default real KIND=', kind(test)
      print*,'======================================================================='
    end if
    first_time=.false.
  end if
  self%datadir   = datadir
  self%id_debug  = id_debug
  self%verbose   = verbose
  self%do_debug  = do_debug
  self%do_trace  = do_trace
  self%do_output = do_output
  self%do_flags  = do_flags
  self%do_legacy = do_legacy
  self%do_direct = do_direct
  self%omp_trace = omp_trace
  self%levelmax  = levelmax
  self%inputname = inputname
  self%top       = top
  self%method    = method
  if (do_legacy) self%method = 'legacy'
  if (do_direct) self%method = 'direct'
  io_unit%top    = top
  io_unit%do_validate = do_validate
  call print_id (id)
  call timer%init
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
  if (io_unit%do_validate) return
  if (id .ne. ' ' .or. mpi_trace) then
    if (io%master) then
      call print_hl
      print'(1x,a,f12.2)', trim(id) !, wallclock()
    end if
    if (.not.mpi_trace) id = ' '
  end if
END SUBROUTINE print_id

!> =============================================================================
!> Check debug flag files
!> =============================================================================
SUBROUTINE check_flags (self)
  class (io_t)      :: self
  real(8), save     :: last_checked=0d0
  real(8)           :: flag_time, checked
  logical           :: exists
  character(len=80) :: file
  integer, save     :: itimer=0
  !.............................................................................
  call timer%begin ('io_t%check_flags', itimer)
  flag_time = merge(self%flag_time, self%flag_max, self%do_flags)
  !$omp atomic read
  checked = last_checked
  if (wallclock() > checked+flag_time) then
    !$omp critical (flag_cr)
    last_checked = wallclock()
    file = trim(io%outputname)//'flag'
    inquire (file=trim(file), exist=exists)
    if (exists) then
      open (io_unit%flag, file=trim(file), form='formatted', status='old')
      close (io_unit%flag, status='delete')
      ! logical
      call check ('do_trace.flag'       , lvalue=self%do_trace)
      call check ('do_flags.flag'       , lvalue=self%do_flags)
      call check ('do_debug.flag'       , lvalue=self%do_debug)
      call check ('do_output.flag'      , lvalue=self%do_output)
      call check ('stop.flag'           , lvalue=self%do_stop)
      ! integer
      call check ('id_debug.flag'       , ivalue=self%id_debug)
      call check ('id_track.flag'       , ivalue=self%id_track)
      call check ('verbose.flag'        , ivalue=self%verbose)
      ! single prec
      call check ('flag_time.flag'      , rvalue=self%flag_time)
      call check ('flag_max.flag'       , rvalue=self%flag_max)
      call check ('ampl.flag'           , rvalue=self%ampl)
      ! double prec
      call check ('out_time.flag'       , dvalue=self%out_time)
      call check ('print_time.flag'     , dvalue=self%print_time)
      call check ('end_time.flag'       , dvalue=self%end_time)
      call check ('sec_per_report.flag' , dvalue=timer%sec_per_report)
    end if
    !$omp end critical (flag_cr)
  end if
  call timer%end (itimer)
contains
  !=============================================================================
  subroutine check (file, lvalue, ivalue, rvalue, dvalue)
  character(len=*)  :: file
  logical, optional :: lvalue
  integer, optional :: ivalue
  real   , optional :: rvalue
  real(8), optional :: dvalue
  logical           :: exists
  integer           :: iostat
  character(len=80) :: f
  !.............................................................................
  f = trim(io%outputname)//trim(file)
  inquire (file=trim(f), exist=exists)
  if (exists) then
    open (io_unit%flag, file=trim(f), form='formatted', status='old', iostat=iostat)
    if (present(ivalue)) read (io_unit%flag,*) ivalue
    if (present(rvalue)) read (io_unit%flag,*) rvalue
    if (present(dvalue)) read (io_unit%flag,*) dvalue
    if (present(lvalue).and.iostat==0) read (io_unit%flag,*,iostat=iostat) lvalue
    if (iostat/=0) lvalue = .not.lvalue
    close (io_unit%flag, status='delete')
    if (present(lvalue)) write (*,*) 'flag: ', file(1:index(file,'.')-1), ' =', lvalue
    if (present(ivalue)) write (*,*) 'flag: ', file(1:index(file,'.')-1), ' =', ivalue
    if (present(rvalue)) write (*,*) 'flag: ', file(1:index(file,'.')-1), ' =', rvalue
    if (present(dvalue)) write (*,*) 'flag: ', file(1:index(file,'.')-1), ' =', dvalue
  end if
  end subroutine
END SUBROUTINE check_flags

!> =============================================================================
!> Issue message for missing namelist in input file
!> =============================================================================
SUBROUTINE namelist_warning (namelist)
  character(len=*):: namelist
  !.............................................................................
  if (io%master) then
    print *, '*************************************************************************************'
    print *, '*************************************************************************************'
    print *, '           WARNING: namelist '//trim(namelist)//' had read error'
    print *, '*************************************************************************************'
    print *, '*************************************************************************************'
  end if
END SUBROUTINE namelist_warning

!> =============================================================================
!> Interface to mpi%abort, for modules that do not otherwise USE mpi_mod
!> =============================================================================
SUBROUTINE abort (error)
  character(len=*), optional:: error
  !.............................................................................
  call mpi%abort (error)
END SUBROUTINE abort

END MODULE io_mod
