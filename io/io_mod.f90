!===============================================================================
!> $Id$
!===============================================================================
MODULE io_mod
  USE iso_fortran_env
  USE omp_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE mpi_mod
  USE io_unit_mod
  USE timer_mod
  USE os_mod
  USE trace_mod
  implicit none
  private
  logical:: mpi_trace=.false.
  type input_t
    integer:: mpi_dims(3), mpi_size, format, ntotal, time_derivs, nml_version
    integer:: nv, mv
    real(8):: out_time
    logical:: guard_zones, no_mans_land, do_generic
    character(16):: method
  end type
  type io_t
    integer:: verbose=0, input, output=6, data_unit, namelist
    logical:: master=.true.
    logical:: do_legacy, do_direct, omp_trace
    character(len=16):: method
    integer:: levelmax
    integer:: word_size=4
    integer:: iodir=-1
    integer:: format=0
    integer:: nml_version=4
    integer:: restart=-1
    integer:: ntask=0           ! number of tasks on this rank (?)
    integer:: nwrite=0          ! number of tasks to output on this rank
    integer:: ntotal=0          ! number of tasks to output on all ranks
    !integer:: dims(3)=0         ! FIXME
    integer:: mpi_dims(3)=1     ! FIXME
    integer:: print_every=0, iprint=0       ! print cadence
    character(len=512):: inputname          ! input namelist file (e.g. run.nml)
    character(len=512):: inputdir='data/'   ! input directory (restart_mod)
    character(len=512):: run_name=''        ! run part of input file
    character(len=512):: data_dir='data/'   ! data
    character(len=512):: run_dir='data/'    ! data/run
    character(len=512):: outputname         ! same as run_dir, for compatibily
    character(len=512):: top                ! base name
    character(len=1):: sep = '/'
    character(len=72):: &
    hl='-------------------------------------------------------------------------'
    character(len=72):: &
    hs='*************************************************************************'
    ! flag values
    logical:: do_trace=.false.
    logical:: do_output=.true.
    logical:: do_debug=.false.
    logical:: do_flags=.false.
    logical:: do_stop=.false.
    logical:: do_generic=.false.
    logical:: guard_zones=.false.
    logical:: namelist_errors=.true.
    logical:: needs_check=.false.
    logical:: halt=.false.      ! used in data_io and amr_io
    integer:: task_logging=0
    integer:: log_sent=0
    integer:: time_derivs=0
    integer:: id_debug=0, key_debug(4)=0
    integer:: id_track=0
    integer:: nv=0
    integer:: id=1              ! the task to print a one-liner for
    integer:: max_files=1000
    real(8):: start_time=0d0
    real(8):: end_time=1d30
    real(8):: out_time=1d0
    real(8):: out_next=1d0
    real(8):: print_time=0d0
    real(8):: print_seconds=0d0
    real(8):: print_next=0d0
    real(8):: dtime=huge(1d0)   ! the shortest time step
    real(8):: job_seconds=1d30
    real(8):: processing=0d0
    real(8):: log_interval
    real:: flag_time=1.
    real:: flag_max=10.
    real:: ampl=0.0
    real:: grace=0.0
    real:: smallr=0.0
    real:: courant=0.0
    real:: dmin=1.0
    real:: dmax=1.0
    real:: gb_out=0.0
    type(mpi_counter_t):: counter
    type(lock_t):: lock
    type(input_t):: in
  contains
    procedure:: init
    procedure:: debug
    procedure:: check_flags
    procedure:: bits_mem
    procedure:: namelist_warning
    procedure, nopass:: open_rank_file
    procedure, nopass:: warning
    procedure, nopass:: abort
    procedure, nopass:: assert
    procedure, nopass:: header
    procedure, nopass:: line
    procedure, nopass:: print_hl
    procedure, nopass:: ensure_dirname
  end type
  type(io_t), public:: io
  public mch, io_unit, stderr, stdout, stdin, stdlog
CONTAINS

!===============================================================================
!> Function suitable for testing in debug print statements
!===============================================================================
LOGICAL FUNCTION debug (self, verbose)
  class(io_t):: self
  integer:: verbose
  debug = self%master .and. (self%verbose >= verbose)
END FUNCTION debug

!===============================================================================
!> Initialize I/O related parameters.  The first thread that enters here blocks
!> access in a critical region, reads namelist parameters, and stores them in a
!> global, public instance, which we normally refer to.  Other calls pick up the
!> global values and store them in other instances (which may be the same).
!===============================================================================
SUBROUTINE init (self, name)
  class(io_t):: self
  integer, save:: verbose=0, levelmax=10
  character(len=*), optional:: name
  character(len=256):: filename, datadir='data'
  character(len=256), save:: top='../../'
  character(len=16), save:: method='background'
  integer, save:: id_debug=-1, restart=-9, format=0, time_derivs=0
  integer, save:: log_sent=0, task_logging=0, nml_version=4
  integer, save:: print_every=0, max_files=1000, os_verbose=0
  logical, save:: first_time=.true., omp_trace=.false., do_validate=.false.
  logical, save:: do_debug=.false., do_trace=.false., do_output=.false., exist
  logical, save:: do_legacy=.false., do_flags=.false., do_direct=.false.
  logical ,save:: namelist_errors=.true., guard_zones=.false.
  logical ,save:: do_generic=.false., flush_logs=.false.
  real(8):: log_interval=5d0, print_seconds=1d0, gb_step=1d0
  real(8), save :: start_time=0d0, end_time=1.0, out_time=0.1, print_time=0.0
  real:: test
  integer:: i, iostat, i1, i2
  character(len=10), save:: hash=''
  character(len=16), save:: solver=''
  character(len=32), save:: branch='', opts=''
  character(len=128), save:: compiler=''
  character(len=120):: ids = &
  '$Id$'
  namelist /io_params/ verbose, do_debug, do_trace, do_output, do_flags, &
    do_validate, do_legacy, do_direct, do_generic, levelmax, omp_trace, &
    id_debug, method, restart, format, guard_zones, time_derivs, log_sent, &
    task_logging, namelist_errors, nml_version, log_interval, top, datadir, &
    flush_logs, print_seconds, max_files, print_time, print_every, &
    start_time, out_time, end_time, os_verbose, gb_step
  namelist /out_params/ start_time, out_time, end_time, print_every
  namelist /code_params/ branch, hash, compiler, solver
  !-----------------------------------------------------------------------------
  ! I/O unit numbers, for convenience and backwards compatibility.  All other
  ! unit numbers should be taken directly from io_unit_mod
  !-----------------------------------------------------------------------------
  self%input     = io_unit%input
  self%namelist  = io_unit%nml
  !-----------------------------------------------------------------------------
  ! Use serial region to avoid threads using the input file simultaneously
  !-----------------------------------------------------------------------------
  if (first_time) then
    !---------------------------------------------------------------------------
    ! The input filename is input.nml, or as given on the command line
    !---------------------------------------------------------------------------
    call getarg(1,filename)                                     ! command line
    if (filename/=' '.and.trim(filename)/='input.nml') then
      self%inputname = filename                                 ! run.nml
      i = index(filename,'.',back=.true.) - 1                   ! last char
      self%run_name = filename(1:i)                             ! run
    else
      self%inputname = 'input.nml'
      self%run_name = ''
    end if
    
    if (mpi%rank == 0) &
      write(*,'(a,i2)') &
      ' This version was compiled with default real KIND =', kind(test)
    if (mpi%rank > 0) then
      io%output = io_unit%mpi
    end if
    open (self%input, file=self%inputname, form='formatted', status='old')
    !---------------------------------------------------------------------------
    ! Read out_params (obsolete)
    !---------------------------------------------------------------------------
    rewind (io%input)
    read (io%input, out_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, out_params)
    end if
    call io%namelist_warning ('out_params', iostat)
    !---------------------------------------------------------------------------
    ! Read io_params
    !---------------------------------------------------------------------------
    rewind (io%input)
    read (io%input, io_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, io_params)
    end if
    if (print_every > 0) &
      print_seconds = 0d0
    io%start_time  = start_time
    io%end_time    = end_time
    io%out_time    = out_time
    io%print_time  = print_time
    io%print_every = print_every                ! print 1st task
    io%iprint      = print_every                ! print 1st task
    os%gb_step     = gb_step
    call self%namelist_warning ('io_params', iostat)
    !---------------------------------------------------------------------------
    ! `mkdir -p` doesn't work from inside the code, so two `mkdir` are required
    call ensure_dirname (datadir)
    self%run_dir = trim(datadir)//trim(self%run_name)
    call ensure_dirname (self%run_dir)
    call os%mkdir (trim(datadir))
    call os%mkdir (trim(self%run_dir))
    call os%copy_file (trim(self%run_dir)//'params.nml', &
                       trim(self%run_dir)//'restart.nml')
    if (mpi%master) then
      write (io_unit%hash,'(a)') trim(ids)
    end if
    !---------------------------------------------------------------------------
    open (io_unit%nml, file=trim(self%run_dir)//'params.nml', form='formatted', &
          delim='apostrophe', status='unknown')
    open (io_unit%hash, file=trim(self%run_dir)//'hash.log', form='formatted', &
          status='unknown')
    open (self%input, file=self%inputname, form='formatted', status='old')
    !---------------------------------------------------------------------------
    self%task_logging = task_logging
    self%log_interval = log_interval
    self%format = format
    self%nml_version = nml_version
    self%outputname = self%run_dir                               ! compatibility
    self%do_generic = do_generic
    self%print_seconds = print_seconds
    !---------------------------------------------------------------------------
    ! Initialize the modules needed for comprehensive I/O
    !---------------------------------------------------------------------------
    call omp%init ()
    !$omp parallel
    !call omp%info ()
    call omp_lock%init ('io')
    !$omp end parallel
    call os%init (os_verbose)
    call mpi%init ()
    mpi%master = mpi%rank==0
    io%master = mpi%master
    !---------------------------------------------------------------------------
    ! Make the run directory
    !---------------------------------------------------------------------------
    write (filename,'(a,i5.5,"/")') trim(self%run_dir), 0
    if (io%do_output) call os%mkdir (trim(filename))
    !---------------------------------------------------------------------------
    ! OMP specific setup; need to set io_unit in all threadprivate instances
    !---------------------------------------------------------------------------
    !$omp parallel
    !$omp critical (open_cr)
    if (omp%master) then
      io_unit%verbose = verbose
    else if (.not. omp_trace) then
      io_unit%verbose = -2
    end if
    io_unit%do_trace = do_trace
    io_unit%omp_trace = omp_trace
    io_unit%run_dir = self%run_dir
    io_unit%do_validate = do_validate
    io_unit%flush_logs = flush_logs
    !---------------------------------------------------------------------------
    ! Log-file names:
    ! -- io_unit%output     : stdout on rank 0, same as io_unit%mpi on rank 1-n
    ! -- io_unit%mpi        : data/run/rank_rrrrr.log
    ! -- io_unit%queue      : data/run/rank_rrrrr.nq
    ! -- io_unit%dispatcher : data/run/rank_rrrrr.disp
    ! -- io_unit%log        : data/run/thread_rrrrr_ttt.log
    !---------------------------------------------------------------------------
    call open_rank_file (io_unit%queue     , ".que"  , 'formatted')
    call open_rank_file (io_unit%mpi       , ".log"  , 'formatted')
    call open_rank_file (io_unit%dbg       , ".dbg"  , 'unformatted')
    call open_rank_file (io_unit%dump      , ".dump" , 'unformatted')
    call open_rank_file (io_unit%validate  , ".val"  , 'unformatted')
    call open_rank_file (io_unit%dispatcher, ".disp" , 'formatted')
    call open_rank_file (io_unit%particle_d, ".pdat" , 'unformatted')
    call open_rank_file (io_unit%particle_i, ".pinf" , 'unformatted')
    if (self%task_logging > 0) &
    call open_rank_file (io_unit%task      , ".task1", 'formatted')
    !---------------------------------------------------------------------------
    ! Open MPI rank and OMP thread-specific log files if omp_trace is true,
    ! else fall back to io_unit%log being just MPI rank-specific
    !---------------------------------------------------------------------------
    if (mpi%master) then
      io%namelist = io_unit%nml
    else
      io%namelist = io_unit%mpi
    end if
    if (omp_trace .and. omp%nthreads>1) then
      io_unit%log = 0
      call self%open_rank_file (io_unit%log,".log",'formatted',omp%mythread())
    else
      io_unit%log = io_unit%mpi
    end if
    stdlog = io_unit%log
    call omp_lock%init_unit ()
    !---------------------------------------------------------------------------
    ! For non-master MPI ranks, redirect log file output to run_dir/rank_rrrrr.log
    !---------------------------------------------------------------------------
    if (mpi%rank > 0) then
      io_unit%output = io_unit%mpi
    end if
    call io%lock%init ('io')
    !$omp end critical (open_cr)
    io_unit%master = io%master .and. omp%master
    write (io_unit%log,'(a,2i4,l,i4)') &
      'io_mod::init io_unit%log, omp%thread, io_unit%master, io_unit%verbose:', &
      io_unit%log, omp%thread, io_unit%master, io_unit%verbose
    !$omp end parallel
    !---------------------------------------------------------------------------
    ! stdout is the same as the the %output unit, and is not threadprivate
    !---------------------------------------------------------------------------
    stdout = io_unit%output
    io%output = io_unit%output
    write (self%namelist, out_params)
    write (self%namelist, io_params)
    !---------------------------------------------------------------------------
    ! Backwards compatibility
    !---------------------------------------------------------------------------
    io%data_unit = io_unit%data
    first_time=.false.
  end if
  !-----------------------------------------------------------------------------
  ! Code information
  !-----------------------------------------------------------------------------
#ifdef BRANCH
  if (mpi%rank == 0) then
    print '(a)', &
' ==============================================================================='
    write (stdout,*)  'DISPATCH code framework, branch:', BRANCH, ', hash:', HASH
    write (stdout,*)  'SOLVER=', SOLVER, ' OPTS=', OPTS
    solver = SOLVER
    branch = BRANCH
    hash = HASH
    opts = OPTS
#ifdef __NVCOMPILER
#if __NVCOMPILER_MAJOR__ < 22
    compiler = 'nvfortran'
#endif
#else
    compiler = compiler_version()
#endif
    i1 = index(compiler,'(R)')+3
    i2 = index(compiler,'Version')
    if (i1>0 .and. i2>0) then
      write (stdout,*) compiler(1:i1)//compiler(i2:)
    else
      write (stdout,*) compiler
    end if
    write (io_unit%nml, code_params)
  end if
#endif
  !-----------------------------------------------------------------------------
  ! MPI/OMP info
  !-----------------------------------------------------------------------------
  call io%header ('mpi_t%info:')
  call mpi%info
  if (mpi%master) then
    if (.not. do_validate) then
      print'(a,i4)', ' n_sockets =', omp%nsockets
      print'(a,i4)', '   n_cores =', omp%ncores
      print'(a,i4)', ' n_threads =', omp%nthreads
      print'(a,i4)', '  n_places =', omp%nplaces
      print'(a,i4)', '   n_ranks =', mpi%size
      print'(a,i6,a2)', ' Stacksize =', int(omp%stacksize/1024.**2), 'MB'
    end if
    write(*,'(a,i2)') &
    ' This version was compiled with default real KIND =', kind(test)
    call io%header ('io_t%init:')
    write (io%output,*) 'parameters from: '//trim(self%inputname)
    write (io%output,*) 'output directory: '//trim(self%run_dir)
  end if
  self%guard_zones = guard_zones
  self%time_derivs = time_derivs
  self%restart   = restart
  self%data_dir  = datadir
  self%id_debug  = id_debug
  self%verbose   = verbose
  self%do_debug  = do_debug
  self%do_trace  = do_trace
  self%do_output = do_output
  self%do_flags  = do_flags
  self%do_legacy = do_legacy
  self%do_direct = do_direct
  self%omp_trace = omp_trace
  self%log_sent  = log_sent
  self%levelmax  = levelmax
  self%top       = top
  self%method    = method
  self%namelist_errors = namelist_errors
  if (do_legacy) self%method = 'legacy'
  if (do_direct) self%method = 'direct'
  io_unit%top    = top
  io_unit%do_validate = do_validate
  stdout = io_unit%output
  io%output = io_unit%output
  call io%counter%init (mpi%size)
  call omp_lock%init ('io')
  call timer%init
END SUBROUTINE init

!===============================================================================
!> Make sure directory name ends with exactly one slash
!===============================================================================
SUBROUTINE ensure_dirname (s)
  character(len=*):: s
  integer:: l
  l = len(trim(s))
  if (s(l:l) == '/') then
    if (s(l-1:l-1) == '/') s(l:l) = ' '
  else
    if (s(l:l) /= '/') s(l+1:l+1)='/'
  end if
END SUBROUTINE

!===============================================================================
!> Generic text file open
!===============================================================================
SUBROUTINE open_rank_file (unit, ext, form, thread)
  integer:: unit
  integer, optional:: thread
  character(len=*):: ext, form
  character(len=256):: filename
  !-----------------------------------------------------------------------------
  if (present(thread)) then
    write (io_unit%threadbase,'(a,"thread_",i5.5,"_",i3.3)') &
      trim(io%run_dir), mpi%rank, thread
    filename = trim(io_unit%threadbase)//ext
  else
    write (io_unit%rankbase,'(a,"rank_",i5.5)') &
      trim(io%run_dir), mpi%rank
    filename = trim(io_unit%rankbase)//ext
  end if
  if (unit==0) then
    open (newunit=unit, file=trim(filename), form=form, status='unknown')
  else
    open (unit=unit, file=trim(filename), form=form, status='unknown')
  end if
END SUBROUTINE open_rank_file

!===============================================================================
!> Register increase of memory, counted in storage_size bits, times a count
!===============================================================================
SUBROUTINE bits_mem (self, bits, count, label)
  class (io_t)      :: self
  integer:: bits, count
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  call os%mem_register (bits, count, label)
END SUBROUTINE bits_mem

!===============================================================================
!> Check debug flag files.  To avoid having to synchronize when using MPI, we
!> give the ranks 5 seconds to detect and read flag files, and then we let one
!> rank remove them.
!===============================================================================
SUBROUTINE check_flags (self)
  class (io_t)      :: self
  real(8), save     :: last_checked=0d0
  real(8)           :: flag_time
  logical           :: exists
  character(len=80) :: file
  integer, save     :: itimer=0
  !.............................................................................
  !$omp master
  call timer%begin ('io_t%check_flags', itimer)
  flag_time = merge(self%flag_time, self%flag_max, self%do_flags)
  !------------------------------------------------------------------------------
  ! If ready to detect flag files, try.
  !------------------------------------------------------------------------------
  if (io%processing==0d0) then
    if (wallclock() > last_checked+flag_time) then
      last_checked = wallclock()
      call check_all (.false.)
    end if
  !-----------------------------------------------------------------------------
  ! After processing for 10 sec, all ranks turn off processing
  !-----------------------------------------------------------------------------
  else if (wallclock() > io%processing+10d0) then
    io%processing = 0d0
  !-----------------------------------------------------------------------------
  ! After 5 sec, the master thread on the master rank removes the flag files
  !-----------------------------------------------------------------------------
  else if (wallclock() > io%processing+5d0) then
    if (io%master) call check_all (.true.)
  end if
  call timer%end (itimer)
  !$omp end master
contains

  !-----------------------------------------------------------------------------
  ! Thee OMP master checks for all valid types of flag files
  !-----------------------------------------------------------------------------
  subroutine check_all (remove)
  logical:: remove
  !.............................................................................
  !$omp critical (flag_cr)
  file = trim(io%run_dir)//'flag'
  inquire (file=trim(file), exist=exists)
  if (exists) then
    !---------------------------------------------------------------------------
    ! If a flag file is detected, start the processing
    !---------------------------------------------------------------------------
    io%processing = last_checked + flag_time
    open (io_unit%flag, file=trim(file), form='formatted', status='old')
    if (remove) then
      close (io_unit%flag, status='delete')
    else
      close (io_unit%flag)
    end if
    ! logical
    call check ('do_trace.flag'       , remove, lvalue=self%do_trace)
    call check ('do_flags.flag'       , remove, lvalue=self%do_flags)
    call check ('do_debug.flag'       , remove, lvalue=self%do_debug)
    call check ('do_output.flag'      , remove, lvalue=self%do_output)
    call check ('stop.flag'           , remove, lvalue=self%do_stop)
    ! integer
    call check ('id_debug.flag'       , remove, ivalue=self%id_debug)
    call check ('id_track.flag'       , remove, ivalue=self%id_track)
    call check ('verbose.flag'        , remove, ivalue=self%verbose)
    ! single prec
    call check ('flag_time.flag'      , remove, rvalue=self%flag_time)
    call check ('flag_max.flag'       , remove, rvalue=self%flag_max)
    call check ('ampl.flag'           , remove, rvalue=self%ampl)
    call check ('grace.flag   '       , remove, rvalue=self%grace)
    call check ('smallr.flag  '       , remove, rvalue=self%smallr)
    call check ('courant.flag  '      , remove, rvalue=self%courant)
    ! double prec
    call check ('out_time.flag'       , remove, dvalue=self%out_time)
    call check ('out_next.flag'       , remove, dvalue=self%out_next)
    call check ('print_time.flag'     , remove, dvalue=self%print_time)
    call check ('end_time.flag'       , remove, dvalue=self%end_time)
    call check ('sec_per_report.flag' , remove, dvalue=timer%sec_per_report)
    !-------------------------------------------------------------------------
    if (io%do_stop .and. mpi%master) then
      write (stdout,*) file
      open (io_unit%flag, file=file, form='formatted', status='unknown')
      close (io_unit%flag, status='delete')
      call mpi%abort ('stop.flag detected')
    end if
  end if
  !$omp end critical (flag_cr)
  end subroutine

  !=============================================================================
  subroutine check (file, remove, lvalue, ivalue, rvalue, dvalue)
  character(len=*)  :: file
  logical           :: remove
  logical, optional :: lvalue
  integer, optional :: ivalue
  real   , optional :: rvalue
  real(8), optional :: dvalue
  logical           :: exists
  integer           :: iostat
  character(len=80) :: f
  !.............................................................................
  f = trim(io%run_dir)//trim(file)
  if (remove) then
    inquire (file=trim(f), exist=exists)
    if (exists) then
      open (io_unit%flag, file=trim(f), form='formatted', status='old', iostat=iostat)
      close (io_unit%flag, status='delete')
    end if
  else
    inquire (file=trim(f), exist=exists)
    if (exists) then
      open (io_unit%flag, file=trim(f), form='formatted', status='old', iostat=iostat)
      if (present(ivalue)) read (io_unit%flag,*) ivalue
      if (present(rvalue)) read (io_unit%flag,*) rvalue
      if (present(dvalue)) read (io_unit%flag,*) dvalue
      if (present(lvalue).and.iostat==0) read (io_unit%flag,*,iostat=iostat) lvalue
      if (iostat/=0) lvalue = .not.lvalue
      close (io_unit%flag)
      if (present(lvalue)) write (io_unit%output,*) &
        'flag: ', file(1:index(file,'.',back=.true.)-1), ' =', lvalue
      if (present(ivalue)) write (io_unit%output,*) &
        'flag: ', file(1:index(file,'.',back=.true.)-1), ' =', ivalue
      if (present(rvalue)) write (io_unit%output,*) &
        'flag: ', file(1:index(file,'.',back=.true.)-1), ' =', rvalue
      if (present(dvalue)) write (io_unit%output,*) &
        'flag: ', file(1:index(file,'.',back=.true.)-1), ' =', dvalue
    end if
  end if
  end subroutine
END SUBROUTINE check_flags

!===============================================================================
!> Issue message for missing namelist in input file
!===============================================================================
SUBROUTINE namelist_warning (self, namelist, iostat, error)
  class(io_t):: self
  character(len=*):: namelist
  integer, optional:: iostat
  logical, optional:: error
  integer:: iostat_l
  !.............................................................................
  if (io%master) then
    if (present(iostat)) then
      iostat_l = iostat
    else
      iostat_l = 0
    end if
    if (.not.present(error) .and. iostat_l < 0) return
    if (iostat_l == 0) return
    write(stdout,'(a)') ''
    write(stdout,'(a)') '*************************************************************************************'
    write(stdout,'(a)') '*************************************************************************************'
    if (iostat_l < 0) then
      if (self%namelist_errors) then
        write(stdout,'(a)') '           ERROR: namelist '//trim(namelist)//' missing'
      else
        write(stdout,'(a)') '           WARNING: namelist '//trim(namelist)//' missing'
      end if
    else
      if (self%namelist_errors) then
        write(stdout,'(a)') '           ERROR: namelist '//trim(namelist)//' had read error'
      else
        write(stdout,'(a)') '           WARNING: namelist '//trim(namelist)//' had read error'
      end if
    end if
    write(stdout,'(a)') '*************************************************************************************'
    write(stdout,'(a)') '*************************************************************************************'
    if (self%namelist_errors) then
      if (present(error)) then
        if (error) &
          call mpi%abort('Namelist error')
      else
        call mpi%abort('Namelist error')
      end if
    end if
  end if
END SUBROUTINE namelist_warning

!===============================================================================
!> Issue warnings, to stdout if io%verbose > 1, otherwise to rank log
!===============================================================================
SUBROUTINE warning (text)
  character(len=*), optional:: text
  integer:: unit, i, n
  !.............................................................................
  if (io%verbose >= 0) then
    unit = stdout
  else
    unit = io_unit%mpi
  end if
  if (present(text)) then
    n = 80
    write (unit,'(/1x,80a)') ('*',i=1,n)
    write (unit,*) ' io_t%warning: '//trim(text)
    write (unit,'(1x,80a/)') ('*',i=1,n)
    if (unit == stdout) then
      write (stdlog,*) 'io_t%warning: '//trim(text)
      flush (stdlog)
    end if
  else
    write (unit,*) 'io_t%warning: called without specific text'
    if (unit == stdout) then
      write (stdlog,*) 'io_t%warning: called without specific text'
      flush (stdlog)
    end if
  end if
END SUBROUTINE warning

!===============================================================================
!> Interface to mpi%abort, for modules that do not otherwise USE mpi_mod
!===============================================================================
SUBROUTINE abort (error)
  character(len=*), optional:: error
  !.............................................................................¨
  call trace%back ()
  call mpi%abort (error)
END SUBROUTINE abort

!===============================================================================
!> Interface to mpi%abort, for modules that do not otherwise USE mpi_mod
!===============================================================================
SUBROUTINE assert (ok, message)
  logical:: ok
  character(len=*):: message
  !.............................................................................
  if (.not.ok) &
    call io%abort (message)
END SUBROUTINE assert

!===============================================================================
!> Print a centered message, preceded and followed by full === lines, unless the
!> first character is a '-', in which case print a singe ---- message ---- line
!===============================================================================
SUBROUTINE header (str, left)
  character(len=*):: str
  integer, optional:: left
  !.............................................................................
!  call timer%print()
  if (io%master) then
    if (str(1:1) == '-') then
      call line (str, left)
    else
      !call line ('=', left)
      call line (str, left)
    end if
  end if
END SUBROUTINE header

!===============================================================================
!> Print a centered or left adjusted message, surrounded by repeated c chars
!===============================================================================
SUBROUTINE line (s, left)
  character(len=*):: s
  integer, optional:: left
  !...........................................................................
  character(len=1):: c
  integer, parameter:: w=80
  character(len=w):: buf
  integer:: i, j, l1, l2, l
  !---------------------------------------------------------------------------
  l = len_trim(s)
  !---------------------------------------------------------------------------
  ! Optionally, use a fixed left start, else center
  !---------------------------------------------------------------------------
  if (s(1:1)=='-') then
    c = '-'
    j = 2
  else
    c = '='
    j = 1
  end if
  if (present(left)) then
    l1 = left
  else
    l1 = (w-l)/2-1
  end if
  l2 = l1 + l
  do i=1,w
    if (l <= 1) then
      ! --- no message ---
      buf(i:i) = c
    else if (i==l1 .or. i==l2) then
      ! --- space before and after message ---
      buf(i:i) = ' '
    else if (i > l1 .and. i < l2 .and. j <= l) then
      ! --- copy message ---
      buf(i:i) = s(j:j)
      j = j+1
    else
      ! --- surrounding char ---
      buf(i:i) = c
    end if
  end do
  write (stdout,'(1x,a)') buf
END SUBROUTINE line

!===============================================================================
SUBROUTINE print_hl
  character(len=120), save:: hl= &
    '--------------------------------------------------------------------------------'
  !..............................................................................
  if (io%master) then
    write (stdout,'(a)') trim(hl)
  end if
END SUBROUTINE print_hl

END MODULE io_mod
