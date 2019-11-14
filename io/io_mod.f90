!===============================================================================
!> $Id$
!===============================================================================
MODULE io_mod
  USE omp_mod
  USE omp_timer_mod
  USE mpi_mod
  USE io_unit_mod
  USE timer_mod
  USE os_mod
  implicit none
  private
  logical:: mpi_trace=.false.
  type io_t
    integer:: verbose=0, input, output, data_unit
    logical:: master=.true.
    logical:: do_legacy, do_direct, omp_trace
    character(len=16):: method='legacy'
    integer:: levelmax
    integer:: word_size=4
    integer:: iodir=-1
    integer:: format=0
    integer:: nml_version=0
    integer:: restart=-1
    integer:: nwrite=0, ntask=0, ntotal=0
    integer:: dims(3)=0, mpi_dims(3), mpi_odims(3)=1
    character(len=64):: inputname, outputname, datadir='data', top
    character(len=64):: rundir='data', inputdir='data'
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
    logical:: namelist_errors=.true.
    logical:: needs_check=.false.
    logical:: halt=.false.      ! used in data_io and amr_io
    integer:: task_logging=0
    integer:: log_sent=0
    integer:: time_derivs=0
    integer:: id_debug=0
    integer:: id_track=0
    integer:: nv=0
    integer:: id=1              ! the task to print a one-liner for
    real(8):: end_time=1d30
    real(8):: out_time=1d0
    real(8):: out_next=1d0
    real(8):: print_time=0d0
    real(8):: print_next=0d0
    real(8):: dtime=huge(1d0)   ! the shortest time step
    real(8):: gb=0.0
    real(8):: job_seconds=1d30
    real(8):: processing=0d0
    real:: flag_time=1.
    real:: flag_max=10.
    real:: ampl=0.0
    real:: grace=0.0
    real:: smallr=0.0
    real:: courant=0.0
    real:: dmin=1.0
    real:: dmax=1.0
    real:: gb_out=0.0
  contains
    procedure:: init
    procedure:: debug
    procedure:: check_flags
    procedure:: bits_mem
    procedure:: gb_mem
    procedure:: gb_print
    procedure:: namelist_warning
    procedure, nopass:: abort
    procedure, nopass:: assert
    procedure, nopass:: header
    procedure, nopass:: print_hl
  end type
  character(len=64), save:: inputname, outputname
  type(io_t), public:: io
  public mch, io_unit, stderr, stdout, stdin
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
  character(len=120), save:: id = &
  'io_mod.f90 $Id$'
  integer, save:: verbose=0, levelmax=10
  character(len=*), optional:: name
  character(len=64):: filename, datadir='data', inputdir='data'
  character(len=64), save:: top='../../'
  character(len=8), save:: method='legacy'
  integer, save:: id_debug=-1, restart=-9, format=0, time_derivs=0, nml_version=1
  integer, save:: log_sent=0, task_logging=0
  logical, save:: first_time=.true., omp_trace=.false., do_validate=.false.
  logical, save:: do_debug=.false., do_trace=.false., do_output=.false., exist, &
    do_legacy=.false., do_flags=.false., do_direct=.false., guard_zones=.false.
  logical ,save:: namelist_errors=.true.
  namelist /io_params/ verbose, do_debug, do_trace, do_output, do_flags, &
    do_validate, do_legacy, do_direct, levelmax, omp_trace, id_debug, top, &
    datadir, inputdir, method, restart, format, guard_zones, time_derivs, log_sent, &
    task_logging, namelist_errors, nml_version
  real:: test
  integer:: iostat
  character(len=120):: ids = &
  '$Id$ io/io_mod.f90'
  !.............................................................................
  if (mpi%master) then
    print'(a)', self%hl
    print'(a)', trim(ids)
  end if
  self%master = mpi%master
  io%master   = mpi%master
  !$omp parallel
  !$omp end parallel
  !-----------------------------------------------------------------------------
  ! I/O unit numbers, for convenience and backwards compatibility.  All other
  ! unit numbers should be taken directly from io_unit_mod
  !-----------------------------------------------------------------------------
  self%input     = io_unit%input
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
      self%rundir = trim(datadir)//self%sep// &                 ! data/run_name/
        trim(filename(1:index(filename,'.')-1))//self%sep
    else
      self%inputname = 'input.nml'
      self%rundir = trim(datadir)//self%sep                     ! data/
    end if
    call os%mkdir (trim(self%rundir))
    self%outputname = self%rundir                               ! compatibility
    write (filename,'(a,i5.5,"/")') trim(self%outputname), 0
    call os%mkdir (trim(filename))
    !---------------------------------------------------------------------------
    open (io_unit%nml, file=trim(self%rundir)//'params.nml', form='formatted', status='unknown')
    open (io_unit%hash, file=trim(self%rundir)//'hash.log', form='formatted', status='unknown')
    open (self%input, file=self%inputname, form='formatted', status='old')
    !---------------------------------------------------------------------------
    ! Read io_params namelist
    !---------------------------------------------------------------------------
    inputdir = self%rundir
    rewind (io_unit%input)
    read (io_unit%input, io_params, iostat=iostat)
    if (iostat > 0) call self%namelist_warning ('io_params')
    if (mpi%master .and. .not. do_validate) then
      print'(a,i4)', ' n_socket =', omp%nsockets
      print'(a,i4)', '   n_core =', omp%ncores
      print'(a,i4)', ' n_thread =', omp_nthreads
    end if
    self%task_logging = task_logging
    self%format = format
    self%nml_version = nml_version
    self%inputdir = inputdir
    call ensure_dirname (self%inputdir)
    call ensure_dirname (self%outputname)
    !---------------------------------------------------------------------------
    ! OMP specific setup; need to set io_unit in all threadprivate instances
    !---------------------------------------------------------------------------
    !$omp parallel
    !$omp critical (open_cr)
    if (omp%master) then
      io_unit%verbose = verbose
    else
      io_unit%verbose = -2
    end if
    io_unit%rundir = self%rundir
    io_unit%inputdir = self%inputdir    
    io_unit%outputname = self%outputname
    io_unit%do_validate = do_validate
    !---------------------------------------------------------------------------
    ! Log-file names:
    ! -- io_unit%output     : stdout on rank 0, same as io_unit%mpi on rank 1-n
    ! -- io_unit%mpi        : data/run/rank_rrrrr.log
    ! -- io_unit%queue      : data/run/rank_rrrrr.nq
    ! -- io_unit%dispatcher : data/run/rank_rrrrr.disp
    ! -- io_unit%log        : data/run/thread_rrrrr_ttt.log
    !---------------------------------------------------------------------------
    call open_rank_file (io_unit%queue     , ".nq"   , 'formatted')
    call open_rank_file (io_unit%mpi       , ".log"  , 'formatted')
    call open_rank_file (io_unit%task      , ".task" , 'formatted')
    call open_rank_file (io_unit%dbg       , ".dbg"  , 'unformatted')
    call open_rank_file (io_unit%dump      , ".dump" , 'unformatted')
    call open_rank_file (io_unit%validate  , ".val"  , 'unformatted')
    call open_rank_file (io_unit%dispatcher, ".disp" , 'formatted')
    !---------------------------------------------------------------------------
    ! Open MPI rank and OMP thread-specific log files if omp_trace is true,
    ! else fall back to io_unit%log being just MPI rank-specific
    !---------------------------------------------------------------------------
    if (omp_trace .and. omp%nthreads>1) then
      io_unit%log = 110 + omp%mythread()
      call open_rank_file (io_unit%log     , ".log"  , 'formatted', omp%mythread())
    else
      io_unit%log = io_unit%mpi
    end if
    !---------------------------------------------------------------------------
    ! For non-master MPI ranks, redirect log file output to rundir/rank_rrrrr.log
    !---------------------------------------------------------------------------
    if (mpi%rank > 0) then
      io_unit%output = io_unit%mpi
    end if
    !$omp end critical (open_cr)
    io_unit%master = io%master .and. omp%master
    write (io_unit%log,'(a,2i4,3l4)') &
      'io_mod::init io_unit%log, omp%thread, io_unit%master, io_unit%verbose:', &
      io_unit%log, omp%thread, io_unit%master, io_unit%verbose
    !$omp end parallel
    !---------------------------------------------------------------------------
    ! stdout is the same as the the %output unit, and is not threadprivate
    !---------------------------------------------------------------------------
    stdout = io_unit%output
    write (stdout,io_params)
    !---------------------------------------------------------------------------
    ! Backwards compatibility
    !---------------------------------------------------------------------------
    io%output    = io_unit%output
    io%data_unit = io_unit%data
    write(io%output,*) 'logfile:', filename
    write(io%output,*) '======================================================================='
    write(io%output,*) 'NOTE: Reading parameters from '//trim(self%inputname)
    write(io%output,*) ' This version was compiled with default real KIND=', kind(test)
    write(io%output,*) '======================================================================='
    first_time=.false.
  end if
  self%guard_zones = guard_zones
  self%time_derivs = time_derivs
  self%restart   = restart
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
  self%log_sent  = log_sent
  self%levelmax  = levelmax
  self%inputname = inputname
  self%top       = top
  self%method    = method
  self%namelist_errors = namelist_errors
  if (do_legacy) self%method = 'legacy'
  if (do_direct) self%method = 'direct'
  io_unit%top    = top
  io_unit%do_validate = do_validate
  stdout = io_unit%output
  call timer%init
contains
  !-----------------------------------------------------------------------------
  subroutine ensure_dirname (s)
    character(len=*):: s
    integer:: l
    l = len(trim(s))
    if (s(l:l) /= '/') s(l+1:l+1)='/'
  end subroutine
  subroutine open_rank_file (unit, ext, form, thread)
    integer:: unit
    integer, optional:: thread
    character(len=*):: ext, form
    character(len=64):: filename
    !-----------------------------------------------------------------------------
    if (present(thread)) then
      open (unit=unit, file=trim(filename), form=form, status='unknown')
      write (io_unit%threadbase,'(a,"thread_",i5.5,"_",i3.3)') &
        trim(self%outputname), mpi%rank, thread
      filename = trim(io_unit%threadbase)//ext
    else
      write (io_unit%rankbase,'(a,"rank_",i5.5)') &
        trim(self%outputname), mpi%rank
      filename = trim(io_unit%rankbase)//ext
    end if
    open (unit=unit, file=trim(filename), form=form, status='unknown')
  end subroutine
END SUBROUTINE init

!===============================================================================
!> Register increase of memory, counted in storage_size bits, times a count
!===============================================================================
SUBROUTINE bits_mem (self, bits, count, label)
  class (io_t)      :: self
  integer:: bits, count
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (io%verbose>2) then
    if (present(label)) then
      write(io%output,1) bits,' bits per word', count,' words, for', &
        (bits/(8.*1024.**3))*count, ' GB '//trim(label)
      1 format(i3,a,i8,a,f6.3,a)
    else
      write(io%output,1) bits,' bits per word', count,' words, for', &
        (bits/(8.*1024.**3))*count, ' GB'
    end if
  end if
  call self%gb_mem ((bits/(8.*1024.**3))*count)
END SUBROUTINE bits_mem

!===============================================================================
!> Register increase of memory, counted in gigabytes
!===============================================================================
SUBROUTINE gb_mem (self, gb)
  class (io_t):: self
  real:: gb, lgb
  real, save:: gb_next=1.0
  !-----------------------------------------------------------------------------
  if (mpi%master) then
    !$omp atomic capture
    self%gb = self%gb + gb
    lgb = self%gb
    !$omp end atomic
    if (abs(lgb - gb_next) > 1.0) then
      !$omp critical (gb_cr)
      if (abs(lgb - gb_next) > 1.0) then
        do while (abs(lgb - gb_next) > 1.0)
          call self%gb_print(lgb)
          gb_next = gb_next + sign(1.0,lgb - gb_next)
        end do
      end if
      !$omp end critical (gb_cr)
    end if
  end if
END SUBROUTINE gb_mem

!===============================================================================
!> Register increase of memory, counted in gigabytes
!===============================================================================
SUBROUTINE gb_print (self, gb)
  class (io_t):: self
  real :: gb
  !.............................................................................
  print '(1x,a,f8.3,a)', 'process memory allocated:', gb,' GB'
END SUBROUTINE gb_print


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
  file = trim(io%outputname)//'flag'
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
      call io%abort ('stop.flag detected')
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
  f = trim(io%outputname)//trim(file)
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
      if (present(lvalue)) write (io_unit%output,*) 'flag: ', file(1:index(file,'.')-1), ' =', lvalue
      if (present(ivalue)) write (io_unit%output,*) 'flag: ', file(1:index(file,'.')-1), ' =', ivalue
      if (present(rvalue)) write (io_unit%output,*) 'flag: ', file(1:index(file,'.')-1), ' =', rvalue
      if (present(dvalue)) write (io_unit%output,*) 'flag: ', file(1:index(file,'.')-1), ' =', dvalue
    end if
  end if
  end subroutine
END SUBROUTINE check_flags

!===============================================================================
!> Issue message for missing namelist in input file
!===============================================================================
SUBROUTINE namelist_warning (self, namelist, error)
  class(io_t):: self
  character(len=*):: namelist
  logical, optional:: error
  !.............................................................................
  if (io%master) then
    write(stdout,'(a)') ''
    write(stdout,'(a)') '*************************************************************************************'
    write(stdout,'(a)') '*************************************************************************************'
    if (self%namelist_errors) then
      write(stdout,'(a)') '           ERROR: namelist '//trim(namelist)//' had read error'
    else
      write(stdout,'(a)') '           WARNING: namelist '//trim(namelist)//' had read error'
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
!> Interface to mpi%abort, for modules that do not otherwise USE mpi_mod
!===============================================================================
SUBROUTINE abort (error)
  character(len=*), optional:: error
  !.............................................................................
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
    call mpi%abort (message)
END SUBROUTINE assert

!===============================================================================
!> Print a centered message, preceded and followed by full === lines, unless the
!> first character is a '-', in which case print a singe ---- message ---- line
!===============================================================================
SUBROUTINE header (str, left)
  character(len=*):: str
  integer, optional:: left
  !.............................................................................
  call timer%print()
  if (io%master) then
    if (str(1:1) == '-') then
      call line (str, left)
    else
      call line ('=', left)
      call line (str, left)
      call line ('=', left)
    end if
  end if
!===============================================================================
!> Print a centered or left adjusted message, surrounded by repeated c chars
!===============================================================================
contains
  subroutine line (s, left)
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
    print '(a)', buf
  end subroutine line
END SUBROUTINE header

!===============================================================================
SUBROUTINE print_hl
  character(len=120), save:: hl= &
    '--------------------------------------------------------------------------------'
  !..............................................................................
  if (io%master) then
    write (io_unit%output,'(a)') trim(hl)
  end if
END SUBROUTINE print_hl

END MODULE io_mod
