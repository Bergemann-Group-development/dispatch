!> =============================================================================
!> $Id$
!> =============================================================================
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
    integer:: verbose=0, input, output, data_unit, log_unit
    logical:: master=.true.
    logical:: do_legacy, do_direct, omp_trace
    character(len=16):: method='legacy'
    integer:: levelmax
    integer:: word_size=4
    integer:: iodir=-1
    integer:: format=0
    integer:: restart=-1
    integer:: nwrite=0, ntask=0, ntotal=0
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
    logical:: time_derivs=.false.
    integer:: id_debug=0
    integer:: id_track=0
    integer:: nv=0
    real(8):: end_time=1d30
    real(8):: out_time=1d0
    real(8):: print_time=0d0
    real:: flag_time=1.
    real:: flag_max=10.
    real:: ampl=0.0
    real(8):: gb=0.0
    real(8):: job_seconds=1d30
  contains
    procedure:: init
    procedure:: debug
    procedure:: check_flags
    procedure:: bits_mem
    procedure:: gb_mem
    procedure, nopass:: namelist_warning
    procedure, nopass:: abort
    procedure, nopass:: assert
    procedure, nopass:: header
  end type
  character(len=64), save:: inputname, outputname
  type(io_t), public:: io
  public mch, io_unit, log_unit
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
  character(len=64):: filename, datadir='data', inputdir='data'
  character(len=64), save:: top='../../'
  character(len=8), save:: method='legacy'
  integer, save:: id_debug=-1, restart=-9
  logical, save:: first_time=.true., omp_trace=.false., do_validate=.false.
  logical, save:: do_debug=.false., do_trace=.false., do_output=.false., exist, &
    do_legacy=.false., do_flags=.false., do_direct=.false.
  namelist /io_params/ verbose, do_debug, do_trace, do_output, do_flags, &
    do_validate, do_legacy, do_direct, levelmax, omp_trace, id_debug, top, &
    datadir, inputdir, method, restart
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
      self%rundir = trim(datadir)//self%sep// &                 ! data/run_name/
        trim(filename(1:index(filename,'.')-1))//self%sep
    else
      self%inputname = 'input.nml'
      self%rundir = trim(datadir)//self%sep                     ! data/
    end if
    self%outputname = self%rundir                               ! compatibility
    call os%mkdir (trim(self%rundir))
    open (io_unit%nml, file=trim(self%rundir)//'params.nml', form='formatted', status='unknown')
    open (self%input, file=self%inputname, form='formatted', status='old')
    !---------------------------------------------------------------------------
    ! Read io_params namelist
    !---------------------------------------------------------------------------
    inputdir = self%rundir
    rewind (io_unit%input)
    read (io_unit%input, io_params, iostat=iostat)
    if (self%master.and.first_time) write (*,io_params)
    if (mpi%master .and. .not. do_validate) then
      print'(a,i4)', ' n_socket =', omp%nsockets
      print'(a,i4)', '   n_core =', omp%ncores
      print'(a,i4)', ' n_thread =', omp_nthreads
    end if
    self%inputdir = inputdir
    !---------------------------------------------------------------------------
    ! OMP specific setup
    !---------------------------------------------------------------------------
    !$omp parallel
    !$omp critical (open_cr)
    io_unit%rundir = self%rundir
    io_unit%inputdir = self%inputdir    
    io_unit%outputname = self%outputname
    if (omp%master) then
      io_unit%verbose = verbose
    else
      io_unit%verbose = -2
    end if
    !---------------------------------------------------------------------------
    ! Open MPI rank and OMP thread-specific log files
    !---------------------------------------------------------------------------
    io_unit%do_validate = do_validate
    write (filename,'(a,i5.5,"/")') trim(io%outputname), 0
    call os%mkdir (trim(filename))
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
  self%levelmax  = levelmax
  self%inputname = inputname
  self%top       = top
  self%method    = method
  if (do_legacy) self%method = 'legacy'
  if (do_direct) self%method = 'direct'
  io_unit%top    = top
  io_unit%do_validate = do_validate
  call timer%init
END SUBROUTINE init

!> =============================================================================
!> Register increase of memory, counted in storage_size bits, times a count
!> =============================================================================
SUBROUTINE bits_mem (self, bits, count, label)
  class (io_t)      :: self
  integer:: bits, count
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (io%verbose>2) then
    if (present(label)) then
      print *, bits,' bits per word', count,' words, for', &
        (bits/(8.*1024.**3))*count, 'GB '//trim(label)
    else
      print *, bits,' bits per word', count,' words, for', &
        (bits/(8.*1024.**3))*count, 'GB'
    end if
  end if
  call self%gb_mem ((bits/(8.*1024.**3))*count)
END SUBROUTINE bits_mem

!> =============================================================================
!> Register increase of memory, counted in gigabytes
!> =============================================================================
SUBROUTINE gb_mem (self, gb)
  class (io_t)      :: self
  real:: gb
  real, save:: gb_next=1.0
  !-----------------------------------------------------------------------------
  if (mpi%master) then
    self%gb = self%gb + gb
    if (abs(self%gb - gb_next) > 1.0) then
      !$omp critical (gb_cr)
      do while (abs(self%gb - gb_next) > 1.0)
        print '(1x,a,f8.3,a)', 'process memory allocated:', self%gb,' GB'
        gb_next = gb_next + sign(1d0,self%gb - gb_next)
      end do
      !$omp end critical (gb_cr)
    end if
  end if
END SUBROUTINE gb_mem

!> =============================================================================
!> Check debug flag files.  To avoid having to synchronize when using MPI, we
!> give the ranks 5 seconds to detect and read flag files, and then we let one
!> rank remove them.
!> =============================================================================
SUBROUTINE check_flags (self)
  class (io_t)      :: self
  real(8), save     :: last_checked=0d0
  real(8)           :: flag_time, processing=0d0
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
  if (processing==0d0) then
    if (wallclock() > last_checked+flag_time) then
      last_checked = wallclock()
      call check_all (.false.)
    end if
  !-----------------------------------------------------------------------------
  ! After processing for 10 sec, all ranks turn off processing
  !-----------------------------------------------------------------------------
  else if (wallclock() > processing+10d0) then
    processing = 0d0
  !-----------------------------------------------------------------------------
  ! After 5 sec, the master thread on the master rank removes the flag files
  !-----------------------------------------------------------------------------
  else if (wallclock() > processing+5d0) then
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
    processing = last_checked + flag_time
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
    ! double prec
    call check ('out_time.flag'       , remove, dvalue=self%out_time)
    call check ('print_time.flag'     , remove, dvalue=self%print_time)
    call check ('end_time.flag'       , remove, dvalue=self%end_time)
    call check ('sec_per_report.flag' , remove, dvalue=timer%sec_per_report)
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
      if (present(lvalue)) write (*,*) 'flag: ', file(1:index(file,'.')-1), ' =', lvalue
      if (present(ivalue)) write (*,*) 'flag: ', file(1:index(file,'.')-1), ' =', ivalue
      if (present(rvalue)) write (*,*) 'flag: ', file(1:index(file,'.')-1), ' =', rvalue
      if (present(dvalue)) write (*,*) 'flag: ', file(1:index(file,'.')-1), ' =', dvalue
    end if
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
    print *, ''
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

!> =============================================================================
!> Interface to mpi%abort, for modules that do not otherwise USE mpi_mod
!> =============================================================================
SUBROUTINE assert (ok, message)
  logical:: ok
  character(len=*):: message
  !.............................................................................
  if (.not.ok) &
    call mpi%abort (message)
END SUBROUTINE assert

!> =============================================================================
!> Interface to mpi%abort, for modules that do not otherwise USE mpi_mod
!> =============================================================================
SUBROUTINE header (message)
  character(len=*):: message
  !.............................................................................
  call line('')
  call line(message)
  call line('')
!> =============================================================================
contains
  subroutine line (s)
    character(len=*):: s
    integer, parameter:: w=80
    character(len=w):: buf
    integer:: i, j, l1, l2, l
    !-----------------------------------------------------------------------------
    l = len_trim(s)
    l1 = (w-l)/2-1
    l2 = l1 + l + 1
    j = 1
    do i=1,w
      buf(i:i) = '='
      if (l==0) then
        buf(i:i) = '='
      else if (i==l1 .or. i==l2) then
        buf(i:i) = ' '
      else if (i > l1 .and. i < l2 .and. j <= l) then
        buf(i:i) = s(j:j)
        j = j+1
      end if
    end do
    print '(a)', buf
  end subroutine
END SUBROUTINE header

END MODULE io_mod
