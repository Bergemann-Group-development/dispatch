!===============================================================================
!> The run's record of itself: `<run_dir>/job.nml`, written by the MPI master
!> as soon as the run directory exists and rewritten on a normal exit.
!>
!> The record makes a run discoverable and attributable whoever started it (a
!> shell, a notebook, the workflow daemon, a Slurm batch script): the
!> process id and host say whether it still runs, the input file and working
!> directory say what it computes, the Slurm job id lets a scheduler-side tool
!> poll or cancel it. A crash leaves `state='running'` behind; a reader then
!> resolves the state from the process id and the completed outputs (abund.txt).
!> The Python side (`tumult.synthesis.jobrecord`) reads this file and keeps the
!> launcher's own knowledge (queued, cancelled, exit status) beside it.
!===============================================================================
MODULE job_record_mod
  USE iso_c_binding
  USE io_mod
  USE mpi_mod
  implicit none
  private
  !===============================================================================
  !> Writer of the run's job.nml.
  !===============================================================================
  type, public:: job_record_t
  contains
    procedure:: start
    procedure:: finish
  end type
  type(job_record_t), public:: job_record
  character(len=*), parameter:: file_name = 'job.nml' !< name of the record in the run directory
  !=============================================================================
  !> getpid, gethostname and getcwd of the C library.
  !=============================================================================
  interface
    function c_getpid() bind(C, name='getpid')
      import c_int
      integer(c_int):: c_getpid
    end function
    function c_gethostname(name, len) bind(C, name='gethostname')
      import c_char, c_int, c_size_t
      character(kind=c_char):: name(*)
      integer(c_size_t), value:: len
      integer(c_int):: c_gethostname
    end function
    function c_getcwd(buf, size) bind(C, name='getcwd')
      import c_char, c_ptr, c_size_t
      character(kind=c_char):: buf(*)
      integer(c_size_t), value:: size
      type(c_ptr):: c_getcwd
    end function
  end interface
CONTAINS

!===============================================================================
!> Write the record with state='running'. Called once, after io%init.
!===============================================================================
SUBROUTINE start (self)
  class(job_record_t):: self
  !-----------------------------------------------------------------------------
  if (.not. mpi%master) return
  call write_record ('running', '')
END SUBROUTINE start

!===============================================================================
!> Rewrite the record with state='finished'. Only a normal exit gets here.
!===============================================================================
SUBROUTINE finish (self)
  class(job_record_t):: self
  !-----------------------------------------------------------------------------
  if (.not. mpi%master) return
  call write_record ('finished', timestamp())
END SUBROUTINE finish

!===============================================================================
!> The record itself. `started` is kept from the first write, so a finish does
!> not lose it; every other field is what this process knows now.
!===============================================================================
SUBROUTINE write_record (new_state, ended_now)
  !$ use omp_lib
  character(len=*), intent(in):: new_state !< 'running' or 'finished'
  character(len=*), intent(in):: ended_now !< end time stamp, '' while running
  character(len=16):: state
  character(len=32):: started, ended
  character(len=256):: host
  character(len=4096):: cwd, input
  character(len=32):: slurm_job_id
  character(len=16):: dispatch_hash
  integer:: pid, mpi_size, omp_threads
  integer:: lunit, iostat
  logical:: exist
  namelist /job_record/ state, pid, host, cwd, input, started, ended, &
    mpi_size, omp_threads, slurm_job_id, dispatch_hash
  !-----------------------------------------------------------------------------
  started = ''
  inquire (file=trim(io%run_dir)//file_name, exist=exist)
  if (exist) then
    ! keep the recorded start; a rewrite only advances the state
    open (newunit=lunit, file=trim(io%run_dir)//file_name, form='formatted', &
          status='old', action='read')
    pid = 0; host = ''; cwd = ''; input = ''; ended = ''
    mpi_size = 0; omp_threads = 0; slurm_job_id = ''; dispatch_hash = ''
    read (lunit, job_record, iostat=iostat)
    close (lunit)
    if (iostat /= 0) started = ''
  end if
  if (started == '') started = timestamp()
  state = new_state
  ended = ended_now
  pid = int(c_getpid())
  host = c_string(gethostname_())
  cwd = c_string(getcwd_())
  input = io%inputname
  mpi_size = mpi%size
  omp_threads = 1
  !$ omp_threads = omp_get_max_threads()
  call get_environment_variable ('SLURM_JOB_ID', slurm_job_id, status=iostat)
  if (iostat /= 0) slurm_job_id = ''
  dispatch_hash = ''
#ifdef HASH
  dispatch_hash = HASH
#endif
  open (newunit=lunit, file=trim(io%run_dir)//file_name, form='formatted', &
        status='replace', action='write', iostat=iostat)
  if (iostat /= 0) then
    ! the record is an observability convenience; it must not stop the run
    write (io%output,*) 'job_record: could not write '//trim(io%run_dir)//file_name
    return
  end if
  write (lunit, job_record)
  close (lunit)
END SUBROUTINE write_record

!===============================================================================
!> The local time in ISO 8601 with the UTC offset, e.g. 2026-11-03T12:00:00+01:00.
!===============================================================================
FUNCTION timestamp () result (text)
  character(len=32):: text
  character(len=8):: date
  character(len=10):: time
  character(len=5):: zone
  !-----------------------------------------------------------------------------
  call date_and_time (date=date, time=time, zone=zone)
  text = date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'T'// &
         time(1:2)//':'//time(3:4)//':'//time(5:6)//zone(1:3)//':'//zone(4:5)
END FUNCTION timestamp

!===============================================================================
!> The host name as a NUL-terminated buffer, empty if gethostname fails.
!===============================================================================
FUNCTION gethostname_ () result (buf)
  character(kind=c_char, len=256):: buf
  integer(c_int):: status
  !-----------------------------------------------------------------------------
  buf = repeat(c_null_char, len(buf))
  status = c_gethostname (buf, int(len(buf), c_size_t))
  if (status /= 0) buf = c_null_char
END FUNCTION gethostname_

!===============================================================================
!> The working directory as a NUL-terminated buffer, empty if getcwd fails.
!===============================================================================
FUNCTION getcwd_ () result (buf)
  character(kind=c_char, len=4096):: buf
  type(c_ptr):: ptr
  !-----------------------------------------------------------------------------
  buf = repeat(c_null_char, len(buf))
  ptr = c_getcwd (buf, int(len(buf), c_size_t))
  if (.not. c_associated(ptr)) buf = c_null_char
END FUNCTION getcwd_

!===============================================================================
!> A NUL-terminated C buffer as a Fortran string.
!===============================================================================
FUNCTION c_string (buf) result (text)
  character(kind=c_char, len=*), intent(in):: buf !< C string, up to its first NUL
  character(len=len(buf)):: text
  integer:: n
  !-----------------------------------------------------------------------------
  n = index(buf, c_null_char) - 1
  if (n < 0) n = len(buf)
  text = buf(1:n)
END FUNCTION c_string

END MODULE job_record_mod
