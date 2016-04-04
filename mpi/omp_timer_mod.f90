!> *****************************************************************************
!> Support tic/toc timing, as in MATLAB, and accurate wallclock() function.
!> The timing is generally much more accurate if this module is compiled with
!> OMP active.
!> *****************************************************************************
MODULE omp_timer_mod
  USE io_unit_mod
  USE omp_mod
  USE mpi_mod
  implicit none
  private
  integer, save:: n_core = 0
  real(8), save:: wt=0d0
  interface toc
    module procedure toc_i4, toc_i8, toc_1
  end interface
  public tic, toc, wallclock
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE tic
  implicit none
  !.............................................................................
  real(8):: omp_get_wtime
  wt = omp_get_wtime()
END SUBROUTINE tic

!> =============================================================================
!> =============================================================================
SUBROUTINE toc_1 (label)
  implicit none
  character(len=*):: label
  integer:: n
  !.............................................................................
  call toc_r (label, 1.0)
END SUBROUTINE toc_1

!> =============================================================================
!> =============================================================================
SUBROUTINE toc_i4 (label, n)
  implicit none
  character(len=*):: label
  integer:: n
  !.............................................................................
  call toc_r (label, real(n))
END SUBROUTINE toc_i4

!> =============================================================================
!> =============================================================================
SUBROUTINE toc_i8 (label, n)
  implicit none
  character(len=*):: label
  integer(8):: n
  !.............................................................................
  call toc_r (label, real(n))
END SUBROUTINE toc_i8

!> =============================================================================
FUNCTION str_f4 (t) RESULT (out)
  real(8):: t
  character(len=8):: out
  !-----------------------------------------------------------------------------
  write (out,'(f8.3)') t
  out = adjustl(out)
  out = trim(out(1:4))
END FUNCTION

!> =============================================================================
FUNCTION time_str (t) RESULT (out)
  real(8):: t
  character(len=8):: out
  real, parameter:: ns=1e-9, mu=1e-6, ms=1e-3, s=1e0, mi=6e1, hr=36e2, dy=24*hr
  real, parameter:: wk=7*dy, yr=365.*dy, mo=yr/12.
  !-----------------------------------------------------------------------------
  if      (t/mu < 1d0) then
    out=trim(str_f4(t/ns))//' ns'
    1 format(f7.3,1x,a2)
  else if (t/ms < 1d0) then
    out=trim(str_f4(t/mu))//' mus'
  else if (   t < 1d0) then
    out=trim(str_f4(t/ms))//' ms'
  else if (t/mi < 1d0) then
    out=trim(str_f4(   t))//' s '
  else if (t/hr < 1d0) then
    out=trim(str_f4(t/mi))//' mn'
  else if (t/dy < 1d0) then
    out=trim(str_f4(t/hr))//' hr'
  else if (t/wk < 1d0) then
    out=trim(str_f4(t/dy))//' dy'
  else if (t/mo < 1d0) then
    out=trim(str_f4(t/wk))//' wk'
  else if (t/yr < 1d0) then
    out=trim(str_f4(t/mo))//' mo'
  else
    out=trim(str_f4(t/yr))//' yr'
  end if
END FUNCTION

!> =============================================================================
!> =============================================================================
SUBROUTINE toc_r (label, n)
  implicit none
  character(len=*):: label
  real:: n, n_tmp
  real(8):: omp_get_wtime
  !.............................................................................
  call detect_cores
  wt = (omp_get_wtime() - wt)
  n_tmp = merge(omp_nthreads, min(n_core,omp_nthreads), n_core==0)
  write (io_unit%log,1) label, trim(time_str(wt)), &
    trim(time_str(wt*n_tmp/max(n,1.))), '/upd,', n, 'updates'
  1 format (1x, a, ':  ', a, ',  ', a, a, 1p, e10.2, 1x, a)
  flush (io_unit%log)
  if (mpi%rank==0) &
    print 1, label, trim(time_str(wt)), &
      trim(time_str(wt*n_tmp/max(n,1.))), '/upd,', n, 'updates'
  flush (6)
  wt = omp_get_wtime()
END SUBROUTINE toc_r

!> =============================================================================
!> Attempt to determine the actual number of cores
!> =============================================================================
SUBROUTINE detect_cores
  character(len=32):: file="/proc/cpuinfo"
  character(len=32):: line
  logical:: exist
  integer:: i, iostat
  !.............................................................................
  inquire (file=file, exist=exist)
  if (exist) then
    open (file=file, unit=io_unit%os, form='formatted', status='old')
    do while (.true.)
      read (io_unit%os,'(a32)',iostat=iostat) line
      if (iostat/=0) then
        exit
      else
        if (line(1:8)=='cpu core') then
          i = index(line,':')
          read (line(i+1:),*) n_core
          if (mpi%rank==0) print'(a,i4)', '   n_core =', n_core
          exit
        end if
      end if
    end do
    close (io_unit%os)
  end if
  if (mpi%rank==0) print'(a,i4)', ' n_thread =', omp_nthreads
END SUBROUTINE detect_cores

!> =============================================================================
FUNCTION wallclock() result (delta)
  real(kind=8):: time, delta
  real(kind=8), save:: offset=-1.
  real(8):: omp_get_wtime
!...............................................................................
  time = omp_get_wtime()
  if (offset == -1.) then
    offset=time
  end if
  delta = time-offset
END FUNCTION wallclock

END MODULE omp_timer_mod
