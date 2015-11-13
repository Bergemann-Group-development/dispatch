! $Id$
!*******************************************************************************
MODULE mpi_basic
  implicit none
  integer:: mpi_err             ! MPI error status
  logical:: mpi_master          ! true on MPI master (rank 0)
  logical:: master              ! true on OMP master thread on MPI master process
  logical:: do_trace=.false.    ! debug trace option
PUBLIC
  integer:: mpi_size            ! number of MPI processes
  integer:: mpi_rank            ! process number (zero-based)
  integer, parameter:: mch=120  ! max character string size
CONTAINS

!*******************************************************************************
SUBROUTINE init_mpi
  implicit none
  character(len=mch):: id='nompi.f90 $Id$'

!-------------------------------------------------------------------------------
!  Start up MPI and get number of nodes and our node rank
!-------------------------------------------------------------------------------
  mpi_size = 1
  mpi_rank = 0
  mpi_master = (mpi_rank == 0)
  master = mpi_master
  call print_id (id)
END SUBROUTINE init_mpi

!*******************************************************************************
SUBROUTINE end_mpi
  implicit none
!..............................................................................
  call exit        ! avoid "stop" which may produce one output line per process
END SUBROUTINE

!*******************************************************************************
SUBROUTINE abort_mpi
  implicit none
!...............................................................................
  call exit
END SUBROUTINE abort_mpi

!*******************************************************************************
SUBROUTINE barrier_mpi (label)
  implicit none
  character(len=*) label
!...............................................................................
END SUBROUTINE barrier_mpi

!*******************************************************************************
REAL(kind=8) FUNCTION wallclock()
  implicit none
  integer, save:: count, count_rate=0, count_max
  real(kind=8), save:: offset=0.
!
  if (count_rate == 0) then
    call system_clock(count=count, count_rate=count_rate, count_max=count_max)
    offset = count/real(count_rate)
  else
    call system_clock(count=count)
  end if
  wallclock = count/real(count_rate) - offset
END FUNCTION wallclock

!*******************************************************************************
SUBROUTINE print_id (id)
  implicit none
  character(len=mch) id
  character(len=mch), save:: hl= &
    '------------------------------------------------------------------------------------------'
!..............................................................................
  if (id .ne. ' ' .or. do_trace) then
    if (master) then
      print'(a)',hl
      print'(1x,a,f12.2)', trim(id), wallclock()
    end if
    if (.not.do_trace) id = ' '
  end if
END SUBROUTINE print_id

END MODULE mpi_basic
