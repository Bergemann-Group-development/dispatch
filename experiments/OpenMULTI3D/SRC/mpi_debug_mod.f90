#include "../HOP/read_namelist.inc"
!===============================================================================
!> Per-rank process ids and an optional halt for attaching a debugger to MPI runs.
!===============================================================================
MODULE mpi_debug_mod
  USE io_mod
  USE namelist_check_mod, only: namelist_error
  USE mpi_mod
#ifdef __INTEL_COMPILER
  USE IFPORT
#endif
  implicit none

contains

  !===============================================================================
  !> In a run with more than one MPI rank, writes each rank's process id to
  !> `./data/rank_NNNNN.pid` and reads &mpi_debug_params. With debug set, all ranks
  !> then wait after a barrier in an endless loop, until a debugger attached to the
  !> process sets halt to a nonzero value.
  !===============================================================================
  SUBROUTINE check_mpi_debug()
    integer :: lunit
    integer :: halt=0
    logical :: debug=.false.
    character(len=128) :: pid_file
    namelist /mpi_debug_params/ debug
    !-----------------------------------------------------------------------------
    if (mpi%size==1) return

    write(pid_file,'(a,I5.5,a)') './data/rank_', mpi%rank,'.pid'
    open(newunit=lunit, file=pid_file, form='formatted')
    write(lunit, '(I0)') getpid()
    close(lunit)

    READ_NAMELIST(mpi_debug_params)
    write (io%namelist, mpi_debug_params)

    if (.not. debug) return

    call mpi%barrier ('debug')
    print *, 'Programm is halted in mpi_debug_mod. Waiting...'
    do while (halt==0)
      call sleep(5)
    end do

  END SUBROUTINE check_mpi_debug

END MODULE mpi_debug_mod
