MODULE mpi_debug_mod
  USE io_mod
  USE mpi_mod
  USE IFPORT
  implicit none
    
contains

!===============================================================================
SUBROUTINE check_mpi_debug()
  integer :: lunit, iostat
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

  rewind (io%input)
  read(io%input, mpi_debug_params, iostat=iostat)
  call io%namelist_warning ('mpi_debug_params', iostat)
  write (io%namelist, mpi_debug_params)

  if (.not. debug) return

  call mpi%barrier ('debug')
  print *, 'Programm is halted in mpi_debug_mod. Waiting...'
  do while (halt==0)
    call sleep(5)
  end do

END SUBROUTINE check_mpi_debug
    
END MODULE mpi_debug_mod