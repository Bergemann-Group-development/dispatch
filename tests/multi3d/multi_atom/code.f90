!===============================================================================
!> This mockup executes the same code, using a thread-specific input file
!===============================================================================
SUBROUTINE code
  USE omp_mod
  implicit none
  character(80):: file, line
  integer:: unit
  !-----------------------------------------------------------------------------
  ! Open a thread-specific input file.  To make sure the unit number is unique,
  ! let one thread at-a-time request a unique unit number.
  !-----------------------------------------------------------------------------
  write (file,'("input",i1,"/input.dat")') omp%thread
  !$omp critical (open_cr)
  open (newunit=unit, file=file, form='formatted', status='old')
  !$omp end critical (open_cr)
  print *,'thread',omp%thread,' is reading '//trim(file)//' with unit', unit
  !-----------------------------------------------------------------------------
  ! Read and print the 1st line
  !-----------------------------------------------------------------------------
  read (unit,'(a)') line
  print *,'thread',omp%thread,' read this: ', trim(line)
  close (unit)
END SUBROUTINE code
