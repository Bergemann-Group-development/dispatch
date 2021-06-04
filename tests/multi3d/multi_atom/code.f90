!===============================================================================
!> This mockup executes the same code, using a thread-specific input file
!===============================================================================
SUBROUTINE code
  USE omp_mod
  implicit none
  character(4):: dir
  character(80):: line
  integer:: unit
  !-----------------------------------------------------------------------------
  ! Open a thread-specific input file.  To make sure the unit number is unique,
  ! let one thread at-a-time request a unique unit number.
  !-----------------------------------------------------------------------------
  write (dir,'("run",i1)') omp%thread
  !$omp critical (open_cr)
  call chdir (dir)
  open (newunit=unit, file='input.dat', form='formatted', status='old')
  print *,'thread',omp%thread,' is reading in directory '//trim(dir)
  !-----------------------------------------------------------------------------
  ! Read and print the 1st line
  !-----------------------------------------------------------------------------
  read (unit,'(a)') line
  print *,'thread',omp%thread,' read this: ', trim(line)
  close (unit)
  call chdir ('..')
  !$omp end critical (open_cr)
  !-----------------------------------------------------------------------------
  ! Do some work
  !-----------------------------------------------------------------------------
  ! ...
  !-----------------------------------------------------------------------------
  ! Open a thread-specific output file.  To make sure the unit number is unique,
  ! let one thread at-a-time request a unique unit number.
  !-----------------------------------------------------------------------------
  write (dir,'("run",i1)') omp%thread
  !$omp critical (open_cr)
  call chdir (dir)
  open (newunit=unit, file='output.dat', form='formatted', status='unknown')
  print *,'thread',omp%thread,' is writing in directory '//trim(dir)
  !-----------------------------------------------------------------------------
  ! Write the line
  !-----------------------------------------------------------------------------
  write (unit,'(a)') line
  print *,'thread',omp%thread,' wrote this: ', trim(line)
  close (unit)
  call chdir ('..')
  !$omp end critical (open_cr)
END SUBROUTINE code
