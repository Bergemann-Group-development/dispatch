!===============================================================================
!===============================================================================
MODULE os_mod
  USE io_unit_mod
  USE omp_mod
  USE mpi_mod
  USE omp_timer_mod
  implicit none
  private
  type, public:: os_t
  contains
    procedure, nopass:: mkdir
    procedure, nopass:: mkdir_next_output
    procedure, nopass:: copy_file
  end type
  type(os_t), public:: os
  real, parameter:: ms=100.
CONTAINS

!===============================================================================
!> Create a new directory
!===============================================================================
SUBROUTINE mkdir (dir)
  character(len=*):: dir
  real(8):: start
  logical:: exists
  !.............................................................................
  if (mpi%rank==0) then
    call mkdir_real (dir)
  else
    start = wallclock()
    do while (wallclock()-start < 10.0)
#ifdef __INTEL_COMPILER
      inquire (directory=dir, exist=exists)
#else
      inquire (file=dir, exist=exists)
#endif
      if (exists) exit
      call omp%delay (ms=ms)
    end do
    if (.not.exists) then
      write (io_unit%output,*) 'WARNING: rank',mpi%rank,' cannot see',dir
    end if
  end if
END SUBROUTINE mkdir

!===============================================================================
!> Utility to create a new directory
!===============================================================================
SUBROUTINE mkdir_real (dir)
#ifdef __INTEL_COMPILER
  USE ifport
#else
  use iso_c_binding
  interface
    function fmkdir(path,mode) bind(c,name="mkdir")
      use iso_c_binding
      integer(c_int) :: fmkdir
      character(kind=c_char,len=1) :: path(*)
      integer(c_int16_t), value :: mode
    end function fmkdir
  end interface
  integer :: i, iter
#endif
  character(len=*):: dir
  logical:: exists
  !.............................................................................
#ifdef __INTEL_COMPILER
  inquire (directory=dir, exist=exists)
  if (exists) return
  if (.not.io_unit%do_validate) then
    exists = MakeDirQQ (dir)
    if (exists) then
      write(stderr,*) ' Intel created directory '//dir
    else
      write(stderr,'(2(a,i5,2x),a)') 'rank:', mpi%rank, 'thread:', omp%thread, &
        ' WARNING: Intel failed to create directory '//dir
    end if
  end if
#else
  !$omp critical (system_cr)
  inquire (file=dir, exist=exists)
  if (.not.exists) then
    do iter=1,10
      i = fmkdir(dir, int(o'772',c_int16_t))
      call omp%delay (ms=ms)
      inquire (file=dir, exist=exists)
      if (exists) then
        if (.not.io_unit%do_validate) &
          write(stderr,*) ' C-binding call created directory '//dir
        exit
      else
        if (.not.io_unit%do_validate) &
          write(stderr,'(2(a,i5,2x),a)') 'rank:', mpi%rank, 'thread:', omp%thread, &
            ' WARNING: C-binding call failed to create directory '//dir
      end if
      call system ('mkdir -p '//dir)
      call omp%delay (ms=ms)
      inquire (file=dir, exist=exists)
      if (exists) then
        if (.not.io_unit%do_validate) &
          write(stderr,*) ' system call created directory '//dir
        exit
      else
        if (.not.io_unit%do_validate) &
          write(stderr,'(2(a,i5,2x),a)') 'rank:', mpi%rank, 'thread:', omp%thread, &
          ' WARNING: system call failed to create directory '//dir
      end if
    end do
  end if
  !$omp end critical (system_cr)
#endif
END SUBROUTINE mkdir_real

!===============================================================================
!> Make sure the *next* directory is created ahead of time; failure to do so
!> can sometimes result in a thread trying to write output to a non-existent
!> directory!
!===============================================================================
SUBROUTINE mkdir_next_output (iout)
  integer:: iout
  character(len=120):: filename
  integer, save:: checked=-1
  !-----------------------------------------------------------------------------
  !$omp critical (mkdir_no_cr)
  if (iout > checked) then
    write (filename,'(a,i5.5,"/")') trim(io_unit%run_dir), iout
    if (io_unit%iodir/=iout) call os%mkdir (trim(filename))
    io_unit%iodir = iout
    checked = iout
  end if
  !$omp end critical (mkdir_no_cr)
END SUBROUTINE mkdir_next_output

!===============================================================================
!===============================================================================
SUBROUTINE copy_file (file1, file2)
  character(len=*):: file1, file2
  !.............................................................................
  logical:: exists
  integer:: u1, u2, n, iostat
  character(len=120):: line
  !-----------------------------------------------------------------------------
  call mpi%barrier ('copy_file')
  if (mpi%master .and. omp%master) then
    print *, 'copying '//trim(file1)//' to '//trim(file2)
    inquire (file=trim(file1), exist=exists)
    if (exists) then
      open (newunit=u1,file=trim(file1), form='formatted', status='old')
      open (newunit=u2,file=trim(file2), form='formatted', status='unknown')
      n = 0
      do while (.true.)
        read (u1,'(a)', iostat=iostat) line
        if (iostat /= 0) then
          exit
        end if
        write (u2,'(a)') trim(line)
        n = n+1
      end do
      close (u1)
      close (u2)
      print *, n,' lines copied'
    end if
  end if
  call mpi%barrier ('copy_file')
END SUBROUTINE copy_file

END MODULE os_mod
