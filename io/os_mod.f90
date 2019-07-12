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
    procedure, nopass:: mkdir_no
  end type
  type(os_t), public:: os
  real, parameter:: ms=100.
CONTAINS

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
      call mpi%delay (ms=ms)
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
      call mpi%delay (ms=ms)
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
      call mpi%delay (ms=ms)
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
!> Make sure the next directory is created ahead of time
!===============================================================================
SUBROUTINE mkdir_no (iout)
  integer:: iout
  character(len=120):: filename
  integer, save:: checked=-1
  !-----------------------------------------------------------------------------
  !$omp critical (mkdir_no_cr)
  if (iout > checked) then
    write (filename,'(a,i5.5,"/")') trim(io_unit%outputname), iout
    if (io_unit%iodir/=iout) call os%mkdir (trim(filename))
    io_unit%iodir = iout
    checked = iout
  end if
  !$omp end critical (mkdir_no_cr)
END SUBROUTINE mkdir_no
END MODULE os_mod
