!> =============================================================================
!> =============================================================================
MODULE os_mod
  USE omp_mod
  USE io_unit_mod
  implicit none
  private
  type, public:: os_t
  contains
    procedure, nopass:: mkdir
    procedure, nopass:: mkdir_no
  end type
  type(os_t), public:: os

CONTAINS

!> =============================================================================
!> Utility to create a new directory
!> =============================================================================
SUBROUTINE mkdir (dir)
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
  real, parameter:: ms=100.
  !.............................................................................
#ifdef __INTEL_COMPILER
  inquire (directory=dir, exist=exists)
  if (exists) return
  if (.not.io_unit%do_validate) then
    exists = MakeDirQQ (dir)
    if (exists) then
      print *,omp%thread,' Intel created directory '//dir
    else
      print *,omp%thread,' WARNING: Intel failed to create directory '//dir
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
          print *, mpi%rank, omp%thread,' C-binding call created directory '//dir
        exit
      else
        if (.not.io_unit%do_validate) &
          print *, mpi%rank, omp%thread,' WARNING: C-binding call failed to create directory '//dir
      end if
      call system ('mkdir -p '//dir)
      call mpi%delay (ms=ms)
      inquire (file=dir, exist=exists)
      if (exists) then
        if (.not.io_unit%do_validate) &
          print *, mpi%rank, omp%thread,' system call created directory '//dir
        exit
      else
        if (.not.io_unit%do_validate) &
          print *, mpi%rank, omp%thread,' WARNING: system call failed to create directory '//dir
      end if
    end do
  end if
  !$omp end critical (system_cr)
#endif
END SUBROUTINE mkdir

!> =============================================================================
!> Make sure the next directory is created ahead of time
!> =============================================================================
SUBROUTINE mkdir_no (iout)
  integer:: iout
  character(len=120):: filename
  !-----------------------------------------------------------------------------
  !$omp master
  write (filename,'(a,i5.5,"/")') trim(io_unit%outputname), iout
  if (io_unit%iodir/=iout) call os%mkdir (trim(filename))
  io_unit%iodir = iout
  !$omp end master
END SUBROUTINE mkdir_no
END MODULE os_mod
