!===============================================================================
!===============================================================================
MODULE os_mod
  USE io_unit_mod
  USE omp_mod
  USE mpi_mod
  USE omp_timer_mod
  USE kinds_mod
  implicit none
  private
  type, public:: os_t
    real(8):: gb=0.0
  contains
    procedure, nopass:: mkdir
    procedure, nopass:: mkdir_next_output
    procedure, nopass:: copy_file
    procedure:: allocate2
    procedure:: allocate3
    procedure:: allocate4
    procedure:: allocate5
    generic:: allocate => allocate2,allocate3,allocate4,allocate5
    procedure:: deallocate2
    procedure:: deallocate3
    procedure:: deallocate4
    procedure:: deallocate5
    generic:: deallocate => deallocate2,deallocate3,deallocate4,deallocate5
    procedure:: gb_mem
    procedure:: gb_print
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
      if (exists) then
        write (stdout,*) 'os%mkdir: rank,wait=', mpi%rank, wallclock()-start
        exit
      end if
      call omp%delay (ms=ms)
    end do
    if (.not.exists) then
      write (stderr,*) 'WARNING: rank',mpi%rank,' cannot see ',dir
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
  !-----------------------------------------------------------------------------
  ! Note that the calling code makes sure only one thread executes this
  !-----------------------------------------------------------------------------
#ifdef __INTEL_COMPILER
  inquire (directory=dir, exist=exists)
  if (exists) return
  if (.not.io_unit%do_validate) then
    exists = MakeDirQQ (dir)
    if (exists) then
      write(stderr,*) 'Intel created directory '//dir
    else
      write(stderr,'(2(a,i5,2x),a)') 'rank:', mpi%rank, 'thread:', omp%thread, &
        ' WARNING: Intel failed to create directory '//dir
    end if
  end if
#else
  inquire (file=dir, exist=exists)
  if (.not.exists) then
    !---------------------------------------------------------------------------
    ! First try the C-binding, wating up to 1s for a result
    !---------------------------------------------------------------------------
    i = fmkdir(dir, int(o'772',c_int16_t))
    do iter=1,10
      inquire (file=dir, exist=exists)
      if (exists) then
        if (.not.io_unit%do_validate) &
          write(stderr,*) 'C-binding call created directory '//dir, iter
        exit
      end if
      i = fmkdir(dir, int(o'772',c_int16_t))
      call omp%delay (ms=ms)
    end do
    inquire (file=dir, exist=exists)
    if (.not. exists) then
      if (.not.io_unit%do_validate) then
        write(stderr,'(2(a,i5,2x),a)') 'rank:', mpi%rank, 'thread:', omp%thread, &
          ' WARNING: C-binding call failed to create directory '//dir
      end if
      !-------------------------------------------------------------------------
      ! If the dir does not exist, fall back to a system call
      !-------------------------------------------------------------------------
      call system ('mkdir -p '//dir)
      do iter=1,10
        inquire (file=dir, exist=exists)
        if (exists) then
          if (.not.io_unit%do_validate) &
            write(stderr,*) 'system call created directory '//dir, iter
          exit
        end if
        call system ('mkdir -p '//dir)
        call omp%delay (ms=ms)
      end do
      if (.not. exists) then
        if (.not.io_unit%do_validate) &
          write(stderr,'(2(a,i5,2x),a)') 'rank:', mpi%rank, 'thread:', omp%thread, &
          ' WARNING: system call failed to create directory '//dir
      end if
    end if
  end if
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

SUBROUTINE allocate2 (self, a, m)
  class(os_t):: self
  real(var_k), pointer, contiguous:: a(:,:)
  integer:: m(2)
  !.............................................................................
  allocate (a(m(1),m(2)))
  call os%gb_mem (storage_size(a), product(shape(a)))
END SUBROUTINE allocate2

SUBROUTINE allocate3 (self, a, m)
  class(os_t):: self
  real(var_k), pointer, contiguous:: a(:,:,:)
  integer:: m(3)
  !.............................................................................
  allocate (a(m(1),m(2),m(3)))
  call os%gb_mem (storage_size(a), product(shape(a)))
END SUBROUTINE allocate3

SUBROUTINE allocate4 (self, a, m)
  class(os_t):: self
  real(var_k), pointer, contiguous:: a(:,:,:,:)
  integer:: m(4)
  !.............................................................................
  allocate (a(m(1),m(2),m(3),m(4)))
  call os%gb_mem (storage_size(a), product(shape(a)))
END SUBROUTINE allocate4

SUBROUTINE allocate5 (self, a, m)
  class(os_t):: self
  real(var_k), pointer, contiguous:: a(:,:,:,:,:)
  integer:: m(5)
  !.............................................................................
  allocate (a(m(1),m(2),m(3),m(4),m(5)))
  call os%gb_mem (storage_size(a), product(shape(a)))
END SUBROUTINE allocate5

SUBROUTINE deallocate2 (self, a)
  class(os_t):: self
  real(var_k), pointer, contiguous:: a(:,:)
  !.............................................................................
  call os%gb_mem (-storage_size(a), product(shape(a)))
  deallocate (a)
END SUBROUTINE deallocate2

SUBROUTINE deallocate3 (self, a)
  class(os_t):: self
  real(var_k), pointer, contiguous:: a(:,:,:)
  !.............................................................................
  call os%gb_mem (-storage_size(a), product(shape(a)))
  deallocate (a)
END SUBROUTINE deallocate3

SUBROUTINE deallocate4 (self, a)
  class(os_t):: self
  real(var_k), pointer, contiguous:: a(:,:,:,:)
  !.............................................................................
  call os%gb_mem (-storage_size(a), product(shape(a)))
  deallocate (a)
END SUBROUTINE deallocate4

SUBROUTINE deallocate5 (self, a)
  class(os_t):: self
  real(var_k), pointer, contiguous:: a(:,:,:,:,:)
  !.............................................................................
  call os%gb_mem (-storage_size(a), product(shape(a)))
  deallocate (a)
END SUBROUTINE deallocate5

!===============================================================================
!> Register increase of memory, counted in gigabytes
!===============================================================================
SUBROUTINE gb_mem (self, bits, count)
  class (os_t):: self
  integer(4):: bits, count
  real(8):: gb, lgb
  real, save:: gb_next=1.0
  !-----------------------------------------------------------------------------
  gb = (bits/(8.*1024.**3))*count
  !$omp atomic capture
  self%gb = self%gb + gb
  lgb = self%gb
  !$omp end atomic
  if (abs(lgb - gb_next) > 1.0) then
    !$omp critical (gb_cr)
    if (abs(lgb - gb_next) > 1.0) then
      do while (abs(lgb - gb_next) > 1.0)
        call self%gb_print (lgb)
        gb_next = gb_next + sign(1d0,lgb-gb_next)
      end do
    end if
    !$omp end critical (gb_cr)
  end if
END SUBROUTINE gb_mem

!===============================================================================
!> Print memory use, gigabytes
!===============================================================================
SUBROUTINE gb_print (self, gb)
  class (os_t):: self
  real(8):: gb
  !.............................................................................
  write (stdout,'(1x,a,f8.3,a)') 'process memory allocated:', gb,' GB'
END SUBROUTINE gb_print

END MODULE os_mod
