!===============================================================================
!===============================================================================
MODULE os_mod
  USE iso_fortran_env, only: int8
  USE io_unit_mod
  USE omp_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE mpi_mod
  USE kinds_mod
  USE shared_mod
  implicit none
  private
  integer, parameter:: maxlabels=100
  type, public:: os_t
    integer:: verbose=0
    real(8):: gb=0.0
    real(8):: max_mem=0.0
    real(8):: gb_step=1.0
    character(len=64), dimension(maxlabels):: labels
    real(8), dimension(maxlabels):: labelled_mem
    integer, dimension(maxlabels):: n_allocated
    integer:: nlabels=0
  contains
    procedure:: mkdir
    procedure:: exists_dir
    procedure, nopass:: mkdir_next_output
    procedure, nopass:: copy_file
    procedure:: gb_mem
    procedure:: gb_print
    procedure, nopass:: allocate_1d_r
    procedure, nopass:: allocate_2d_r
    procedure, nopass:: allocate_3d_r
    procedure, nopass:: allocate_4d_r
    procedure, nopass:: allocate_5d_r
    procedure, nopass:: allocate_6d_r
    procedure, nopass:: allocate_1d_d
    procedure, nopass:: allocate_2d_d
    procedure, nopass:: allocate_3d_d
    procedure, nopass:: allocate_4d_d
    procedure, nopass:: allocate_1d_i
    procedure, nopass:: allocate_2d_i
    procedure, nopass:: allocate_3d_i
    procedure, nopass:: allocate_3d_i8
    procedure, nopass:: allocate_4d_i
    procedure, nopass:: allocate_1d_l
    procedure, nopass:: allocate_3d_l
    generic:: allocate => &
      allocate_1d_r, allocate_2d_r, allocate_3d_r, allocate_4d_r, &
      allocate_1d_d, allocate_2d_d, allocate_3d_d, allocate_4d_d, &
      allocate_1d_i, allocate_2d_i, allocate_3d_i, allocate_4d_i, &
      allocate_5d_r, allocate_6d_r, allocate_1d_l, allocate_3d_i8, &
      allocate_3d_l
    procedure, nopass:: deallocate_1d_r
    procedure, nopass:: deallocate_2d_r
    procedure, nopass:: deallocate_3d_r
    procedure, nopass:: deallocate_4d_r
    procedure, nopass:: deallocate_5d_r
    procedure, nopass:: deallocate_6d_r
    procedure, nopass:: deallocate_1d_d
    procedure, nopass:: deallocate_2d_d
    procedure, nopass:: deallocate_3d_d
    procedure, nopass:: deallocate_4d_d
    procedure, nopass:: deallocate_1d_i
    procedure, nopass:: deallocate_2d_i
    procedure, nopass:: deallocate_3d_i
    procedure, nopass:: deallocate_3d_i8
    procedure, nopass:: deallocate_4d_i
    procedure, nopass:: deallocate_1d_l
    procedure, nopass:: deallocate_3d_l
    generic:: deallocate => &
      deallocate_1d_r, deallocate_2d_r, deallocate_3d_r, deallocate_4d_r, &
      deallocate_1d_d, deallocate_2d_d, deallocate_3d_d, deallocate_4d_d, &
      deallocate_1d_i, deallocate_2d_i, deallocate_3d_i, deallocate_4d_i, &
      deallocate_5d_r, deallocate_6d_r, deallocate_1d_l, deallocate_3d_i8, &
      deallocate_3d_l
    procedure, nopass:: mem_register
    procedure:: memory_info
    procedure:: track_mem
    procedure:: print_tracked_mem
  end type
  logical, save:: initialized = .false.
  real, parameter:: ms=30.
  type(lock_t), public:: lock
  type(os_t), public:: os
CONTAINS

!===============================================================================
!> Set parameters
!===============================================================================
SUBROUTINE set_params (self, verbose)
  class(os_t):: self
  integer, optional:: verbose
  !-----------------------------------------------------------------------------
  if (present(verbose)) then
    self%verbose = verbose
  end if
END SUBROUTINE set_params

!===============================================================================
!> Create a new directory
!===============================================================================
SUBROUTINE mkdir (self, name)
  class(os_t):: self
  character(len=*):: name
  real(8):: start
  logical:: exists
  !.............................................................................
  exists = self%exists_dir (name)
  if (.not. exists) then
    if (mpi%rank==0) then
      call mkdir_real (self, name)
    else
      start = wallclock()
      do while (wallclock()-start < 10.0)
        if (exists) then
          if (self%verbose > 0) &
            write (stdout,*) 'os%mkdir: rank,wait=', mpi%rank, wallclock()-start
          exit
        end if
        call omp%delay (ms=ms)
        exists = self%exists_dir (name)
      end do
      if (.not.exists) then
        write (stderr,*) 'WARNING: rank', mpi%rank, ' cannot see ', name
      end if
    end if
  end if
END SUBROUTINE mkdir

!===============================================================================
!> Check if a directory exists
!===============================================================================
FUNCTION exists_dir (self, name) RESULT (exists)
  class(os_t):: self
  character(len=*):: name
  logical:: exists
  !-----------------------------------------------------------------------------
#ifdef __INTEL_COMPILER
  inquire (directory=name, exist=exists)
#else
  inquire (file=name, exist=exists)
#endif
END FUNCTION exists_dir

!===============================================================================
!> Utility to create a new directory
!===============================================================================
SUBROUTINE mkdir_real (self, dir)
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
  class(os_t):: self
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
      write (io_unit%output,'(a,i4,a)') &
        ' copied',n,' lines from '//trim(file1)//' to '//trim(file2)
    end if
  end if
  call mpi%barrier ('copy_file')
END SUBROUTINE copy_file

!===============================================================================
!> Register increase of memory, counted in gigabytes
!===============================================================================
SUBROUTINE gb_mem (self, bits, count, label)
  class (os_t):: self
  integer(4):: bits, count
  character(len=*), optional:: label
  !.............................................................................
  real(8):: gb, lgb
  real, save:: gb_next=0d0
  !-----------------------------------------------------------------------------
  if (gb_next==0d0) gb_next = self%gb_step
  gb = (bits/(8.*1024.**3))*count
  !$omp atomic capture
  self%gb = self%gb + gb
  lgb = self%gb
  !$omp end atomic
  if (present(label)) then
    call self%track_mem(gb, label)
  endif
  !-----------------------------------------------------------------------------
  ! When passing a GB limit, print and incremnt
  !-----------------------------------------------------------------------------
  if (abs(lgb - gb_next) > self%gb_step) then
    !$omp critical (gb_cr)
    if (abs(lgb - gb_next) > self%gb_step) then
      do while (abs(lgb - gb_next) > self%gb_step)
        call self%gb_print (real(gb_next,kind=8))
        gb_next = gb_next + sign (self%gb_step, lgb-gb_next)
      end do
      if (self%gb > self%max_mem) call self%print_tracked_mem()
    end if
    !$omp end critical (gb_cr)
  end if
  if (self%gb > self%max_mem) then
    !$omp atomic write
    self%max_mem = self%gb
  endif
END SUBROUTINE gb_mem

!===============================================================================
SUBROUTINE track_mem (self, gb, label)
  class (os_t):: self
  real(8):: gb
  character(len=*):: label
  integer :: i
  !.............................................................................
  !$omp critical (track_mem_cr)
  i = findloc(self%labels(1:self%nlabels), label, 1)
  if (i == 0) then
    if (gb > 0) then
      i = findloc(self%n_allocated, 0, 1)
      if (i > maxlabels) call shared%abort ('os_mod: increase maxlabels')
      self%nlabels = max(i, self%nlabels)
      self%labels(i) = label
      self%labelled_mem(i) = self%labelled_mem(i) + gb
      self%n_allocated(i) = self%n_allocated(i) + sign(1d0, gb)
    else
      call shared%abort ('os_mod: memory with unknown label='//trim(label)//' is deallocated.')
    endif
  endif
  self%labelled_mem(i) = self%labelled_mem(i) + gb
  self%n_allocated(i) = self%n_allocated(i) + sign(1d0, gb)
  !$omp end critical (track_mem_cr)
END SUBROUTINE track_mem

!===============================================================================
SUBROUTINE print_tracked_mem (self)
  class (os_t):: self
  real(8) :: total
  integer :: i
  !.............................................................................
  total = sum(self%labelled_mem(1:self%nlabels))
  do i=1, self%nlabels
    if (self%labelled_mem(i) / total > 0.05) then
      write(stdout, '(A, 1x, I6, f10.4)') &
            adjustr(self%labels(i)), &
            self%n_allocated(i), &
            self%labelled_mem(i)
    endif
  end do
  write(stdout, '(56x, A, 1x, f10.4)') 'rest_mem', self%gb - total
  write(stdout, '(56x, A, 1x, f10.4)') 'total', self%gb
END SUBROUTINE print_tracked_mem

!===============================================================================
!> Print memory use, gigabytes
!===============================================================================
SUBROUTINE gb_print (self, gb)
  class (os_t):: self
  real(8):: gb
  !.............................................................................
  write (stdout,'(1x,a,f8.3,a)') 'process memory allocated:', gb,' GB'
END SUBROUTINE gb_print

!===============================================================================
!> Allocate a 1-D real array
!===============================================================================
SUBROUTINE allocate_1d_r (a, m, zero, label)
  real(4), pointer, contiguous:: a(:)
  integer:: m
  logical, optional:: zero
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  allocate (a(m))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_1d_r

!===============================================================================
!> Allocate a 2-D real array
!===============================================================================
SUBROUTINE allocate_2d_r (a, m, zero, label)
  real(4), pointer:: a(:,:)
  integer:: m(:)
  logical, optional:: zero
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_2d_r

!===============================================================================
!> Allocate a 3-D real array
!===============================================================================
SUBROUTINE allocate_3d_r (a, m, zero, label)
  real(4), pointer, contiguous:: a(:,:,:)
  integer:: m(:)
  logical, optional:: zero
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_3d_r

!===============================================================================
!> Allocate a 4-D real array
!===============================================================================
SUBROUTINE allocate_4d_r (a, m, zero, label)
  real(4), pointer, contiguous:: a(:,:,:,:)
  integer:: m(:)
  logical, optional:: zero
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_4d_r

!===============================================================================
!> Allocate a 5-D real array
!===============================================================================
SUBROUTINE allocate_5d_r (a, m, zero, label)
  real(4), pointer, contiguous:: a(:,:,:,:,:)
  integer:: m(:)
  logical, optional:: zero
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4),m(5)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_5d_r

!===============================================================================
!> Allocate a 6-D real array
!===============================================================================
SUBROUTINE allocate_6d_r (a, m, zero, label)
  real(4), pointer, contiguous:: a(:,:,:,:,:,:)
  integer:: m(:)
  logical, optional:: zero
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4),m(5),m(6)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_6d_r

!===============================================================================
!> Allocate a 1-D logical array
!===============================================================================
SUBROUTINE allocate_1d_l (a, m, zero, label)
  logical, pointer, contiguous:: a(:)
  integer:: m
  logical, optional:: zero
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  allocate (a(m))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = .false.
END SUBROUTINE allocate_1d_l

!===============================================================================
!> Allocate a 3-D logical array
!===============================================================================
SUBROUTINE allocate_3d_l (a, m, zero, label)
  logical, pointer, contiguous:: a(:,:,:)
  integer:: m(:)
  logical, optional:: zero
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  allocate (a(m(1), m(2), m(3)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = .false.
END SUBROUTINE allocate_3d_l

!===============================================================================
!> Allocate a 1-D integer array
!===============================================================================
SUBROUTINE allocate_1d_i (a, m, zero, label)
  integer(4), pointer, contiguous:: a(:)
  integer:: m
  logical, optional:: zero
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  allocate (a(m))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_1d_i

!===============================================================================
!> Allocate a 2-D integer array
!===============================================================================
SUBROUTINE allocate_2d_i (a, m, zero, label)
  integer(4), pointer:: a(:,:)
  integer:: m(:)
  logical, optional:: zero
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_2d_i

!===============================================================================
!> Allocate a 3-D integer array
!===============================================================================
SUBROUTINE allocate_3d_i (a, m, zero, label)
  integer(4), pointer, contiguous:: a(:,:,:)
  integer:: m(:)
  logical, optional:: zero
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_3d_i

!===============================================================================
!> Allocate a 3-D integer array
!===============================================================================
SUBROUTINE allocate_3d_i8 (a, m, zero, label)
  integer(kind=int8), pointer, contiguous:: a(:,:,:)
  integer:: m(:)
  logical, optional:: zero
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_3d_i8

!===============================================================================
!> Allocate a 4-D integer array
!===============================================================================
SUBROUTINE allocate_4d_i (a, m, zero, label)
  integer(4), pointer, contiguous:: a(:,:,:,:)
  integer:: m(:)
  logical, optional:: zero
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_4d_i

!===============================================================================
!> Allocate a 1-D double array
!===============================================================================
SUBROUTINE allocate_1d_d (a, m, zero, label)
  real(8), pointer, contiguous:: a(:)
  integer:: m
  logical, optional:: zero
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  allocate (a(m))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_1d_d

!===============================================================================
!> Allocate a 2-D double array
!===============================================================================
SUBROUTINE allocate_2d_d (a, m, zero, label)
  real(8), pointer:: a(:,:)
  integer:: m(:)
  logical, optional:: zero
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_2d_d

!===============================================================================
!> Allocate a 3-D double array
!===============================================================================
SUBROUTINE allocate_3d_d (a, m, zero, label)
  real(8), pointer, contiguous:: a(:,:,:)
  integer:: m(:)
  logical, optional:: zero
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_3d_d

!===============================================================================
!> Allocate a 4-D double array
!===============================================================================
SUBROUTINE allocate_4d_d (a, m, zero, label)
  real(8), pointer, contiguous:: a(:,:,:,:)
  integer:: m(:)
  logical, optional:: zero
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_4d_d

!===============================================================================
!> Dealloccate a 1-D integer array
!===============================================================================
SUBROUTINE deallocate_1d_r (a, label)
  real(4), pointer, contiguous:: a(:)
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_1d_r

!===============================================================================
!> Dealloccate a 2-D integer array
!===============================================================================
SUBROUTINE deallocate_2d_r (a, label)
  real(4), pointer, contiguous:: a(:,:)
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_2d_r

!===============================================================================
!> Dealloccate a 3-D integer array
!===============================================================================
SUBROUTINE deallocate_3d_r (a, label)
  real(4), pointer, contiguous:: a(:,:,:)
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_3d_r

!===============================================================================
!> Dealloccate a 4-D integer array
!===============================================================================
SUBROUTINE deallocate_4d_r (a, label)
  real(4), pointer, contiguous:: a(:,:,:,:)
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_4d_r

!===============================================================================
!> Dealloccate a 5-D integer array
!===============================================================================
SUBROUTINE deallocate_5d_r (a, label)
  real(4), pointer, contiguous:: a(:,:,:,:,:)
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_5d_r

!===============================================================================
!> Dealloccate a 6-D integer array
!===============================================================================
SUBROUTINE deallocate_6d_r (a, label)
  real(4), pointer, contiguous:: a(:,:,:,:,:,:)
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_6d_r

!===============================================================================
!> Dealloccate a 1-D logical array
!===============================================================================
SUBROUTINE deallocate_1d_l (a, label)
  logical, pointer, contiguous:: a(:)
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_1d_l

!===============================================================================
!> Dealloccate a 3-D logical array
!===============================================================================
SUBROUTINE deallocate_3d_l (a, label)
  logical, pointer, contiguous:: a(:,:,:)
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_3d_l

!===============================================================================
!> Dealloccate a 1-D integer array
!===============================================================================
SUBROUTINE deallocate_1d_i (a, label)
  integer(4), pointer, contiguous:: a(:)
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_1d_i

!===============================================================================
!> Dealloccate a 2-D integer array
!===============================================================================
SUBROUTINE deallocate_2d_i (a, label)
  integer(4), pointer, contiguous:: a(:,:)
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_2d_i

!===============================================================================
!> Dealloccate a 3-D integer array
!===============================================================================
SUBROUTINE deallocate_3d_i (a, label)
  integer(4), pointer, contiguous:: a(:,:,:)
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_3d_i

!===============================================================================
!> Dealloccate a 3-D integer array
!===============================================================================
SUBROUTINE deallocate_3d_i8 (a, label)
  integer(kind=int8), pointer, contiguous:: a(:,:,:)
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_3d_i8

!===============================================================================
!> Dealloccate a 4-D integer array
!===============================================================================
SUBROUTINE deallocate_4d_i (a, label)
  integer(4), pointer, contiguous:: a(:,:,:,:)
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_4d_i

!===============================================================================
!> Dealloccate a 1-D double array
!===============================================================================
SUBROUTINE deallocate_1d_d (a, label)
  real(8), pointer, contiguous:: a(:)
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_1d_d

!===============================================================================
!> Dealloccate a 2-D double array
!===============================================================================
SUBROUTINE deallocate_2d_d (a, label)
  real(8), pointer, contiguous:: a(:,:)
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_2d_d

!===============================================================================
!> Dealloccate a 3-D double array
!===============================================================================
SUBROUTINE deallocate_3d_d (a, label)
  real(8), pointer, contiguous:: a(:,:,:)
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_3d_d

!===============================================================================
!> Dealloccate a 4-D double array
!===============================================================================
SUBROUTINE deallocate_4d_d (a, label)
  real(8), pointer, contiguous:: a(:,:,:,:)
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_4d_d

!===============================================================================
!> For verbose==0 no checking is done
!> For verbose==1 or missing a log message is printed if memory has increased.
!> For verbose==2 or larger a log message is always printed
!===============================================================================
SUBROUTINE memory_info (self, label, verbose)
  class(os_t):: self
  real(8):: dgb
  real(8), save:: gb=0d0
  integer, optional:: verbose
  character(*), optional:: label
  integer:: verbose_L
  !.............................................................................
  if (present(verbose)) then
    verbose_L = verbose
  else
    verbose_L = 1
  end if
  if (verbose_L > 0) then
    call lock%set
    if (gb > 0d0) then
      dgb = self%gb-gb
      if (verbose_L > 1 .or. dgb > 0d0) then
        if (present(label)) then
          write (stdout,'(1x,a,f15.9,a,f15.9,a)') &
            'memory allocated:', self%gb,' GB, change:', dgb, ' GB, at '//trim(label)
        else
          write (stdout,'(1x,a,f15.9,a,f15.9,a)') &
            'memory allocated:', self%gb,' GB, change:', dgb, ' GB'
        end if
      end if
    end if
    gb = self%gb
    call lock%unset
  end if
END SUBROUTINE memory_info

!===============================================================================
!> Register increase of memory, counted in storage_size bits, times a count
!===============================================================================
SUBROUTINE mem_register (bits, count, label)
  integer:: bits, count
  character(len=*), optional:: label
  integer:: unit
  !-----------------------------------------------------------------------------
  if (os%verbose > 0) then
    unit = merge (stdout, io_unit%mpi, os%verbose > 1)
    if (present(label)) then
      write(unit,1) bits,' bits per word', count,' words, for', &
        (bits/(8.*1024.**3))*count, os%gb, ' GB '//trim(label)
      1 format("mem_reg: ",i4,a,i9,a,2f8.3,a)
    else
      write(unit,1) bits,' bits per word', count,' words, for', &
        (bits/(8.*1024.**3))*count, os%gb, ' GB'
    end if
  end if
  call os%gb_mem (bits, count, label)
END SUBROUTINE mem_register

!===============================================================================
!> Local utility
!===============================================================================
LOGICAL FUNCTION option (o)
  logical, optional:: o
  if (present(o)) then
    option = o
  else
    option = .false.
  end if
END FUNCTION option

END MODULE os_mod
