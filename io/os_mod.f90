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
    real(8), allocatable:: gb_thread(:)
    character(len=64), dimension(maxlabels):: labels
    real(8), dimension(maxlabels):: labelled_mem
    integer, dimension(maxlabels):: n_allocated
    integer:: nlabels=0
  contains
    procedure:: init
    procedure:: set_verbose
    procedure:: mkdir
    procedure:: exists
    procedure:: exists_dir
    procedure:: check_file
    procedure:: get_filenames
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
    procedure, nopass:: allocate_5d_d
    procedure, nopass:: allocate_1d_ra
    procedure, nopass:: allocate_2d_ra
    procedure, nopass:: allocate_3d_ra
    procedure, nopass:: allocate_4d_ra
    procedure, nopass:: allocate_5d_ra
    procedure, nopass:: allocate_1d_da
    procedure, nopass:: allocate_2d_da
    procedure, nopass:: allocate_3d_da
    procedure, nopass:: allocate_4d_da
    procedure, nopass:: allocate_5d_da
    procedure, nopass:: allocate_1d_i
    procedure, nopass:: allocate_2d_i
    procedure, nopass:: allocate_3d_i
    procedure, nopass:: allocate_2d_i2
    procedure, nopass:: allocate_3d_i2
    procedure, nopass:: allocate_3d_int8
    procedure, nopass:: allocate_1d_l
    procedure, nopass:: allocate_3d_l
    procedure, nopass:: allocate_4d_l
    procedure, nopass:: allocate_4d_i
    generic:: allocate => &
      allocate_1d_r, allocate_2d_r, allocate_3d_r, allocate_4d_r, &
      allocate_1d_d, allocate_2d_d, allocate_3d_d, allocate_4d_d, &
      allocate_5d_d, &
      allocate_1d_i, allocate_2d_i, allocate_3d_i, allocate_4d_i, &
      allocate_5d_r, allocate_6d_r, allocate_3d_l, allocate_1d_l, &
      allocate_3d_int8, allocate_4d_l, allocate_3d_i2, allocate_2d_i2, &
      allocate_1d_ra, allocate_2d_ra, allocate_3d_ra, allocate_4d_ra, &
      allocate_5d_ra, &
      allocate_1d_da, allocate_2d_da, allocate_3d_da, allocate_4d_da, &
      allocate_5d_da
    procedure, nopass:: deallocate_1d_r
    procedure, nopass:: deallocate_2d_r
    procedure, nopass:: deallocate_3d_r
    procedure, nopass:: deallocate_4d_r
    procedure, nopass:: deallocate_5d_r
    procedure, nopass:: deallocate_6d_r
    procedure, nopass:: deallocate_1d_ra
    procedure, nopass:: deallocate_2d_ra
    procedure, nopass:: deallocate_3d_ra
    procedure, nopass:: deallocate_4d_ra
    procedure, nopass:: deallocate_5d_ra
    procedure, nopass:: deallocate_1d_d
    procedure, nopass:: deallocate_2d_d
    procedure, nopass:: deallocate_3d_d
    procedure, nopass:: deallocate_4d_d
    procedure, nopass:: deallocate_5d_d
    procedure, nopass:: deallocate_1d_da
    procedure, nopass:: deallocate_2d_da
    procedure, nopass:: deallocate_3d_da
    procedure, nopass:: deallocate_4d_da
    procedure, nopass:: deallocate_5d_da
    procedure, nopass:: deallocate_1d_i
    procedure, nopass:: deallocate_2d_i
    procedure, nopass:: deallocate_3d_i
    procedure, nopass:: deallocate_2d_i2
    procedure, nopass:: deallocate_3d_i2
    procedure, nopass:: deallocate_3d_int8
    procedure, nopass:: deallocate_1d_l
    procedure, nopass:: deallocate_3d_l
    procedure, nopass:: deallocate_4d_l
    procedure, nopass:: deallocate_4d_i
    generic:: deallocate => &
      deallocate_1d_r, deallocate_2d_r, deallocate_3d_r, deallocate_4d_r, &
      deallocate_1d_d, deallocate_2d_d, deallocate_3d_d, deallocate_4d_d, &
      deallocate_5d_d, &
      deallocate_1d_i, deallocate_2d_i, deallocate_3d_i, deallocate_4d_i, &
      deallocate_5d_r, deallocate_6d_r, deallocate_3d_l, deallocate_1d_l, &
      deallocate_3d_int8, deallocate_4d_l, deallocate_3d_i2, deallocate_2d_i2, &
      deallocate_1d_ra, deallocate_2d_ra, deallocate_3d_ra, deallocate_4d_ra, &
      deallocate_5d_ra, &
      deallocate_1d_da, deallocate_2d_da, deallocate_3d_da, deallocate_4d_da, &
      deallocate_5d_da
    procedure, nopass:: mem => mem_register
    procedure, nopass:: mem_register
    procedure:: memory_info
    procedure:: track_mem
    procedure:: print_tracked_mem
    procedure, nopass:: isinf
    procedure, nopass:: isnan4
    procedure, nopass:: isnan8
    generic:: isnan => isnan4, isnan8
  end type
  logical, save:: initialized = .false.
  real, parameter:: ms=30.
  type(lock_t), public:: lock
  type(os_t), public:: os
CONTAINS

!===============================================================================
!> Initialize
!===============================================================================
SUBROUTINE init (self, verbose)
  class(os_t):: self
  integer, optional:: verbose
  !-----------------------------------------------------------------------------
  if (.not. lock%initialized) &
    call lock%init
  call                                                       lock%set('os%init')
  if (.not. allocated(self%gb_thread)) then
    allocate (self%gb_thread(0:omp%nthreads-1))
  end if
  call                                                     lock%unset('os%init')
  if (present(verbose)) then
    self%verbose = verbose
  end if
END SUBROUTINE init

!===============================================================================
!> Set parameters
!===============================================================================
SUBROUTINE set_verbose (self, verbose)
  class(os_t):: self
  integer, optional:: verbose
  !-----------------------------------------------------------------------------
  if (present(verbose)) then
    self%verbose = verbose
  end if
END SUBROUTINE set_verbose

!===============================================================================
!> Create a new directory
!===============================================================================
SUBROUTINE mkdir (self, name)
  class(os_t):: self
  character(*):: name
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
          if (self%verbose > 1) &
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
FUNCTION exists (self, name)
  class(os_t):: self
  character(*):: name
  logical:: exists
  !-----------------------------------------------------------------------------
  inquire (file=name, exist=exists)
END FUNCTION exists

!===============================================================================
!> Check if a directory exists
!===============================================================================
FUNCTION exists_dir (self, name) RESULT (exists)
  class(os_t):: self
  character(*):: name
  logical:: exists
  !-----------------------------------------------------------------------------
#ifdef __INTEL_COMPILER
  inquire (directory=name, exist=exists)
#else
  inquire (file=name, exist=exists)
#endif
END FUNCTION exists_dir

!===============================================================================
!> Check if a file exists and otherwise aborts
!===============================================================================
SUBROUTINE check_file (self, name)
  class(os_t):: self
  character(*):: name
  logical:: exists
  !-----------------------------------------------------------------------------
  inquire (file=name, exist=exists)
  if (.not. exists) then
    write(*,*) 'Error: File not found - ', trim(name)
    error stop 'File not found'
  endif
END SUBROUTINE check_file

!===============================================================================
!> Obtain array of files within a directory
!===============================================================================
subroutine get_filenames(self, directory, filenames, pattern)
  class(os_t):: self
  character(len=*), intent(in) :: directory
  character(len=256), allocatable, intent(out) :: filenames(:)
  character(len=*), intent(in), optional :: pattern
  integer :: num_files
  character(len=512) :: command
  character(len=256) :: line
  integer :: unit, iostat, status, i
  !-----------------------------------------------------------------------------
  num_files = 0
    
  ! Build command with optional pattern filtering
  if (present(pattern)) then
      write(command, '(A)') 'ls -1 "' // trim(directory) // '"/' // &
                            trim(pattern) // ' 2>/dev/null | xargs -n1 basename > .temp_files'
  else
      write(command, '(A)') 'ls -1 "' // trim(directory) // '" > .temp_files'
  end if
  
  call execute_command_line(command, wait=.true., exitstat=status)
  if (status /= 0) return
  
  ! Count files
  open(newunit=unit, file='.temp_files', status='old', action='read')
  do
      read(unit, '(A)', iostat=iostat) line
      if (iostat /= 0) exit
      if (len_trim(line) > 0) num_files = num_files + 1
  end do
  close(unit)
  
  ! Read filenames
  allocate(filenames(num_files))
  open(newunit=unit, file='.temp_files', status='old', action='read')
  do i = 1, num_files
      read(unit, '(A)', iostat=iostat) filenames(i)
      if (iostat /= 0) exit
  end do
  close(unit)
  
  call execute_command_line('rm .temp_files', wait=.true.)
end subroutine get_filenames

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
  character(*):: dir
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
  character(*):: file1, file2
  !.............................................................................
  logical:: exists
  integer:: u1, u2, n, iostat
  character(len=512):: line
  !-----------------------------------------------------------------------------
  call mpi%barrier ('copy_file')
  if (mpi%master .and. omp%master) then
    inquire (file=trim(file1), exist=exists)
    if (exists) then
      write(line, *) 'cp ', trim(file1), ' ', trim(file2)
      call execute_command_line(line)
      !open (newunit=u1,file=trim(file1), form='formatted', status='old')
      !open (newunit=u2,file=trim(file2), form='formatted', status='unknown')
      !n = 0
      !do while (.true.)
      !  read (u1,'(a)', iostat=iostat) line
      !  if (iostat /= 0) then
      !    exit
      !  end if
      !  write (u2,'(a)') trim(line)
      !  n = n+1
      !end do
      !close (u1)
      !close (u2)
      !write (io_unit%output,'(a,i4,a)') &
      !  ' copied',n,' lines from '//trim(file1)//' to '//trim(file2)
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
  ! Increment atomically
  !-----------------------------------------------------------------------------
  gb = (bits/(8.*1024.**3))*count
  self%gb_thread(omp%thread) = self%gb_thread(omp%thread) + gb
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
        gb_next = gb_next + sign (self%gb_step, lgb-gb_next)
        call self%gb_print (real(gb_next,kind=8))
      end do
      if (self%gb > self%max_mem) call self%print_tracked_mem()
    end if
    !$omp end critical (gb_cr)
  end if
  if (self%gb > self%max_mem) then
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
    if (gb >= 0) then
      i = findloc(self%n_allocated, 0, 1)
      if (i > maxlabels) call shared%abort ('os_mod: increase maxlabels')
      self%nlabels = max(i, self%nlabels)
      self%labels(i) = label
    else
      call shared%abort ('os_mod: memory with unknown label = '//trim(label)//' is deallocated.')
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
    if (self%labelled_mem(i) / total > 0.005) then
      write(stdout, '(A, 1x, I6, f10.4)') &
            adjustr(self%labels(i)), &
            self%n_allocated(i), &
            self%labelled_mem(i)
    endif
  end do
  write(stdout, '(56x, A8, 7x, f10.4)') 'rest_mem', self%gb - total
  write(stdout, '(56x, A8, 7x, f10.4)') '   total', self%gb
END SUBROUTINE print_tracked_mem

!===============================================================================
!> Print memory use, gigabytes
!===============================================================================
SUBROUTINE gb_print (self, gb)
  class (os_t):: self
  real(8):: gb
  !.............................................................................
  write (stdout,'(1x,a,f8.3,a)') 'process memory allocated:', gb,' GB'
  write (stdlog,'(1x,a,f8.3,a)') 'process memory allocated:', gb,' GB'
END SUBROUTINE gb_print

!===============================================================================
!> Allocate a 1-D real array
!===============================================================================
SUBROUTINE allocate_1d_r (a, m, label, zero)
  real, pointer:: a(:)
  integer:: m
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_1d_r

!===============================================================================
!> Allocate a 2-D real array
!===============================================================================
SUBROUTINE allocate_2d_r (a, m, label, zero)
  real(4), pointer:: a(:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_2d_r

!===============================================================================
!> Allocate a 3-D real array
!===============================================================================
SUBROUTINE allocate_3d_r (a, m, label, zero)
  real(4), pointer:: a(:,:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_3d_r

!===============================================================================
!> Allocate a 1-D logical array
!===============================================================================
SUBROUTINE allocate_1d_l (a, m, label, zero)
  logical, pointer:: a(:)
  integer:: m
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = .false.
END SUBROUTINE allocate_1d_l

!===============================================================================
!> Allocate a 3-D logical array
!===============================================================================
SUBROUTINE allocate_3d_l (a, m, label, zero)
  logical, pointer:: a(:,:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = .false.
END SUBROUTINE allocate_3d_l

!===============================================================================
!> Allocate a 4-D logical array
!===============================================================================
SUBROUTINE allocate_4d_l (a, m, label, zero)
  logical, pointer:: a(:,:,:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = .false.
END SUBROUTINE allocate_4d_l

!===============================================================================
!> Allocate a 4-D real array
!===============================================================================
SUBROUTINE allocate_4d_r (a, m, label, zero)
  real(4), pointer:: a(:,:,:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_4d_r

!===============================================================================
!> Allocate a 5-D real array
!===============================================================================
SUBROUTINE allocate_5d_r (a, m, label, zero)
  real(4), pointer:: a(:,:,:,:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4),m(5)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_5d_r

!===============================================================================
!> Allocate a 6-D real array
!===============================================================================
SUBROUTINE allocate_6d_r (a, m, label, zero)
  real(4), pointer:: a(:,:,:,:,:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4),m(5),m(6)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_6d_r

!===============================================================================
!> Allocate a 1-D integer array
!===============================================================================
SUBROUTINE allocate_1d_i (a, m, label, zero)
  integer, pointer:: a(:)
  integer:: m
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_1d_i

!===============================================================================
!> Allocate a 2-D integer array
!===============================================================================
SUBROUTINE allocate_2d_i (a, m, label, zero)
  integer(4), pointer:: a(:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_2d_i

!===============================================================================
!> Allocate a 3-D integer array
!===============================================================================
SUBROUTINE allocate_3d_i (a, m, label, zero)
  integer(4), pointer:: a(:,:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_3d_i

!===============================================================================
!> Allocate a 3-D integer array
!===============================================================================
SUBROUTINE allocate_2d_i2 (a, m, label, zero)
  integer(2), pointer:: a(:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_2d_i2

!===============================================================================
!> Allocate a 3-D integer array
!===============================================================================
SUBROUTINE allocate_3d_i2 (a, m, label, zero)
  integer(2), pointer:: a(:,:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_3d_i2

!===============================================================================
!> Allocate a 1-D integer array
!===============================================================================
SUBROUTINE allocate_3d_int8 (a, m, label, zero)
  integer(kind=int8), pointer:: a(:,:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_3d_int8

!===============================================================================
!> Allocate a 4-D integer array
!===============================================================================
SUBROUTINE allocate_4d_i (a, m, label, zero)
  integer(4), pointer:: a(:,:,:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_4d_i

!===============================================================================
!> Allocate a 1-D double array
!===============================================================================
SUBROUTINE allocate_1d_d (a, m, label, zero)
  real(8), pointer:: a(:)
  integer:: m
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_1d_d

!===============================================================================
!> Allocate a 2-D double array
!===============================================================================
SUBROUTINE allocate_2d_d (a, m, label, zero)
  real(8), pointer:: a(:,:)
  integer, intent(in) :: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_2d_d

!===============================================================================
!> Allocate a 3-D double array
!===============================================================================
SUBROUTINE allocate_3d_d (a, m, label, zero)
  real(8), pointer:: a(:,:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_3d_d

!===============================================================================
!> Allocate a 4-D double array
!===============================================================================
SUBROUTINE allocate_4d_d (a, m, label, zero)
  real(8), pointer:: a(:,:,:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_4d_d

!===============================================================================
!> Allocate a 5-D double array
!===============================================================================
SUBROUTINE allocate_5d_d (a, m, label, zero)
  real(8), pointer:: a(:,:,:,:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4),m(5)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_5d_d

!===============================================================================
!> Allocate a 1-D double array
!===============================================================================
SUBROUTINE allocate_1d_ra (a, m, label, zero)
  real(4), allocatable:: a(:)
  integer:: m
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_1d_ra

!===============================================================================
!> Allocate a 2-D double array
!===============================================================================
SUBROUTINE allocate_2d_ra (a, m, label, zero)
  real(4), allocatable:: a(:,:)
  integer, intent(in) :: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_2d_ra

!===============================================================================
!> Allocate a 3-D double array
!===============================================================================
SUBROUTINE allocate_3d_ra (a, m, label, zero)
  real(4), allocatable:: a(:,:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_3d_ra

!===============================================================================
!> Allocate a 4-D double array
!===============================================================================
SUBROUTINE allocate_4d_ra (a, m, label, zero)
  real(4), allocatable:: a(:,:,:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_4d_ra

!===============================================================================
!> Allocate a 5-D double array
!===============================================================================
SUBROUTINE allocate_5d_ra (a, m, label, zero)
  real(4), allocatable:: a(:,:,:,:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4),m(5)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_5d_ra

!===============================================================================
!> Allocate a 1-D double array
!===============================================================================
SUBROUTINE allocate_1d_da (a, m, label, zero)
  real(8), allocatable:: a(:)
  integer:: m
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_1d_da

!===============================================================================
!> Allocate a 2-D double array
!===============================================================================
SUBROUTINE allocate_2d_da (a, m, label, zero)
  real(8), allocatable:: a(:,:)
  integer, intent(in) :: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_2d_da

!===============================================================================
!> Allocate a 3-D double array
!===============================================================================
SUBROUTINE allocate_3d_da (a, m, label, zero)
  real(8), allocatable:: a(:,:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_3d_da

!===============================================================================
!> Allocate a 4-D double array
!===============================================================================
SUBROUTINE allocate_4d_da (a, m, label, zero)
  real(8), allocatable:: a(:,:,:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_4d_da

!===============================================================================
!> Allocate a 5-D double array
!===============================================================================
SUBROUTINE allocate_5d_da (a, m, label, zero)
  real(8), allocatable:: a(:,:,:,:,:)
  integer:: m(:)
  character(*), optional:: label
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4),m(5)))
  call mem_register (storage_size(a), product(shape(a)), label)
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_5d_da

!===============================================================================
!> Dealloccate a 1-D integer array
!===============================================================================
SUBROUTINE deallocate_1d_r (a, label)
  real, pointer:: a(:)
  character(*), optional:: label
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
  real(4), pointer:: a(:,:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_2d_r

!===============================================================================
!> Deallocate a 3-D integer array
!===============================================================================
SUBROUTINE deallocate_3d_r (a, label)
  real(4), pointer:: a(:,:,:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_3d_r

!===============================================================================
!> Dealloccate a 3-D integer array
!===============================================================================
SUBROUTINE deallocate_1d_l (a, label)
  logical, pointer:: a(:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_1d_l

!===============================================================================
!> Dealloccate a 3-D integer array
!===============================================================================
SUBROUTINE deallocate_3d_l (a, label)
  logical, pointer:: a(:,:,:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_3d_l

!===============================================================================
!> Dealloccate a 3-D integer array
!===============================================================================
SUBROUTINE deallocate_4d_l (a, label)
  logical, pointer:: a(:,:,:,:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_4d_l

!===============================================================================
!> Dealloccate a 4-D integer array
!===============================================================================
SUBROUTINE deallocate_4d_r (a, label)
  real(4), pointer:: a(:,:,:,:)
  character(*), optional:: label
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
  real(4), pointer:: a(:,:,:,:,:)
  character(*), optional:: label
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
  real(4), pointer:: a(:,:,:,:,:,:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_6d_r

!===============================================================================
!> Dealloccate a 1-D integer array
!===============================================================================
SUBROUTINE deallocate_1d_i (a, label)
  integer, pointer:: a(:)
  character(*), optional:: label
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
  integer(4), pointer:: a(:,:)
  character(*), optional:: label
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
  integer(4), pointer:: a(:,:,:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_3d_i

!===============================================================================
!> Dealloccate a 3-D integer array
!===============================================================================
SUBROUTINE deallocate_2d_i2 (a, label)
  integer(2), pointer:: a(:,:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_2d_i2

!===============================================================================
!> Dealloccate a 3-D integer array
!===============================================================================
SUBROUTINE deallocate_3d_i2 (a, label)
  integer(2), pointer:: a(:,:,:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_3d_i2

!===============================================================================
!> Dealloccate a 3-D integer array
!===============================================================================
SUBROUTINE deallocate_3d_int8 (a, label)
  integer(kind=int8), pointer:: a(:,:,:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_3d_int8

!===============================================================================
!> Dealloccate a 4-D integer array
!===============================================================================
SUBROUTINE deallocate_4d_i (a, label)
  integer(4), pointer:: a(:,:,:,:)
  character(*), optional:: label
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
  real(8), pointer:: a(:)
  character(*), optional:: label
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
  real(8), pointer:: a(:,:)
  character(*), optional:: label
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
  real(8), pointer:: a(:,:,:)
  character(*), optional:: label
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
  real(8), pointer:: a(:,:,:,:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_4d_d

!===============================================================================
!> Dealloccate a 5-D double array
!===============================================================================
SUBROUTINE deallocate_5d_d (a, label)
  real(8), pointer:: a(:,:,:,:,:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_5d_d

!===============================================================================
!> Dealloccate a 1-D double array
!===============================================================================
SUBROUTINE deallocate_1d_da (a, label)
  real(8), allocatable:: a(:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (allocated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_1d_da

!===============================================================================
!> Dealloccate a 2-D double array
!===============================================================================
SUBROUTINE deallocate_2d_da (a, label)
  real(8), allocatable:: a(:,:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (allocated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_2d_da

!===============================================================================
!> Dealloccate a 3-D double array
!===============================================================================
SUBROUTINE deallocate_3d_da (a, label)
  real(8), allocatable:: a(:,:,:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (allocated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_3d_da

!===============================================================================
!> Dealloccate a 4-D double array
!===============================================================================
SUBROUTINE deallocate_4d_da (a, label)
  real(8), allocatable:: a(:,:,:,:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (allocated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_4d_da

!===============================================================================
!> Dealloccate a 5-D double array
!===============================================================================
SUBROUTINE deallocate_5d_da (a, label)
  real(8), allocatable:: a(:,:,:,:,:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (allocated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_5d_da

!===============================================================================
!> Dealloccate a 1-D double array
!===============================================================================
SUBROUTINE deallocate_1d_ra (a, label)
  real(4), allocatable:: a(:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (allocated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_1d_ra

!===============================================================================
!> Dealloccate a 2-D double array
!===============================================================================
SUBROUTINE deallocate_2d_ra (a, label)
  real(4), allocatable:: a(:,:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (allocated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_2d_ra

!===============================================================================
!> Dealloccate a 3-D double array
!===============================================================================
SUBROUTINE deallocate_3d_ra (a, label)
  real(4), allocatable:: a(:,:,:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (allocated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_3d_ra

!===============================================================================
!> Dealloccate a 4-D double array
!===============================================================================
SUBROUTINE deallocate_4d_ra (a, label)
  real(4), allocatable:: a(:,:,:,:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (allocated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_4d_ra

!===============================================================================
!> Dealloccate a 5-D double array
!===============================================================================
SUBROUTINE deallocate_5d_ra (a, label)
  real(4), allocatable:: a(:,:,:,:,:)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  if (allocated(a)) then
    call mem_register (storage_size(a), -product(shape(a)), label)
    deallocate (a)
  end if
END SUBROUTINE deallocate_5d_ra

!===============================================================================
!> For verbose==0 no checking is done
!> For verbose==1 mem size is written as part of task%log lines
!> For verbose==2 or missing a log message is printed if memory has increased.
!> For verbose==3 or larger a log message is always printed
!===============================================================================
SUBROUTINE memory_info (self, label, verbose, action)
  class(os_t):: self
  character(*), optional:: label, action
  integer, optional:: verbose
  !.............................................................................
  real(8):: dgb
  real(8), save:: gb=0d0, gb_thread=0d0
  !$omp threadprivate(gb_thread)
  integer:: verbose_L
  !-----------------------------------------------------------------------------
  if (present(verbose)) then
    verbose_L = verbose
  else
    verbose_L = max(1,self%verbose)
  end if
  if (verbose_L > 0) then
    if (present(action)) then
      if (action == 'begin') then
        gb_thread = self%gb_thread(omp%thread)
      else if (action == 'end') then
        dgb = self%gb_thread(omp%thread) - gb_thread
        write (stdlog,'(a,1p,g12.3)') &
          'os_t%mem_info: thread,dgb=', omp%thread, dgb
      end if
    end if
    call                                           lock%set ('os_t%memory_info')
    if (gb > 0d0) then
      dgb = self%gb-gb
      if (verbose_L > 1 .or. dgb /= 0d0) then
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
    call                                         lock%unset ('os_t%memory_info')
  end if
END SUBROUTINE memory_info

!===============================================================================
!> Register increase of memory, counted in storage_size bits, times a count
!===============================================================================
SUBROUTINE mem_register (bits, count, label)
  integer:: bits, count
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  call os%gb_mem (bits, count, label)
  if (os%verbose > 1) then
    if (present(label)) then
      write(stdlog,1) bits,' bits per item', count,' items, for', &
        (bits/(8.*1024.**3))*count, os%gb_thread(omp%thread),' GB '//trim(label)
      1 format(" mem_reg: ",i6,a,i9,a,2f11.6,a)
    else
      write(stdlog,1) bits,' bits per item', count,' items, for', &
        (bits/(8.*1024.**3))*count, os%gb_thread(omp%thread), ' GB'
    end if
  end if
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

!===============================================================================
!> Detect inf and nan
!===============================================================================
LOGICAL FUNCTION isinf(a)
  use ieee_arithmetic, only : ieee_is_nan
  real(var_k):: a
  isinf = ieee_is_nan(a) .or. (abs(a) > huge(1.0_KindScalarVar))
END FUNCTION isinf
LOGICAL FUNCTION isnan4(a)
  use ieee_arithmetic, only : ieee_is_nan
  real(4):: a
  isnan4 = ieee_is_nan(a)
END FUNCTION isnan4
LOGICAL FUNCTION isnan8(a)
  use ieee_arithmetic, only : ieee_is_nan
  real(8):: a
  isnan8 = ieee_is_nan(a)
END FUNCTION isnan8

END MODULE os_mod
