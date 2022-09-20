!===============================================================================
!> Simple MPI global array data type.  The mpi_global_t%update() function
!> returns the current accumulated sum of a double precision array, whose value
!> is maintained with MPI_Get_accumulate.  Calling sequence:
!>
!>   type(mpi_global_t):: global
!>   ...
!>   call global%init (4)
!>   ...
!>   real(8):: accumulated_sum(4), contrib(4)
!>   ...
!>   contrib(:) = ...
!>   accumulated_sum = global%update (contrib)
!>
!> To maintain for example a set of global averages, use the local value of the
!> the contribution in the first update call, and continue with the contribution
!> set to the diffence between the current value and the previous one.
!===============================================================================
MODULE mpi_global_mod
  USE mpi_lib
  USE mpi_mod
  USE trace_mod
  USE io_mod
  implicit none
  private
  type, public:: mpi_global_t
    integer:: window_d1=0, window_i2=0
    integer:: n=0
    real(8), pointer:: a(:)    => null()
    integer, pointer:: ii(:,:) => null()
  contains
    procedure:: init_d1
    procedure:: init_i2
    generic:: init => init_d1, init_i2
    procedure, private:: update4
    procedure, private:: update8
    generic:: update => update4, update8
    procedure:: put_d1
    procedure:: put_i2
    generic:: put => put_d1, put_i2
    procedure:: get_d1
    procedure:: get_i2
    generic:: get => get_d1, get_i2
  end type
CONTAINS

!===============================================================================
!> Create and initialize a global double precision array(n).
!> NOTE: This is a COLLECTIVE call!
!===============================================================================
SUBROUTINE init_d1 (self, n, d)
  class(mpi_global_t):: self
  integer:: n
  real(8), optional:: d(:)
  !.............................................................................
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND):: nbytes
#endif
  integer:: err
  !----------------------------------------------------------------------------
  call trace%begin('mpi_global_t%init8')
  if (self%window_d1==0) then
    self%n = n
    allocate (self%a(n))
    if (present(d)) then
      self%a = d
    else
      self%a = 0d0
    end if
#ifdef MPI
    nbytes = 8*n
    !$omp single
    call MPI_Win_create (self%a, nbytes, 8, MPI_INFO_NULL, MPI_COMM_WORLD, &
                         self%window_d1, err)
    call mpi%assert ('MPI_Win_create', err)
    !$omp end single nowait
#else
    self%a = 0d0
#endif
  end if
  call trace%end()
END SUBROUTINE init_d1

!===============================================================================
!> Create and initialize a global double precision array(n).
!> NOTE: This is a COLLECTIVE call!
!===============================================================================
SUBROUTINE init_i2 (self, ii)
  class(mpi_global_t):: self
  integer:: ii(:,:)
  !.............................................................................
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND):: nbytes
#endif
  integer:: err
  !----------------------------------------------------------------------------
  call trace%begin('mpi_global_t%init8')
  if (self%window_i2==0) then
    allocate (self%ii(size(ii,1),size(ii,2)))
    self%ii = ii
#ifdef MPI
    nbytes = 4*product(shape(ii))
    !$omp single
    call MPI_Win_create (self%ii, nbytes, 4, MPI_INFO_NULL, MPI_COMM_WORLD, &
                         self%window_i2, err)
    call mpi%assert ('MPI_Win_create', err)
    !$omp end single nowait
#else
    self%ii = 0d0
#endif
  end if
  call trace%end()
END SUBROUTINE init_i2

!===============================================================================
!> Add a value a(:) to the counter on master.
!===============================================================================
FUNCTION update4 (self, a) RESULT (b)
  class(mpi_global_t):: self
  real:: a(self%n), b(self%n)
  !.............................................................................
  b= update8 (self, real(a,kind=8))
END FUNCTION update4
!===============================================================================
FUNCTION update8 (self, a) RESULT (b)
  class(mpi_global_t):: self
  real(8):: a(self%n), b(self%n)
  integer:: i
  !.............................................................................
  call trace%begin('mpi_global_t%update8')
#ifdef MPI
  if (mpi%size == 1) then
    do i=1,self%n
      !$omp atomic update
      self%a(i) = self%a(i) + a(i)
    end do
    b = self%a
  else
    b = update8_real (self, a)
  end if
#else
  do i=1,self%n
    !$omp atomic update
    self%a(i) = self%a(i) + a(i)
  end do
  b = self%a
#endif
  call trace%end()
END FUNCTION update8

!===============================================================================
FUNCTION update8_real (self, a) RESULT (b)
  class(mpi_global_t):: self
  real(8):: a(self%n), b(self%n)
  !.............................................................................
  integer:: master=0, err
  !-----------------------------------------------------------------------------
#ifdef MPI
  if (mpi%mode == MPI_THREAD_MULTIPLE) then
    call update_l
  else
    call update_l
  end if
#endif
  b = b + a
contains

!===============================================================================
  subroutine update_l
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND):: offset=0
  integer, save:: itimer=0
  call trace%begin ('mpi_global%update', itimer=itimer)
  !$omp critical (global_cr)
  call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, master, 0, self%window_d1, err)
  call mpi%assert ('MPI_Win_lock', err)
  call MPI_Get_accumulate (a, self%n, MPI_REAL8, &
                           b, self%n, MPI_REAL8, master, offset, &
                              self%n, MPI_REAL8, MPI_SUM, self%window_d1, err)
  call mpi%assert ('MPI_Accumulate', err)
  call MPI_Win_unlock (master, self%window_d1, err)
  !$omp end critical (global_cr)
  call mpi%assert ('MPI_Win_unlock', err)
  call trace%end (itimer)
#endif
  end subroutine update_l
END FUNCTION update8_real

!===============================================================================
SUBROUTINE put_d1 (self, d)
  class(mpi_global_t):: self
  real(8):: d(:)
  !.............................................................................
  integer:: master=0, err, n
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND):: offset=0
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_global%put_i2_real', itimer=itimer)
  !$omp critical (global_cr)
  n = 8*size(d)
  call MPI_Put (d, n, MPI_REAL8, master, offset, n, MPI_REAL8, &
                self%window_d1, err)
  call mpi%assert ('MPI_Put', err)
  call MPI_Win_unlock (master, self%window_d1, err)
  !$omp end critical (global_cr)
  call mpi%assert ('MPI_Win_unlock', err)
  call trace%end (itimer)
#endif
END SUBROUTINE put_d1 

!===============================================================================
SUBROUTINE get_d1 (self, d)
  class(mpi_global_t):: self
  real(8):: d(:)
  !.............................................................................
  integer:: master=0, err, n
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND):: offset=0
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_global%put_d1_real', itimer=itimer)
  !$omp critical (global_cr)
  call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, master, 0, self%window_d1, err)
  call mpi%assert ('MPI_Win_lock', err)
  n = 8*size(d)
  call MPI_Get (d, n, MPI_REAL8, master, offset, n, MPI_REAL8, &
                self%window_d1, err)
  call mpi%assert ('MPI_Get', err)
  !$omp end critical (global_cr)
  call trace%end (itimer)
#endif
END SUBROUTINE get_d1 

!===============================================================================
SUBROUTINE put_i2 (self, ii)
  class(mpi_global_t):: self
  integer:: ii(:,:)
  !.............................................................................
  integer:: master=0, err, n, win
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND):: offset=0
  !-----------------------------------------------------------------------------
  !$omp critical (global_cr)
  n = product(shape(self%ii))
  win = self%window_i2
  call MPI_Put (ii, n, MPI_INTEGER, master, offset, n, MPI_INTEGER, win, err)
  call mpi%assert ('MPI_Put', err)
  call MPI_Win_unlock (master, win, err)
  !$omp end critical (global_cr)
  call mpi%assert ('MPI_Win_unlock', err)
#endif
END SUBROUTINE put_i2 

!===============================================================================
SUBROUTINE get_i2 (self, ii)
  class(mpi_global_t):: self
  integer:: ii(:,:)
  !.............................................................................
  integer:: master=0, err, n, win
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND):: offset=0
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_global%get_i2')
  !$omp critical (global_cr)
  win = self%window_i2
  call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, master, 0, win, err)
  call mpi%assert ('MPI_Win_lock', err)
  n = product(shape(self%ii))
  call MPI_Get (ii, n, MPI_INTEGER, master, offset, n, MPI_INTEGER, win, err)
  call mpi%assert ('MPI_Get', err)
  !$omp end critical (global_cr)
  call trace%end ()
#endif
END SUBROUTINE get_i2 

END MODULE mpi_global_mod
