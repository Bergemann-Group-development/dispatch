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
!===============================================================================
MODULE mpi_global_mod
  USE mpi_lib
  USE mpi_mod
  USE trace_mod
  implicit none
  private
  type, public:: mpi_global_t
    integer:: window=0
    real(8), pointer:: a(:)
    integer:: n=0
  contains
    procedure:: init
    procedure, private:: update4
    procedure, private:: update8
    generic:: update => update4, update8
  end type
CONTAINS

!===============================================================================
!> Create and initialize a global double precision array(n).
!> NOTE: This is a COLLECTIVE call!
!===============================================================================
SUBROUTINE init (self, n)
  class(mpi_global_t):: self
  integer:: n
  !.............................................................................
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND):: nbytes
#endif
  integer:: mpi_err
  !----------------------------------------------------------------------------
  call trace%begin('mpi_global_t%init8')
  if (self%window==0) then
    self%n = n
    allocate (self%a(n))
    self%a(:) = 0d0
#ifdef MPI
    nbytes = 8*n
    if (mpi%mode == MPI_THREAD_MULTIPLE) then
      call MPI_Win_create (self%a, nbytes, 8, MPI_INFO_NULL, MPI_COMM_WORLD, &
                          self%window, mpi_err)
    else
      !$omp critical (mpi_cr)
      call MPI_Win_create (self%a, nbytes, 8, MPI_INFO_NULL, MPI_COMM_WORLD, &
                          self%window, mpi_err)
      !$omp end critical (mpi_cr)
    end if
    call mpi%assert ('MPI_Win_create', mpi_err)
#else
    self%a = 0d0
#endif
  end if
  call trace%end()
END SUBROUTINE init

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
  integer:: master=0, mpi_err
  !-----------------------------------------------------------------------------
#ifdef MPI
  if (mpi%mode == MPI_THREAD_MULTIPLE) then
    call update_l
  else
    !$omp critical (mpi_cr)
    call update_l
    !$omp end critical (mpi_cr)
  end if
#endif
  b = b + a
contains

!===============================================================================
  subroutine update_l
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND):: offset=0
  call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, master, 0, self%window, mpi_err)
  call mpi%assert ('MPI_Win_lock', mpi_err)
  call MPI_Get_accumulate (a, self%n, MPI_REAL8, &
                           b, self%n, MPI_REAL8, master, offset, &
                              self%n, MPI_REAL8, MPI_SUM, self%window, mpi_err)
  call mpi%assert ('MPI_Accumulate', mpi_err)
  call MPI_Win_unlock (master, self%window, mpi_err)
  call mpi%assert ('MPI_Win_unlockx½', mpi_err)
#endif
  end subroutine update_l
END FUNCTION update8_real

END MODULE mpi_global_mod
