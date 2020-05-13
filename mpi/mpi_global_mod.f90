!===============================================================================
!> Simple MPI global array data type
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
#ifdef MPI
    nbytes = 8*n
#endif
    self%a = 0d0
#ifdef __GFORTRAN__
    self%window = -1
#else
#ifdef MPI
  if (mpi%mode == MPI_THREAD_MULTIPLE) then
   call MPI_Win_create (self%a, nbytes, 8, MPI_INFO_NULL, MPI_COMM_WORLD, &
                        self%window, mpi_err)
  else
    !$omp critical (mpi_cr)
    call MPI_Win_create (self%a, nbytes, 8, MPI_INFO_NULL, MPI_COMM_WORLD, &
                         self%window, mpi_err)
    !$omp end critical (mpi_cr)
  end if
#endif
    call mpi%assert ('MPI_Win_create', mpi_err)
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
#ifdef __GFORTRAN__
  do i=1,3
    !$omp atomic update
    self%a(i) = self%a(i) + a(i)
  end do
  b = self%a
#else
  if (mpi%size == 1) then
    do i=1,3
      !$omp atomic update
      self%a(i) = self%a(i) + a(i)
    end do
    b = self%a
  else
    b = update8_real (self, a)
  end if
#endif
END FUNCTION update8

!===============================================================================
FUNCTION update8_real (self, a) RESULT (b)
  class(mpi_global_t):: self
  real(8):: a(self%n), b(self%n)
  !.............................................................................
  integer:: master=0, mpi_err
  !-----------------------------------------------------------------------------
  call trace%begin('mpi_global_t%init')
#ifdef MPI
  if (mpi%mode == MPI_THREAD_MULTIPLE) then
    call update_l
  else
    !$omp critical (mpi_cr)
    call update_l
    !$omp end critical (mpi_cr)
  end if
#endif
  call trace%end()
contains

!===============================================================================
  subroutine update_l
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND):: offset=0
    call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, master, 0, self%window, mpi_err)
    call MPI_Get_accumulate (a, self%n, MPI_REAL8, &
                        self%a, self%n, MPI_REAL8, master, offset, &
                                self%n, MPI_REAL8, MPI_SUM, self%window, mpi_err)
    call mpi%assert ('MPI_Accumulate', mpi_err)
    call MPI_Win_flush (master, self%window, mpi_err)   
    if (mpi%rank == master) then
      b = self%a
    else
      b = self%a + a
    end if
    call MPI_Win_unlock (master, self%window, mpi_err)
#endif
  end subroutine update_l
END FUNCTION update8_real

END MODULE mpi_global_mod
