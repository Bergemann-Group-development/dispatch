!> =============================================================================
!> $Id$
!> MPI reduction operators
!> =============================================================================
MODULE mpi_reduce
  USE mpi_basic
  implicit none
PUBLIC
CONTAINS

!> =============================================================================
!> Return the smallest integer among the MPI ranks
!> =============================================================================
SUBROUTINE min_int_mpi (a)
  implicit none
  integer:: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,1,MPI_INT,MPI_MIN,mpi_comm_world,mpi_err)
  call assert_mpi ('min_int_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the sum of integers over the MPI ranks
!> =============================================================================
SUBROUTINE sum_int_mpi (a)
  implicit none
  integer:: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,1,MPI_INT,MPI_SUM,mpi_comm_world,mpi_err)
  call assert_mpi ('sum_int_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the largest integer among the MPI ranks
!> =============================================================================
SUBROUTINE max_int_mpi (a)
  implicit none
  integer:: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,1,MPI_INT,MPI_MAX,mpi_comm_world,mpi_err)
  call assert_mpi ('max_int_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the smallest real among the MPI ranks
!> =============================================================================
SUBROUTINE min_real_mpi (a)
  implicit none
  real:: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,1,MPI_REAL,MPI_MIN,mpi_comm_world,mpi_err)
  call assert_mpi ('min_real_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the sum of reals over the MPI ranks
!> =============================================================================
SUBROUTINE sum_real_mpi (a)
  implicit none
  real:: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,1,MPI_REAL,MPI_SUM,mpi_comm_world,mpi_err)
  call assert_mpi ('sum_real_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the largest real among the MPI ranks
!> =============================================================================
SUBROUTINE max_real_mpi (a)
  implicit none
  real:: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,1,MPI_REAL,MPI_MAX,mpi_comm_world,mpi_err)
  call assert_mpi ('max_real_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the largest real(kind=8) among the MPI ranks
!> =============================================================================
SUBROUTINE max_real8_mpi (a)
  implicit none
  real(kind=8):: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,1,MPI_REAL8,MPI_MAX,mpi_comm_world,mpi_err)
  call assert_mpi ('max_real_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the smallest integers among the MPI ranks
!> =============================================================================
SUBROUTINE min_ints_mpi (a, n)
  implicit none
  integer n
  integer, dimension(n):: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,n,MPI_INT,MPI_MIN,mpi_comm_world,mpi_err)
  call assert_mpi ('min_ints_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the largest integers among the MPI ranks
!> =============================================================================
SUBROUTINE max_ints_mpi (a, n)
  implicit none
  integer n
  integer, dimension(n):: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,n,MPI_INT,MPI_MAX,mpi_comm_world,mpi_err)
  call assert_mpi ('max_ints_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the smallest reals among the MPI ranks
!> =============================================================================
SUBROUTINE min_reals_mpi (a, n)
  implicit none
  integer n
  real, dimension(n):: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,n,MPI_REAL,MPI_MIN,mpi_comm_world,mpi_err)
  call assert_mpi ('min_reals_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the largest reals among the MPI ranks
!> =============================================================================
SUBROUTINE max_reals_mpi (a, n)
  implicit none
  integer n
  real, dimension(n):: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,n,MPI_REAL,MPI_MAX,mpi_comm_world,mpi_err)
  call assert_mpi ('max_reals_mpi', mpi_err)
  a = a1
END SUBROUTINE

END MODULE mpi_reduce
