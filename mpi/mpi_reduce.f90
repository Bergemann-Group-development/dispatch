!> =============================================================================
!> $Id$
!> MPI reduction operators
!> =============================================================================
MODULE mpi_reduction
  USE mpi_basic
  implicit none
  integer, private:: mpi_err
  include "mpif.h"
  type, extends(mpi_t):: mpi_reduction_t
  contains
    procedure min_int
    procedure sum_int
    procedure max_int
    procedure min_real
    procedure sum_real
    procedure max_real
    procedure max_real8
    procedure min_ints
    procedure max_ints
    procedure min_reals
    procedure max_reals
  end type
PUBLIC
CONTAINS

!> =============================================================================
!> Return the smallest integer among the MPI ranks
!> =============================================================================
SUBROUTINE min_int (mp, a)
  implicit none
  class(mpi_reduction_t):: mp
  integer:: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,1,MPI_INT,MPI_MIN,mp%comm,mpi_err)
  call mp%assert ('min_int_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the sum of integers over the MPI ranks
!> =============================================================================
SUBROUTINE sum_int (mp, a)
  implicit none
  class(mpi_reduction_t):: mp
  integer:: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,1,MPI_INT,MPI_SUM,mp%comm,mpi_err)
  call mp%assert ('sum_int_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the largest integer among the MPI ranks
!> =============================================================================
SUBROUTINE max_int (mp, a)
  implicit none
  class(mpi_reduction_t):: mp
  integer:: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,1,MPI_INT,MPI_MAX,mp%comm,mpi_err)
  call mp%assert ('max_int_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the smallest real among the MPI ranks
!> =============================================================================
SUBROUTINE min_real (mp, a)
  implicit none
  class(mpi_reduction_t):: mp
  real:: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,1,MPI_REAL,MPI_MIN,mp%comm,mpi_err)
  call mp%assert ('min_real_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the sum of reals over the MPI ranks
!> =============================================================================
SUBROUTINE sum_real (mp, a)
  implicit none
  class(mpi_reduction_t):: mp
  real:: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,1,MPI_REAL,MPI_SUM,mp%comm,mpi_err)
  call mp%assert ('sum_real_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the largest real among the MPI ranks
!> =============================================================================
SUBROUTINE max_real (mp, a)
  implicit none
  class(mpi_reduction_t):: mp
  real:: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,1,MPI_REAL,MPI_MAX,mp%comm,mpi_err)
  call mp%assert ('max_real_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the largest real(kind=8) among the MPI ranks
!> =============================================================================
SUBROUTINE max_real8 (mp, a)
  implicit none
  class(mpi_reduction_t):: mp
  real(kind=8):: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,1,MPI_REAL8,MPI_MAX,mp%comm,mpi_err)
  call mp%assert ('max_real_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the smallest integers among the MPI ranks
!> =============================================================================
SUBROUTINE min_ints (mp, a, n)
  implicit none
  class(mpi_reduction_t):: mp
  integer n
  integer, dimension(n):: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,n,MPI_INT,MPI_MIN,mp%comm,mpi_err)
  call mp%assert ('min_ints_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the largest integers among the MPI ranks
!> =============================================================================
SUBROUTINE max_ints (mp, a, n)
  implicit none
  class(mpi_reduction_t):: mp
  integer n
  integer, dimension(n):: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,n,MPI_INT,MPI_MAX,mp%comm,mpi_err)
  call mp%assert ('max_ints_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the smallest reals among the MPI ranks
!> =============================================================================
SUBROUTINE min_reals (mp, a, n)
  implicit none
  class(mpi_reduction_t):: mp
  integer n
  real, dimension(n):: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,n,MPI_REAL,MPI_MIN,mp%comm,mpi_err)
  call mp%assert ('min_reals_mpi', mpi_err)
  a = a1
END SUBROUTINE

!> =============================================================================
!> Return the largest reals among the MPI ranks
!> =============================================================================
SUBROUTINE max_reals (mp, a, n)
  implicit none
  class(mpi_reduction_t):: mp
  integer n
  real, dimension(n):: a, a1
!-----------------------------------------------------------------------
  call MPI_Allreduce(a,a1,n,MPI_REAL,MPI_MAX,mp%comm,mpi_err)
  call mp%assert ('max_reals_mpi', mpi_err)
  a = a1
END SUBROUTINE

END MODULE mpi_reduction
