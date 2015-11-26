!*******************************************************************************
! Module with skeleton data type for MHD.
!*******************************************************************************
MODULE mhd_mod
  USE mesh,    only: m
  USE vectors
  USE io
  type:: mhd_t
    logical:: is_allocated=.false.
    real(8):: t                               ! simulated time
    real(8):: offset(3)=0.0                   ! offset relative cell center
    real, dimension(:,:,:), pointer:: d       ! mass per unit volume
    real, dimension(:,:,:), pointer:: s       ! entropy per unit volume
    type(vector):: p                          ! momentum per unit volume
    type(vector):: B                          ! magnetic flux per unit area
  end type
CONTAINS

!*******************************************************************************
SUBROUTINE allocate_mhd (a)
  implicit none
  type(mhd_t):: a
!...............................................................................
  if (a%is_allocated) return
  call trace_begin ('allocate_mhd')
  call allocate_scalar (a%d)
  call allocate_scalar (a%s)
  call allocate_vector (a%p)
  call allocate_vector (a%B)
  a%is_allocated = .true.
  call trace_end
END SUBROUTINE allocate_mhd

!*******************************************************************************
SUBROUTINE deallocate_mhd (a)
  implicit none
  type(mhd_t):: a
!...............................................................................
  if (.not.a%is_allocated) return
  call deallocate_scalar (a%d)
  call deallocate_scalar (a%s)
  call deallocate_vector (a%p)
  call deallocate_vector (a%B)
  a%is_allocated = .false.
END SUBROUTINE deallocate_mhd
END MODULE mhd_mod
