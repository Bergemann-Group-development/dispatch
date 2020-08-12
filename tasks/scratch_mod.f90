!===============================================================================
!> Module to handle allocation of permanent thread-private memory
!===============================================================================
MODULE scratch_mod
  USE io_mod
  USE kinds_mod
  implicit none
  private
  type, public:: scratch_t
    real(kind=KindScalarVar), pointer, contiguous, dimension(:,:,:):: &
      heating_per_unit_mass => null(), heating_per_unit_volume => null()
    real(kind=KindScalarVar), pointer, contiguous, dimension(:,:,:,:):: &
      force_per_unit_mass => null(), force_per_unit_volume => null()
  contains
    procedure, nopass:: allocate_3d
    procedure, nopass:: allocate_4d
    generic:: allocate => allocate_3d, allocate_4d
  end type
  type(scratch_t), public:: scratch
  !$threadptivate (scratch)
CONTAINS

!===============================================================================
!> Allocate a 3-D real array
!===============================================================================
SUBROUTINE allocate_3d(a, m)
  real(kind=KindScalarVar), pointer, contiguous:: a(:,:,:)
  integer:: m(:)
  !-----------------------------------------------------------------------------
  if (.not.associated(a)) then
    allocate (a(m(1),m(2),m(3)))
    call io%bits_mem (storage_size(a), product(shape(a)), 'fpm')
  end if
END SUBROUTINE allocate_3d

!===============================================================================
!> Allocate a 4-D real array
!===============================================================================
SUBROUTINE allocate_4d(a, m)
  real(kind=KindScalarVar), pointer, contiguous:: a(:,:,:,:)
  integer:: m(:)
  !-----------------------------------------------------------------------------
  if (.not.associated(a)) then
    allocate (a(m(1),m(2),m(3),m(4)))
    call io%bits_mem (storage_size(a), product(shape(a)), 'fpm')
  end if
END SUBROUTINE allocate_4d

END MODULE scratch_mod
