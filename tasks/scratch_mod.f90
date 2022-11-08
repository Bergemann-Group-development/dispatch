!===============================================================================
!> Module to handle allocation of permanent thread-private memory
!===============================================================================
MODULE scratch_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  implicit none
  private
  type, public:: scratch_t
    real(kind=KindScalarVar), pointer, contiguous, dimension(:,:,:):: &
      heating_per_unit_mass => null(), heating_per_unit_volume => null()
    real(kind=KindScalarVar), pointer, contiguous, dimension(:,:,:,:):: &
      force_per_unit_mass => null(), force_per_unit_volume => null()
  contains
    procedure:: pre_update
    !procedure, nopass:: allocate_3d
    !procedure, nopass:: allocate_4d
    !generic:: allocate => allocate_3d, allocate_4d
  end type
  type(scratch_t), public:: scratch
  !$omp threadprivate (scratch)
CONTAINS

!===============================================================================
!> Allocate a 3-D real array
!===============================================================================
SUBROUTINE pre_update (self)
  class(scratch_t):: self
  if (associated(self%heating_per_unit_volume)) self%heating_per_unit_volume=0.0
  if (associated(self%heating_per_unit_mass  )) self%heating_per_unit_mass=0.0
  if (associated(self%force_per_unit_volume)) self%force_per_unit_volume=0.0
  if (associated(self%force_per_unit_mass  )) self%force_per_unit_mass=0.0
END SUBROUTINE pre_update

!!===============================================================================
!!> Allocate a 3-D real array
!!===============================================================================
!SUBROUTINE allocate_3d(a, m)
!  real(kind=KindScalarVar), pointer, contiguous:: a(:,:,:)
!  integer:: m(:)
!  !-----------------------------------------------------------------------------
!  if (.not.associated(a)) then
!    allocate (a(m(1),m(2),m(3)))
!    call io%bits_mem (storage_size(a), product(shape(a)), 'fpm')
!    a = 0.0
!  end if
!END SUBROUTINE allocate_3d
!
!!===============================================================================
!!> Allocate a 4-D real array
!!===============================================================================
!SUBROUTINE allocate_4d(a, m)
!  real(kind=KindScalarVar), pointer, contiguous:: a(:,:,:,:)
!  integer:: m(:)
!  !-----------------------------------------------------------------------------
!  call trace%begin ('scratch_t%allocate_4d')
!  if (.not.associated(a)) then
!    allocate (a(m(1),m(2),m(3),m(4)))
!    call io%bits_mem (storage_size(a), product(shape(a)), 'fpm')
!    a = 0.0
!  end if
!  call trace%end ()
!END SUBROUTINE allocate_4d

END MODULE scratch_mod
