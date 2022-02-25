!===============================================================================
!> Cooling module -- so far only for HLLS
!===============================================================================
MODULE cooling_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE patch_mod
  USE benchmark_mod
  USE scratch_mod
  USE eos_mod
  !
  implicit none
  private
  type, public:: cooling_t
  contains
    procedure:: pre_update
  end type
CONTAINS

!===============================================================================
!> Boundary conditions as specified in the benchmark
!===============================================================================
SUBROUTINE pre_update (self, patch)
  class(cooling_t):: self
  class(patch_t), target:: patch
  !.............................................................................
  real(kind=var_k), dimension(:,:,:), pointer, contiguous:: d, e, s, q
  real(8), pointer:: z(:)
  integer:: id
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('cooling_t%primary', itimer=itimer)
  id    =  patch%id
  !-----------------------------------------------------------------------------
  ! Get quantities for benchmark%cooling
  !-----------------------------------------------------------------------------
  allocate (z(patch%gn(3)), e(patch%gn(1),patch%gn(2),patch%gn(3)))
  z = patch%mesh(3)%r + patch%mesh(3)%p
  d => patch%mem(:,:,:,patch%idx%d,patch%it,1)
  if (patch%idx%s > 0) then
    s => patch%mem(:,:,:,patch%idx%s,patch%it,1)
    e = d**eos%gamma*exp(s*(eos%gamma-1d0)/d)/(eos%gamma-1d0)
  else
    e = patch%mem(:,:,:,patch%idx%e ,patch%it,1) - 0.5*( &
        patch%mem(:,:,:,patch%idx%px,patch%it,1)**2 + &
      + patch%mem(:,:,:,patch%idx%py,patch%it,1)**2 + &
      + patch%mem(:,:,:,patch%idx%pz,patch%it,1)**2)/d
  end if
  !-----------------------------------------------------------------------------
  ! Surface cooling, translated back to entropy per unit volume
  !-----------------------------------------------------------------------------
  call scratch%allocate (scratch%heating_per_unit_mass, patch%gn)
  q => scratch%heating_per_unit_mass
  call benchmark%cooling (z, d, e, q)
  deallocate (z, e)
  call trace%end (itimer)
END SUBROUTINE pre_update

END MODULE cooling_mod
