!===============================================================================
!> Boundary conditions for centered variables, which have the physical boundary
!> in-between mesh points lo and li, or/and in-between mesh points ui and uo.
!===============================================================================
MODULE boundaries_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE bits_mod
  USE mesh_mod
  USE index_mod
  USE standard_boundaries_mod
  USE eos_mod
  USE scaling_mod
  USE units_mod
  USE math_mod
  USE patch_mod
  USE units_mod
  ! ----------------------------------------------------------------------------
  ! Experiment - specific USE modules
  ! ----------------------------------------------------------------------------
  !USE lower_damp_mod
  implicit none
  private
  type, extends(standard_boundaries_t), public:: boundaries_t
  contains
    ! --------------------------------------------------------------------------
    ! Experiment - specific procedures
    ! --------------------------------------------------------------------------
    procedure:: mirror
  end type
  integer, save:: verbose=1
CONTAINS

!===============================================================================
SUBROUTINE mirror(self)
  class(boundaries_t):: self
  integer :: i, b, it
  ! ----------------------------------------------------------------------------
  if (self%patch%ng(3) == 0) return
  it = self%patch%it

  if (self%is_set(bits%zu)) then
      b = ubound(self%patch%mem,3)
      self%patch%mem(:,:,b,:,it,:) = self%patch%mem(:,:,b-1,:,it,:)
  endif

  if (self%is_set(bits%zl)) then
      b = lbound(self%patch%mem,3)
      self%patch%mem(:,:,b,:,it,:) = self%patch%mem(:,:,b+1,:,it,:)
  endif
END SUBROUTINE

END MODULE boundaries_mod
