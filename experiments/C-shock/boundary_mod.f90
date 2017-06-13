!> =============================================================================
!> Boundary conditions for a planet surface
!> =============================================================================
MODULE boundary_mod
  USE io_mod
  USE trace_mod
  USE bits_mod
  USE mesh_mod
  USE index_mod
  implicit none
  private
  type, public:: boundary_t
  contains
    procedure:: xcondition
  end type
  type(boundary_t), public:: boundary
CONTAINS

!> =============================================================================
!> Impose boundary conditions in the x-direction, taking staggering into account
!> Zero-gradient conditions are used.
!> =============================================================================
SUBROUTINE xcondition (self, mem, m, lower)
  class(boundary_t):: self
  type(index_t):: idx
  real, dimension(:,:,:):: mem
  class(mesh_t):: m
  logical, optional:: lower
  !.............................................................................
  integer:: ix
  !-----------------------------------------------------------------------------
  call trace_begin('boundary_t%xcondition')

  if (present(lower)) then
    !do ix=m%lb,m%lo
      mem(m%lb,:,:) = mem(m%li,:,:)
    !end do
  else
    !do ix=m%uo,m%ub
      mem(m%ub,:,:) = mem(m%ui,:,:)
    !end do
  end if

  call trace_end
END SUBROUTINE xcondition

END MODULE boundary_mod