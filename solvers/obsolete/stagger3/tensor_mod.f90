!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE tensor_mod
  USE io_mod
  USE trace_mod
  USE stagger_mod
  USE mesh_mod
  implicit none
  private
  type, public:: tensor_t
    class(mesh_t), dimension(:), pointer:: mesh => null()
  contains
    procedure:: init
    procedure:: grad
    procedure:: strain
    procedure:: reynolds
    procedure:: add_diagonal
    procedure:: edge_scalar
    procedure:: edge_vector
    generic:: edge => edge_scalar, edge_vector
  end type
  type(stagger_t):: stagger
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self, mesh)
  class(tensor_t):: self
  class(mesh_t), dimension(:), pointer:: mesh
  !.............................................................................
  call trace_begin ('tensor_t%init')
  self%mesh => mesh
  call trace_end
END SUBROUTINE init
!===============================================================================

!===============================================================================
!> Return interpolated scalar, with diagonal cell-centered, and off-diagnonals
!> edge-centered.  The 4th index is the vector index, the 5th the shift index.
!===============================================================================
FUNCTION edge_scalar (self, s, log_s) RESULT (t)
  class(tensor_t):: self
  real, dimension(:,:,:):: s
  real, dimension(:,:,:), optional:: log_s
  real, dimension(:,:,:), pointer:: tmp_s
  real, dimension(size(s,1),size(s,2),size(s,3),3,3):: t
  !.............................................................................
  call trace_begin ('tensor_t%edge_scalar')
  t(:,:,:,1,1) =              s(:,:,:)
  t(:,:,:,2,2) =              s(:,:,:)
  t(:,:,:,3,3) =              s(:,:,:)
  if (present(log_s)) then
    t(:,:,:,1,2) = exp(stagger%xdn(stagger%ydn(log_s)))
    t(:,:,:,2,3) = exp(stagger%ydn(stagger%zdn(log_s)))
    t(:,:,:,1,3) = exp(stagger%zdn(stagger%xdn(log_s)))
  else
    if (minval(s) > 0.0) then
      allocate (tmp_s(size(s,1),size(s,2),size(s,3)))
      tmp_s = log(s)
      t(:,:,:,1,2) = exp(stagger%xdn(stagger%ydn(tmp_s)))
      t(:,:,:,2,3) = exp(stagger%ydn(stagger%zdn(tmp_s)))
      t(:,:,:,1,3) = exp(stagger%zdn(stagger%xdn(tmp_s)))
      deallocate (tmp_s)
    else
      t(:,:,:,1,2) = stagger%xdn(stagger%ydn(s))
      t(:,:,:,2,3) = stagger%ydn(stagger%zdn(s))
      t(:,:,:,1,3) = stagger%zdn(stagger%xdn(s))
    end if
  end if
  t(:,:,:,2,1) = t(:,:,:,1,2)
  t(:,:,:,3,2) = t(:,:,:,2,3)
  t(:,:,:,3,1) = t(:,:,:,1,3)
  call trace_end
END FUNCTION edge_scalar
!===============================================================================

!===============================================================================
!> Return interpolated vector, with diagonal cell-centered, and off-diagnonals
!> edge-centered.  The 4th index is the vector index, the 5th the shift index.
!===============================================================================
FUNCTION edge_vector (self, v) RESULT (t)
  class(tensor_t):: self
  real, dimension(:,:,:,:):: v
  real, dimension(size(v,1),size(v,2),size(v,3),3,3):: t
  real:: out
  !.............................................................................
  call trace_begin ('tensor_t%edge_vector')
  t(:,:,:,1,1) = stagger%xup (v(:,:,:,1))
  t(:,:,:,2,2) = stagger%yup (v(:,:,:,2))
  t(:,:,:,3,3) = stagger%zup (v(:,:,:,3))
  t(:,:,:,1,2) = stagger%ydn (v(:,:,:,1))
  t(:,:,:,2,3) = stagger%zdn (v(:,:,:,2))
  t(:,:,:,1,3) = stagger%zdn (v(:,:,:,1))
  t(:,:,:,2,1) = stagger%xdn (v(:,:,:,2))
  t(:,:,:,3,2) = stagger%ydn (v(:,:,:,3))
  t(:,:,:,3,1) = stagger%xdn (v(:,:,:,3))
  call trace_end
END FUNCTION edge_vector

!===============================================================================
!> Return interpolated vector, with diagonal cell-centered, and off-diagnonals
!> edge-centered.  The 4th index is the vector index, the 5th the shift index.
!===============================================================================
FUNCTION reynolds (self, v) RESULT (t)
  class(tensor_t):: self
  real, dimension(:,:,:,:):: v
  real, dimension(size(v,1),size(v,2),size(v,3),3,3):: t
  real:: out
  !.............................................................................
  call trace_begin ('tensor_t%reynolds')
  t = 0.0
  t(:,:,:,1,1) = stagger%xup(v(:,:,:,1))**2
  t(:,:,:,2,2) = stagger%yup(v(:,:,:,2))**2
  t(:,:,:,3,3) = stagger%zup(v(:,:,:,3))**2
  t(:,:,:,1,2) = stagger%ydn(v(:,:,:,1))*stagger%xdn(v(:,:,:,2))
  t(:,:,:,2,3) = stagger%zdn(v(:,:,:,2))*stagger%ydn(v(:,:,:,3))
  t(:,:,:,1,3) = stagger%xdn(v(:,:,:,3))*stagger%zdn(v(:,:,:,1))
  t(:,:,:,2,1) = t(:,:,:,1,2)
  t(:,:,:,3,2) = t(:,:,:,2,3)
  t(:,:,:,3,1) = t(:,:,:,1,3)
  call trace_end
END FUNCTION reynolds

!===============================================================================
!> Return gradient of vector, with diagonal cell-centered, and off-diagnonals
!> edge-centered
!===============================================================================
FUNCTION grad (self, v) RESULT (t)
  class(tensor_t):: self
  real, dimension(:,:,:,:):: v
  real, dimension(size(v,1),size(v,2),size(v,3),3,3):: t
  real:: out
  !.............................................................................
  call trace_begin ('tensor_t%grad')
  t(:,:,:,1,1) = stagger%ddxup (self%mesh, v(:,:,:,1))
  t(:,:,:,2,2) = stagger%ddyup (self%mesh, v(:,:,:,2))
  t(:,:,:,3,3) = stagger%ddzup (self%mesh, v(:,:,:,3))
  t(:,:,:,1,2) = stagger%ddydn (self%mesh, v(:,:,:,1))
  t(:,:,:,2,3) = stagger%ddzdn (self%mesh, v(:,:,:,2))
  t(:,:,:,1,3) = stagger%ddzdn (self%mesh, v(:,:,:,1))
  t(:,:,:,2,1) = stagger%ddxdn (self%mesh, v(:,:,:,2))
  t(:,:,:,3,2) = stagger%ddydn (self%mesh, v(:,:,:,3))
  t(:,:,:,3,1) = stagger%ddxdn (self%mesh, v(:,:,:,3))
  call trace_end
END FUNCTION grad

!===============================================================================
!> Return strain tensor
!===============================================================================
FUNCTION strain (self, v) RESULT (t)
  class(tensor_t):: self
  real, dimension(:,:,:,:):: v
  real, dimension(size(v,1),size(v,2),size(v,3),3,3):: t
  real:: out
  !.............................................................................
  call trace_begin ('tensor_t%strain')
  t(:,:,:,1,1) = stagger%ddxup (self%mesh, v(:,:,:,1))
  t(:,:,:,2,2) = stagger%ddyup (self%mesh, v(:,:,:,2))
  t(:,:,:,3,3) = stagger%ddzup (self%mesh, v(:,:,:,3))
  t(:,:,:,1,2) = 0.5*(stagger%ddxdn (self%mesh, v(:,:,:,2)) + stagger%ddydn (self%mesh, v(:,:,:,1)))
  t(:,:,:,2,3) = 0.5*(stagger%ddydn (self%mesh, v(:,:,:,3)) + stagger%ddzdn (self%mesh, v(:,:,:,2)))
  t(:,:,:,1,3) = 0.5*(stagger%ddzdn (self%mesh, v(:,:,:,1)) + stagger%ddxdn (self%mesh, v(:,:,:,3)))
  t(:,:,:,2,1) = t(:,:,:,1,2)
  t(:,:,:,3,2) = t(:,:,:,2,3)
  t(:,:,:,3,1) = t(:,:,:,1,3)
  call trace_end
END FUNCTION strain

!===============================================================================
!> Subtract a diagonal part to a tensor
!===============================================================================
SUBROUTINE add_diagonal (self, out, in)
  class(tensor_t):: self
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3,3):: out
  !.............................................................................
  out(:,:,:,1,1) = out(:,:,:,1,1) + in
  out(:,:,:,2,2) = out(:,:,:,2,2) + in
  out(:,:,:,3,3) = out(:,:,:,3,3) + in
END SUBROUTINE add_diagonal

END MODULE tensor_mod
