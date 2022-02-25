!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE vector_mod
  USE io_mod
  USE trace_mod
  USE stagger_mod
  USE scalar_mod
  USE mesh_mod
  implicit none
  private
  type, public:: vector_t
    class(mesh_t), dimension(:), pointer:: mesh => null()
    type(scalar_t):: scalar
  contains
    procedure:: init
    procedure:: grad
    procedure:: div
    procedure:: center
    procedure:: face
    procedure:: stats
  end type
  type(stagger_t):: stagger
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self, mesh)
  class(vector_t):: self
  class(mesh_t), dimension(:), pointer:: mesh
  !.............................................................................
  call trace_begin ('vector_t%init')
  self%mesh => mesh
  self%scalar%mesh => mesh
  call trace_end
END SUBROUTINE init
!===============================================================================

!===============================================================================
!> Return gradient of scalar
!===============================================================================
FUNCTION grad (self, s) RESULT (v)
  class(vector_t):: self
  real, dimension(:,:,:):: s
  real, dimension(size(s,1),size(s,2),size(s,3),3):: v
  real:: out
  !.............................................................................
  call trace_begin ('vector_t%grad')
  v(:,:,:,1) = stagger%ddxdn (self%mesh, s)
  v(:,:,:,2) = stagger%ddydn (self%mesh, s)
  v(:,:,:,3) = stagger%ddzdn (self%mesh, s)
  call trace_end
END FUNCTION grad
!===============================================================================

!===============================================================================
!> Return gradient of scalar
!===============================================================================
FUNCTION div (self, t) RESULT (v)
  class(vector_t):: self
  real, dimension(:,:,:,:,:):: t
  real, dimension(size(t,1),size(t,2),size(t,3),3):: v
  real:: out
  !.............................................................................
  call trace_begin ('vector_t%div')
  v(:,:,:,1) = stagger%ddxdn (self%mesh, t(:,:,:,1,1)) &
             + stagger%ddyup (self%mesh, t(:,:,:,1,2)) &
             + stagger%ddzup (self%mesh, t(:,:,:,1,3))
  v(:,:,:,2) = stagger%ddxup (self%mesh, t(:,:,:,2,1)) &
             + stagger%ddydn (self%mesh, t(:,:,:,2,2)) &
             + stagger%ddzup (self%mesh, t(:,:,:,2,3))
  v(:,:,:,3) = stagger%ddxup (self%mesh, t(:,:,:,3,1)) &
             + stagger%ddyup (self%mesh, t(:,:,:,3,2)) &
             + stagger%ddzdn (self%mesh, t(:,:,:,3,3))
  call trace_end
END FUNCTION div
!===============================================================================

!===============================================================================
!> Return gradient of scalar
!===============================================================================
FUNCTION center (self, s) RESULT (v)
  class(vector_t):: self
  real, dimension(:,:,:):: s
  real, dimension(size(s,1),size(s,2),size(s,3),3):: v
  real:: out
  !.............................................................................
  call trace_begin ('vector_t%center')
  v(:,:,:,1) = stagger%xup (s)
  v(:,:,:,2) = stagger%yup (s)
  v(:,:,:,3) = stagger%zup (s)
  call trace_end
END FUNCTION center
!===============================================================================

!===============================================================================
!> Return gradient of scalar
!===============================================================================
FUNCTION face (self, s) RESULT (v)
  class(vector_t):: self
  real, dimension(:,:,:):: s
  real, dimension(size(s,1),size(s,2),size(s,3),3):: v
  real:: out
  !.............................................................................
  call trace_begin ('vector_t%face')
  v(:,:,:,1) = stagger%xdn (s)
  v(:,:,:,2) = stagger%ydn (s)
  v(:,:,:,3) = stagger%zdn (s)
  call trace_end
END FUNCTION face
!===============================================================================

!===============================================================================
!> Return gradient of scalar
!===============================================================================
SUBROUTINE stats (self, v, label)
  class(vector_t):: self
  real, dimension(:,:,:,:):: v
  character(len=*), optional:: label
  character(len=64):: l
  real:: out
  !.............................................................................
  call trace_begin ('vector_t%stats')
  l = trim(label)//'_x'; call self%scalar%stats (v(:,:,:,1), l)
  l = trim(label)//'_y'; call self%scalar%stats (v(:,:,:,2), l)
  l = trim(label)//'_z'; call self%scalar%stats (v(:,:,:,3), l)
  call trace_end
END SUBROUTINE stats
!===============================================================================

END MODULE vector_mod
