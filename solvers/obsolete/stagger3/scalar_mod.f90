!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE scalar_mod
  USE io_mod
  USE trace_mod
  USE stagger_mod
  USE mesh_mod
  implicit none
  private
  type, public:: scalar_t
    class(mesh_t), dimension(:), pointer:: mesh => null()
  contains
    procedure:: init
    procedure:: div
    procedure:: stats
    procedure:: aver
    procedure:: minval => minval_scalar
    procedure:: maxval => maxval_scalar
    procedure:: magnitude_scalar
    procedure:: magnitude_vector
    procedure:: magnitude_tensor
    generic:: magnitude => magnitude_scalar, magnitude_vector, magnitude_tensor
  end type
CONTAINS

!===============================================================================
SUBROUTINE init (self, mesh)
  class(scalar_t):: self
  class(mesh_t), dimension(:), pointer:: mesh
  !.............................................................................
  call trace_begin ('scalar_t%init')
  self%mesh => mesh
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Return divergence of vector
!===============================================================================
FUNCTION div (self, v) RESULT (s)
  class(scalar_t):: self
  real, dimension(:,:,:,:):: v
  real, dimension(size(v,1),size(v,2),size(v,3)):: s
  type(stagger_t):: stagger
  !.............................................................................
  call trace_begin ('scalar_t%div')
  stagger%mesh => self%mesh
  s = stagger%ddxup (self%mesh, v(:,:,:,1)) &
    + stagger%ddyup (self%mesh, v(:,:,:,2)) &
    + stagger%ddzup (self%mesh, v(:,:,:,3))
  call trace_end
END FUNCTION div

!===============================================================================
!> Return the magnitude of a vector
!===============================================================================
FUNCTION magnitude_scalar (self, s) RESULT (out)
  class(scalar_t):: self
  real, dimension(:,:,:):: s
  real, dimension(size(s,1),size(s,2),size(s,3)):: out
  !.............................................................................
  call trace_begin ('scalar_t%magnitude_scalar')
  out = abs(s)
  call trace_end
END FUNCTION magnitude_scalar

!===============================================================================
!> Return the magnitude of a vector
!===============================================================================
FUNCTION magnitude_vector (self, v) RESULT (s)
  class(scalar_t):: self
  real, dimension(:,:,:,:):: v
  real, dimension(size(v,1),size(v,2),size(v,3)):: s
  type(stagger_t):: stagger
  !.............................................................................
  stagger%mesh => self%mesh
  call trace_begin ('scalar_t%magnitude_vector')
  s = stagger%xup(v(:,:,:,1))**2 &
    + stagger%yup(v(:,:,:,2))**2 &
    + stagger%zup(v(:,:,:,3))**2
  call trace_end
END FUNCTION magnitude_vector

!===============================================================================
!> Return the magnitude of a symmetric tensor, where the diagonal elements are
!> cell-centered and the off-diagonal elements are edge-centered
!===============================================================================
FUNCTION magnitude_tensor (self, t) RESULT (s)
  class(scalar_t):: self
  real, dimension(:,:,:,:,:):: t
  real, dimension(size(t,1),size(t,2),size(t,3)):: s
  type(stagger_t):: stagger
  !.............................................................................
  stagger%mesh => self%mesh
  call trace_begin ('scalar_t%magnitude_tensor')
  s = t(:,:,:,1,1)**2 &
    + t(:,:,:,2,2)**2 &
    + t(:,:,:,3,3)**2 &
    + 2.0*stagger%xup(stagger%yup(t(:,:,:,1,2)))**2 &
    + 2.0*stagger%yup(stagger%zup(t(:,:,:,2,3)))**2 &
    + 2.0*stagger%zup(stagger%xup(t(:,:,:,1,3)))**2
  call trace_end
END FUNCTION magnitude_tensor

!===============================================================================
!> Return minval of interior cells
!===============================================================================
FUNCTION minval_scalar (self, s) RESULT (out)
  class(scalar_t):: self
  real, dimension(:,:,:):: s
  real:: out
  integer:: l(3), u(3)
  !.............................................................................
  call trace_begin ('scalar_t%minval')
  l = self%mesh%li
  u = self%mesh%ui
  out = minval(s(l(1):u(1),l(2):u(2),l(3):u(3)))
  call trace_end
END FUNCTION minval_scalar
!===============================================================================

!===============================================================================
!> Return maxval of interior cells
!===============================================================================
FUNCTION maxval_scalar (self, s) RESULT (out)
  class(scalar_t):: self
  real, dimension(:,:,:):: s
  real:: out
  integer:: l(3), u(3)
  !.............................................................................
  call trace_begin ('scalar_t%maxval')
  l = self%mesh%li
  u = self%mesh%ui
  out = maxval(s(l(1):u(1),l(2):u(2),l(3):u(3)))
  call trace_end
END FUNCTION maxval_scalar
!===============================================================================

!===============================================================================
!> Return average of interior cells
!===============================================================================
FUNCTION aver (self, s) RESULT (out)
  class(scalar_t):: self
  real, dimension(:,:,:):: s
  real:: out
  integer:: l(3), u(3)
  !.............................................................................
  call trace_begin ('scalar_t%maxval')
  l = self%mesh%li
  u = self%mesh%ui
  out = sum(s(l(1):u(1),l(2):u(2),l(3):u(3)))/product(u-l+1)
  call trace_end
END FUNCTION aver
!===============================================================================

SUBROUTINE stats (self, s, label)
  class(scalar_t):: self
  real, dimension(:,:,:):: s
  character(len=*), optional:: label
  character(len=20):: l
  real:: mins, aver, maxs
  !.............................................................................
  if (io%verbose<2) return
  mins = self%minval(s)
  maxs = self%maxval(s)
  aver = self%aver(s)
  if (present(label)) then
    l = trim(label)
    print '(1x,a20,1p,5e14.5)', l, &
      minval(s), mins, aver, maxs, maxval(s)
  else
    print '(1p,5e14.5)', &
      minval(s), mins, aver, maxs, maxval(s)
  end if
END SUBROUTINE stats

END MODULE scalar_mod
