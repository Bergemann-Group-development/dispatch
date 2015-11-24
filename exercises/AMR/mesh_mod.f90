!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> A mesh object describes the (approximate in the case of curvilinear patch)
!> size and position of a patch.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE mesh_mod
  USE trace_mod
  implicit none
  private
  type, public:: mesh_t
    real(8):: s(3)  = 1.0_8     ! size (interior)
    real(8):: gs(3)             ! global size
    real(8):: d(3)              ! cell size
    real(8):: p(3)  = 0.5_8     ! mid position
    real(8):: b(3)  = 1.0_8     ! box size
    real(8):: o(3)  = 0.0_8     ! origin
    real(8), dimension(:,:), pointer:: r => null()
    integer:: n(3)  = 16        ! number of interior cells
    integer:: ng(3) = 2         ! number of guards cells
    integer:: gn(3) = 20        ! number of global cells
  contains
    procedure:: init            ! initialize d, gs, gn, and r
    procedure:: place_in_box    ! place the mesh at index ii in box with n
    procedure:: distance        ! return the distance btw two patches
    procedure:: overlaps        ! return true if they overlap (incl. ghost zones)
    procedure:: index_range     ! returns the index ranges btw two points
    procedure:: indices         ! computes indices and fractional indices
    procedure:: index           ! compute the index for a single point
  end type
CONTAINS
!> =============================================================================
SUBROUTINE init (self, size, n)
  class(mesh_t):: self
  real:: size(3)
  integer:: i, j, n(3)
  !.............................................................................
  call trace_begin ('mesh_t%init')
  self%s  = size
  self%n  = n
  self%gn = self%n + 2*self%ng
  self%d  = self%s/self%n
  self%gs = self%s + 2*self%ng*self%d
  if (associated(self%r)) deallocate(self%r)
  do i=1,3
    allocate (self%r(self%gn(i),3))
    do j=1,self%gn(i)
      self%r(j,i) = -(self%gn(i)+1)*0.5_8 + real(j,kind=8)
    end do
  end do
  call trace_end
END SUBROUTINE init
!> =============================================================================
SUBROUTINE place_in_box (self, ii)
  class(mesh_t):: self
  integer:: ii(3)
  !.............................................................................
  self%p = -0.5_8*self%b + (ii-0.5_8)*self%s
END SUBROUTINE place_in_box
!> =============================================================================
FUNCTION index_range (self, p)
  class(mesh_t):: self
  real(8):: p(2)
  !.............................................................................
  integer:: index_range(3,2)
  index_range(:,1) = max(1, min(self%gn, int(modulo(p(1)-self%o, self%b/self%d))))
  index_range(:,2) = max(1, min(self%gn, int(modulo(p(2)-self%o, self%b/self%d))))
END FUNCTION index_range
!> =============================================================================
SUBROUTINE indices (self, np, p, ii, f)
  class(mesh_t):: self
  real(8):: p(np), f(np,3)
  integer:: np, i, ii(np,3)
  !.............................................................................
  do i=1,3
    f(:,i) = modulo(p-self%o(i), self%b(i))/self%d(i)
    ii(:,i) = max(1, min(self%gn(i), int(f(:,i))))
  end do
  f = f-ii
END SUBROUTINE indices
!> =============================================================================
FUNCTION index (self, p)
  class(mesh_t):: self
  real(8), dimension(3):: p, f
  integer:: index(3)
  !.............................................................................
  f = modulo(p-self%o, self%b)/self%d
  index = max(1, min(self%gn, int(f)))
  f = f - index
END FUNCTION index
!> =============================================================================
FUNCTION distance (self, mesh)
  class(mesh_t):: self, mesh
  real(8):: distance(3)
  distance = abs(modulo(self%p-mesh%p+0.5*self%b,self%b)-0.5*self%b)
END FUNCTION distance
!> =============================================================================
FUNCTION overlaps (self, mesh)
  class(mesh_t):: self, mesh
  logical:: overlaps
  !.............................................................................
  overlaps = all(distance(self,mesh) < 0.5000000001_8*(self%gs+mesh%gs))
END FUNCTION
END MODULE mesh_mod
