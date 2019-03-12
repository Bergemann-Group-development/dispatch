!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for mesh
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE mesh_mod
  USE trace_mod
  USE io_mod
  implicit none
  !-----------------------------------------------------------------------------
  ! Mesh for patches
  !-----------------------------------------------------------------------------
  type:: mesh_t
    real(8)::  p(3) = 0d0
    real(8)::  s(3) = 1d0
    real(8)::  d(3) = 1d0
    real(4)::  o(3) = 0.0
    real(4):: lf(3) = 0.0
    real(4):: uf(3) = 0.0
    integer:: n(3)  = 16
    integer:: ng(3) =  2
    integer:: lb(3) =  1
    integer:: lo(3) =  2
    integer:: li(3) =  3
    integer:: ui(3) = 18
    integer:: uo(3) = 19
    integer:: ub(3) = 20
    integer:: gn(3) = 20
    real(4), dimension(:,:), pointer:: r
    contains
    procedure:: init       => init_mesh
    procedure:: deallocate => deallocate_mesh
    procedure:: print      => print_mesh
  end type
PRIVATE
PUBLIC mesh_t
CONTAINS

!===============================================================================
!> The integer mesh is defined by specifying the number of points "n" and the
!> number of guard cells "ng".  The real mesh is then fully specified by givin
!> the size "s".
!===============================================================================
SUBROUTINE init_mesh (self, position, size, n, ng)
  real(8), dimension(3), optional:: position, size
  integer, dimension(3), optional:: n, ng
  class(mesh_t):: self
  integer:: i, dim
  !-----------------------------------------------------------------------------
  call trace_begin ('mesh%init')
  self%n  =  16; if (present(n )) self%n  = n
  self%ng =   2; if (present(ng)) self%ng = ng
  self%s  = 1d0; if (present(size)) self%s  = size
  self%p  = 0d0; if (present(position)) self%p  = position
  print*,self%n
  print*,self%ng
  print*,self%s
  print*,self%p
  !-----------------------------------------------------------------------------
  self%li = self%lb + self%ng
  self%lo = self%li - 1
  self%ui = self%lo + self%n
  self%uo = self%ui + 1
  self%ub = self%ui + self%ng
  self%gn = self%ub
  allocate (self%r(maxval(self%gn),3))
  self%d = self%s/self%n
  self%o = (self%n+1)*self%d*0.5
  do dim=1,3
    do i=1,self%gn(dim)
      self%r(i,dim) = (real(i,kind=8)-(self%gn(dim)+1)*0.5d0)*self%d(dim)
    end do
    self%lf(dim) = self%r(self%li(dim),dim) - self%d(dim)*0.5
    self%uf(dim) = self%r(self%ui(dim),dim) + self%d(dim)*0.5
  end do
  call self%print
  call trace_end
END SUBROUTINE

!===============================================================================
!===============================================================================
SUBROUTINE deallocate_mesh (self)
  class(mesh_t):: self
  !-----------------------------------------------------------------------------
  if (associated(self%r)) deallocate (self%r)
END SUBROUTINE

!===============================================================================
!===============================================================================
SUBROUTINE print_mesh (self)
  class(mesh_t):: self
  integer:: i
  if (io%verbose <= 0) return
  print*,'mesh%print:'
  print'(a,3(3i8,2x))', 'n, ng, gn:', self%n, self%ng, self%gn
  print'(a,6(3i8,2x))', 'lb, lo, li, iu, uo, ub:', self%lb, self%lo, self%li, self%ui, self%uo, self%ub
  print'(a,2(3g16.8))', 'lf, uf:', self%lf, self%uf
  do i=1,self%gn(1)
    print'(i4,3f10.6)',i,self%r(i,:)
  end do
END SUBROUTINE

END MODULE
