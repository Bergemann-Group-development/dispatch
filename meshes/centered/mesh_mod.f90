!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for mesh
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE mesh_mod
  USE trace_mod
  implicit none
  !-----------------------------------------------------------------------------
  ! Mesh for patches
  !-----------------------------------------------------------------------------
  type:: mesh_t
    real(8):: p = 0d0
    real(8):: s = 16d0
    real(8):: d = 1d0
    real(4):: o = 0.0
    real(4):: lf = 0.0
    real(4):: uf = 0.0
    real(4), dimension(:), pointer:: r
    integer:: n  = 16
    integer:: ng =  2
    integer:: lb =  1
    integer:: lo =  2
    integer:: li =  3
    integer:: ui = 18
    integer:: uo = 19
    integer:: ub = 20
    integer:: gn = 20
    integer:: mm(3)
  contains
    procedure:: init       => init_mesh
    procedure:: deallocate => deallocate_mesh
    procedure:: print      => print_mesh
    procedure:: ddx        => ddx_mesh
    procedure:: ddy        => ddy_mesh
    procedure:: ddz        => ddz_mesh
  end type
PRIVATE
PUBLIC mesh_t
CONTAINS

!> =============================================================================
!> The integer mesh is defined by specifying the number of points "n" and the
!> number of guard cells "ng".  The real mesh is then fully specified by givin
!> the size "s".
!> =============================================================================
SUBROUTINE init_mesh (self)
  class(mesh_t):: self
  integer:: i
  !-----------------------------------------------------------------------------
  call trace_begin ('mesh%init')
  self%lb = 1
  self%li = self%lb + self%ng
  self%lo = self%li - 1
  self%ui = self%lo + self%n
  self%uo = self%ui + 1
  self%ub = self%ui + self%ng
  self%gn = self%ub
  self%mm = self%gn
  allocate (self%r(self%gn))
  self%d = self%s/self%n
  self%o = (self%n+1)*self%d*0.5
  do i=1,self%gn
    self%r(i) = (real(i,kind=8)-(self%gn+1)*0.5d0)*self%d
  end do
  self%lf = self%r(self%li) - self%d*0.5
  self%uf = self%r(self%ui) + self%d*0.5
  call self%print
  call trace_end
END SUBROUTINE

!> =============================================================================
!> =============================================================================
SUBROUTINE deallocate_mesh (self)
  class(mesh_t):: self
  !-----------------------------------------------------------------------------
  if (associated(self%r)) deallocate (self%r)
END SUBROUTINE

!> =============================================================================
!> =============================================================================
SUBROUTINE print_mesh (self)
  class(mesh_t):: self
  integer:: i
  if (verbose <= 0) return
  print*,'mesh%print:'
  print'(3i8)', self%n, self%ng, self%gn
  print'(6i8)', self%lb, self%lo, self%li, self%ui, self%uo, self%ub
  print'(4x,2g16.8)', self%uf, self%lf
  do i=1,self%gn
    print*,i,self%r(i)
  end do
END SUBROUTINE

FUNCTION ddx_mesh (self, f)
  class(mesh_t):: self
  real, dimension(self%mm(1),self%mm(2),self%mm(3)):: f, ddx_mesh
  real:: c
  c = 0.5/self%d
  ddx_mesh = (cshift(f,dim=1,shift=1)-cshift(f,dim=1,shift=-1))*c
END FUNCTION

FUNCTION ddy_mesh (self, f)
  class(mesh_t):: self
  real, dimension(self%mm(1),self%mm(2),self%mm(3)):: f, ddy_mesh
  real:: c
  c = 0.5/self%d
  ddy_mesh = (cshift(f,dim=2,shift=1)-cshift(f,dim=2,shift=-1))*c
END FUNCTION

FUNCTION ddz_mesh (self, f)
  class(mesh_t):: self
  real, dimension(self%mm(1),self%mm(2),self%mm(3)):: f, ddz_mesh
  real:: c
  c = 0.5/self%d
  ddz_mesh = (cshift(f,dim=3,shift=1)-cshift(f,dim=3,shift=-1))*c
END FUNCTION

END MODULE
