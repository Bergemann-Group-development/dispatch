!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for mesh
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE mesh_mod
  USE trace_mod
  USE io_mod
  implicit none
  !-----------------------------------------------------------------------------
  ! Mesh for patches
  !-----------------------------------------------------------------------------
  type:: mesh_t
    real(8):: b = 1d0           ! box size (periodic)
    real(8):: p = 0d0           ! position
    real(8):: s = 16d0          ! size
    real(8):: d = 1d0           ! ds
    real(8):: h(8) = 0d0        ! stagger offset for velocity and magnetic field
    real(8):: o = 0.0           ! position offset in index space
    real(8):: lf = 0.0          ! lower floating point boundary
    real(8):: uf = 0.0          ! upper floating point boundary
    real(8), dimension(:), pointer:: r => null()
    integer:: n  = 16           ! number of points
    integer:: ng =  2           ! number of guard zones
    integer:: lb =  1           ! lower index boundary
    integer:: lo =  2           ! lower index outside
    integer:: li =  3           ! lower index inside
    integer:: ui = 18           ! upper index inside
    integer:: uo = 19           ! upper index outside
    integer:: ub = 20           ! upper index boundary
    integer:: gn = 20           ! global number of points
    integer:: nh = 2            ! number of staggered variables
    integer:: id = 0
  contains
    procedure:: init       => init_mesh
    procedure:: deallocate => deallocate_mesh
    procedure:: print      => print_mesh
  end type
PRIVATE
PUBLIC mesh_t
CONTAINS

!> =============================================================================
!> The integer mesh is defined by specifying the number of points "n" and the
!> number of guard cells "ng".  The real mesh is then fully specified by giving
!> the size "s".
!> =============================================================================
SUBROUTINE init_mesh (self)
  class(mesh_t):: self
  integer:: i
  !-----------------------------------------------------------------------------
  call trace_begin ('mesh%init')
  self%lb = 1
  self%li = self%lb + self%ng
  self%gn = self%n + self%ng*2
  self%lo = max(self%li-1,1)
  self%ui = min(self%lo+self%n,self%gn)
  self%uo = min(self%ui+1,self%gn)
  self%ub = self%ui + self%ng
  if (associated(self%r)) deallocate (self%r)
  allocate (self%r(self%gn))
  self%d = self%s/self%n
  self%o = (self%li+self%ui)*0.5
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
  if (.not. io%master) return
  if (io%verbose <= 1) return
  print'(a,2(2x,1p,3e10.2))', 'mesh%print: size, position =', self%s, self%p
  if (io%verbose <= 2) return
  print'(a,3i8)',       '             n, ng, gn:', self%n, self%ng, self%gn
  print'(a,6i8)',       'lb, lo, li, ui, uo, ub:', self%lb, self%lo, self%li, self%ui, self%uo, self%ub
  print'(a,4x,4g16.8)', '        sz, ds, uf, lf:', self%s, self%d, self%uf, self%lf
  if (io%verbose <= 3) return
  do i=1,self%gn
    print*,i,self%r(i)
  end do
END SUBROUTINE

END MODULE
