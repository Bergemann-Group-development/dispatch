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
  USE patch_mod
  implicit none
  private
  type, extends(standard_boundaries_t), public:: boundaries_t
    ! --------------------------------------------------------------------------
    ! Experiment - specific variables
    ! --------------------------------------------------------------------------
  contains
    ! --------------------------------------------------------------------------
    ! Default procedures
    ! --------------------------------------------------------------------------
    procedure:: init
    procedure:: primary
    procedure:: derivative
    ! --------------------------------------------------------------------------
    ! Experiment - specific procedures
    ! --------------------------------------------------------------------------
    procedure:: lower_y
    procedure:: upper_y
  end type
CONTAINS

!===============================================================================
!> Initialize the boundary data type, while testing also if the patch contains
!> a spherical boundary, or possibly is completely immersed in such a boundary.
!===============================================================================
SUBROUTINE init (self)
  class(boundaries_t)     :: self
  !.............................................................................
  call self%standard_boundaries_t%init ()
END SUBROUTINE init

SUBROUTINE primary(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  if (self%mesh(2)%lower_boundary) call self%lower_y(mem)
  if (self%mesh(2)%upper_boundary) call self%upper_y(mem)

END SUBROUTINE primary

SUBROUTINE derivative(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:,:), pointer :: mem
  !
  !if (self%mesh(3)%lower_boundary) call self%lower_damp(mem)

END SUBROUTINE derivative

SUBROUTINE lower_y(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  real, dimension(:,:,:), pointer:: bx, by, bz
  integer:: l,u
  !-----------------------------------------------------------------------------
  d  => mem(:,:,:,idx%d )
  e  => mem(:,:,:,idx%e )
  px => mem(:,:,:,idx%px)
  py => mem(:,:,:,idx%py)
  pz => mem(:,:,:,idx%pz)
  bx => mem(:,:,:,idx%bx)
  by => mem(:,:,:,idx%by)
  bz => mem(:,:,:,idx%bz)
  l = self%patch%mesh(2)%lb
  u = self%patch%mesh(2)%lo

  call self%symmetric(px,idx%px, select=bits%yl)
  call self%antisymmetric(py,idx%py, select=bits%yl)
  call self%symmetric(pz,idx%pz, select=bits%yl)
  call self%symmetric(d ,idx%d , select=bits%yl)
  call self%symmetric(e ,idx%e , select=bits%yl)
  call self%symmetric(bx,idx%bx, select=bits%yl)
  call self%symmetric(by,idx%by, select=bits%yl)
  call self%symmetric(bz,idx%bz, select=bits%yl)

  if (associated(self%patch%filled)) then
    self%patch%filled(:,l:u,:) = 1
  end if

END SUBROUTINE lower_y

SUBROUTINE upper_y(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  real, dimension(:,:,:), pointer:: bx, by, bz
  integer:: l,u
  !-----------------------------------------------------------------------------
  d  => mem(:,:,:,idx%d )
  e  => mem(:,:,:,idx%e )
  px => mem(:,:,:,idx%px)
  py => mem(:,:,:,idx%py)
  pz => mem(:,:,:,idx%pz)
  bx => mem(:,:,:,idx%bx)
  by => mem(:,:,:,idx%by)
  bz => mem(:,:,:,idx%bz)
  l = self%patch%mesh(2)%uo
  u = self%patch%mesh(2)%ub

  call self%symmetric(d,idx%d, select=bits%yu)
  call self%symmetric(e,idx%e, select=bits%yu)
  call self%symmetric(px,idx%px, select=bits%yu)
  call self%antisymmetric(py,idx%py, select=bits%yu)
  call self%symmetric(pz,idx%pz, select=bits%yu)
  call self%symmetric(bx,idx%bx, select=bits%yu)
  call self%symmetric(by,idx%by, select=bits%yu)
  call self%symmetric(bz,idx%bz, select=bits%yu)

  if (associated(self%patch%filled)) then
    self%patch%filled(:,l:u,:) = 1
  end if

END SUBROUTINE upper_y


END MODULE boundaries_mod
