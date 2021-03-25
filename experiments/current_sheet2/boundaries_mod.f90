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
  USE timer_mod
  USE index_mod
  USE math_mod
  USE standard_boundaries_mod
  USE eos_mod
  USE patch_mod
  USE stagger_mod
  implicit none
  private
  type, extends(standard_boundaries_t), public:: boundaries_t
    ! --------------------------------------------------------------------------
    ! Experiment - specific variables
    ! --------------------------------------------------------------------------
    real(kind=KindScalarVar)       :: t0, v0
    real(kind=KindScalarVar), dimension(:,:,:,:),allocatable::ic_copy 
    integer, dimension(3) :: gn
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
    procedure:: lower_x
    procedure:: upper_x
    procedure:: lower_y
    procedure:: upper_y
    procedure:: lower_z
    procedure:: upper_z
    procedure:: copy_ic
    procedure:: static_values
    procedure:: add_flow
    procedure:: vanishing
  end type
  integer, save:: verbose=1
  ! --- add_flow parameters
  real, save:: t0 = 100.
  real, save:: v0 = 0.045
CONTAINS

!===============================================================================
!> Initialize the boundary data type, while testing also if the patch contains
!> a spherical boundary, or possibly is completely immersed in such a boundary.
!===============================================================================
SUBROUTINE init (self, patch)
  class(boundaries_t)   :: self
  class(patch_t), target:: patch
  !.............................................................................
  integer:: iostat
  real, save:: t0, v0
  logical, save:: first_time=.true.
  namelist /boundary_params/ verbose, t0, v0
  !-----------------------------------------------------------------------------
  call trace%begin('boundaries_t%init')
  !.............................................................................
  call io%lock%set('boundaries_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, boundary_params, iostat=iostat)
    if (io%master) write (*, boundary_params)
    first_time = .false.
  end if
  self%t0 = t0
  self%v0 = v0

  self%gn = patch%gn
  call io%lock%unset('boundaries_t%init')
  call self%standard_boundaries_t%init(patch)
END SUBROUTINE init

SUBROUTINE primary(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  if (self%mesh(2)%lower_boundary) call self%lower_y(mem)
  if (self%mesh(2)%upper_boundary) call self%upper_y(mem)
  if (self%mesh(1)%lower_boundary) call self%lower_x(mem)
  if (self%mesh(1)%upper_boundary) call self%upper_x(mem)
  if (self%mesh(3)%lower_boundary) call self%lower_z(mem)
  if (self%mesh(3)%upper_boundary) call self%upper_z(mem)
END SUBROUTINE primary

SUBROUTINE derivative(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:,:), pointer :: mem
END SUBROUTINE derivative

SUBROUTINE lower_x(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  real, dimension(:,:,:), pointer:: bx, by, bz
  !-----------------------------------------------------------------------------
  d  => mem(:,:,:,idx%d )
  e  => mem(:,:,:,idx%e )
  px => mem(:,:,:,idx%px)
  py => mem(:,:,:,idx%py)
  pz => mem(:,:,:,idx%pz)
  bx => mem(:,:,:,idx%bx)
  by => mem(:,:,:,idx%by)
  bz => mem(:,:,:,idx%bz)

  call self%antisymmetric(px,idx%px, select=bits%xl)
  call self%antisymmetric(py,idx%py, select=bits%xl)
  call self%antisymmetric(pz,idx%pz, select=bits%xl)
  call self%static_values(d ,idx%d , select=bits%xl)
  call self%static_values(e ,idx%e , select=bits%xl)
  call self%static_values(bx,idx%bx, select=bits%xl)
  call self%static_values(by,idx%by, select=bits%xl)
  call self%static_values(bz,idx%bz, select=bits%xl)
END SUBROUTINE lower_x  

SUBROUTINE upper_x(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  real, dimension(:,:,:), pointer:: bx, by, bz
  !-----------------------------------------------------------------------------
  d  => mem(:,:,:,idx%d )
  e  => mem(:,:,:,idx%e )
  px => mem(:,:,:,idx%px)
  py => mem(:,:,:,idx%py)
  pz => mem(:,:,:,idx%pz)
  bx => mem(:,:,:,idx%bx)
  by => mem(:,:,:,idx%by)
  bz => mem(:,:,:,idx%bz)

  call self%symmetric(d,idx%d, select=bits%xu)
  call self%symmetric(e,idx%e, select=bits%xu)
  call self%antisymmetric(px,idx%px, select=bits%xu)
  call self%antisymmetric(py,idx%py, select=bits%xu)
  call self%antisymmetric(pz,idx%pz, select=bits%xu)
  call self%extrapolate(bx,idx%bx, select=bits%xu)
  call self%extrapolate(by,idx%by, select=bits%xu)
  call self%extrapolate(bz,idx%bz, select=bits%xu)
  
END SUBROUTINE upper_x

SUBROUTINE lower_y(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  real, dimension(:,:,:), pointer:: bx, by, bz
  real, dimension(:,:,:), allocatable:: vx, vy, vz
  real, dimension(:,:,:), allocatable:: xdnd, ydnd, zdnd
  !-----------------------------------------------------------------------------

  d  => mem(:,:,:,idx%d )
  e  => mem(:,:,:,idx%e )
  px => mem(:,:,:,idx%px)
  py => mem(:,:,:,idx%py)
  pz => mem(:,:,:,idx%pz)
  bx => mem(:,:,:,idx%bx)
  by => mem(:,:,:,idx%by)
  bz => mem(:,:,:,idx%bz)

  allocate(xdnd,ydnd,zdnd, mold=d)
  xdnd=0; ydnd=0; zdnd=0
  ! Convert to velocity
  xdnd = xdn(d); ydnd = ydn(d); zdnd = zdn(d)
  px = px/xdnd; py = py/ydnd; pz = pz/zdnd
  call self%vanishing(px,idx%px, select=bits%yl)
  call self%vanishing(py,idx%py, select=bits%yl)
  call self%vanishing(pz,idx%pz, select=bits%yl)
  ! Back to momentum
  px = px*xdnd; py = py*ydnd; pz = pz*zdnd
  call self%static_values(d ,idx%d , select=bits%yl)
  call self%static_values(e ,idx%e , select=bits%yl)
  call self%static_values(bx,idx%bx, select=bits%yl)
  call self%static_values(by,idx%by, select=bits%yl)
  call self%static_values(bz,idx%bz, select=bits%yl)
  deallocate(xdnd,ydnd,zdnd)

END SUBROUTINE lower_y

SUBROUTINE upper_y(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  real, dimension(:,:,:), pointer:: bx, by, bz
  real, dimension(:,:,:), allocatable:: xdnd
  !-----------------------------------------------------------------------------
  d  => mem(:,:,:,idx%d )
  e  => mem(:,:,:,idx%e )
  px => mem(:,:,:,idx%px)
  py => mem(:,:,:,idx%py)
  pz => mem(:,:,:,idx%pz)
  bx => mem(:,:,:,idx%bx)
  by => mem(:,:,:,idx%by)
  bz => mem(:,:,:,idx%bz)

  call self%symmetric(d,idx%d, select=bits%yu)
  call self%symmetric(e,idx%e, select=bits%yu)
  call self%antisymmetric(px,idx%px, select=bits%yu)
  call self%antisymmetric(py,idx%py, select=bits%yu)
  call self%antisymmetric(pz,idx%pz, select=bits%yu)
  call self%symmetric(bx,idx%bx, select=bits%yu)
  call self%symmetric(by,idx%by, select=bits%yu)
  call self%symmetric(bz,idx%bz, select=bits%yu)

  allocate(xdnd, mold=d)
  xdnd=0
  ! Convert to velocity
  xdnd = xdn(d)
  px = px/xdnd
  call self%add_flow(px, select=bits%yu)
  ! Back to momentum
  px = px*xdnd
  deallocate(xdnd)
END SUBROUTINE upper_y

SUBROUTINE lower_z(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  real, dimension(:,:,:), pointer:: bx, by, bz
  !-----------------------------------------------------------------------------
  d  => mem(:,:,:,idx%d )
  e  => mem(:,:,:,idx%e )
  px => mem(:,:,:,idx%px)
  py => mem(:,:,:,idx%py)
  pz => mem(:,:,:,idx%pz)
  bx => mem(:,:,:,idx%bx)
  by => mem(:,:,:,idx%by)
  bz => mem(:,:,:,idx%bz)
  
  call self%antisymmetric(px,idx%px, select=bits%zl)
  call self%antisymmetric(py,idx%py, select=bits%zl)
  call self%antisymmetric(pz,idx%pz, select=bits%zl)
  call self%static_values(d ,idx%d , select=bits%zl)
  call self%static_values(e ,idx%e , select=bits%zl)
  call self%static_values(bx,idx%bx, select=bits%zl)
  call self%static_values(by,idx%by, select=bits%zl)
  call self%static_values(bz,idx%bz, select=bits%zl)

END SUBROUTINE lower_z

SUBROUTINE upper_z(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  real, dimension(:,:,:), pointer:: bx, by, bz
  !-----------------------------------------------------------------------------
  d  => mem(:,:,:,idx%d )
  e  => mem(:,:,:,idx%e )
  px => mem(:,:,:,idx%px)
  py => mem(:,:,:,idx%py)
  pz => mem(:,:,:,idx%pz)
  bx => mem(:,:,:,idx%bx)
  by => mem(:,:,:,idx%by)
  bz => mem(:,:,:,idx%bz)

  call self%symmetric(d,idx%d, select=bits%zu)
  call self%symmetric(e,idx%e, select=bits%zu)
  call self%antisymmetric(px,idx%px, select=bits%zu)
  call self%antisymmetric(py,idx%py, select=bits%zu)
  call self%antisymmetric(pz,idx%pz, select=bits%zu)
  call self%extrapolate(bx,idx%bx, select=bits%zu)
  call self%extrapolate(by,idx%by, select=bits%zu)
  call self%extrapolate(bz,idx%bz, select=bits%zu)

END SUBROUTINE upper_z
!===============================================================================
!> Copies the initial conditions for use in a static boundary
!===============================================================================
SUBROUTINE copy_ic(self,mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar) :: mem(:,:,:,:)
  ! ----------------------------------------------------------------------------
  allocate(self%ic_copy, mold=mem)
  self%ic_copy = mem
END SUBROUTINE copy_ic
!===============================================================================
!> Static boundary condition
!===============================================================================
SUBROUTINE static_values(self, mem, iv, select)
  class(boundaries_t)      :: self
  class(mesh_t), pointer   :: m(:)
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: mem
  integer :: lb, lo, uo, ub, iv
  integer, optional :: select
  !-----------------------------------------------------------------------------
  m => self%mesh

  if (self%is_set(bits%xl)) then
    lb = m(1)%lb
    lo = m(1)%lo
    mem(lb:lo,:,:) = self%ic_copy(lb:lo,:,:,iv)
  else if (self%is_set(bits%xu)) then
    uo = m(1)%uo
    ub = m(1)%ub
    mem(uo:ub,:,:) = self%ic_copy(uo:ub,:,:,iv)
  else if (self%is_set(bits%yl)) then
    lb = m(2)%lb
    lo = m(2)%lo
    mem(:,lb:lo,:) = self%ic_copy(:,lb:lo,:,iv)
  else if (self%is_set(bits%yu)) then
    uo = m(2)%uo
    ub = m(2)%ub
    mem(:,uo:ub,:) = self%ic_copy(:,uo:ub,:,iv)
  else if (self%is_set(bits%zl)) then
    lb = m(3)%lb
    lo = m(3)%lo
    mem(:,:,lb:lo) = self%ic_copy(:,:,lb:lo,iv)
  else if (self%is_set(bits%zu)) then
    uo = m(3)%uo
    ub = m(3)%ub
    mem(:,:,uo:ub) = self%ic_copy(:,:,uo:ub,iv)
  end if
END SUBROUTINE static_values
!===============================================================================
!> Vanishing boundary condition, more gentle than antisymmetric
!===============================================================================
SUBROUTINE vanishing(self, mem, iv, select)
  class(boundaries_t)      :: self
  class(mesh_t), pointer   :: m(:)
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: mem
  integer :: lb, lo, uo, ub, iv, x, y, z
  integer, optional :: select
  !-----------------------------------------------------------------------------
  m => self%mesh
  if (self%is_set(bits%xl)) then
    lb = m(1)%lb
    lo = m(1)%lo
    do x=lo,lb,-1
      mem(x,:,:) = mem(x+1,:,:)*0.7
    end do
  else if (self%is_set(bits%xu)) then
    ub = m(1)%ub
    uo = m(1)%uo
    do x=uo,ub
      mem(x,:,:) = mem(x-1,:,:)*0.7
    end do
  else if (self%is_set(bits%yl)) then
    lb = m(2)%lb
    lo = m(2)%lo
    do y=lo,lb,-1
      mem(:,y,:) = mem(:,y+1,:)*0.7
    end do
  else if (self%is_set(bits%yu)) then
    ub = m(2)%ub
    uo = m(2)%uo
    do y=uo,ub
      mem(:,y,:) = mem(:,y-1,:)*0.7
    end do
  else if (self%is_set(bits%zl)) then
    lb = m(3)%lb
    lo = m(3)%lo
    do z=lo,lb,-1
      mem(:,:,z) = mem(:,:,z+1)*0.7
    end do
  else if (self%is_set(bits%zu)) then
    ub = m(3)%ub
    uo = m(3)%uo
    do z=uo,ub
      mem(:,:,z) = mem(:,:,z-1)*0.7
    end do
  end if
END SUBROUTINE vanishing
!===============================================================================
!> Add a flow to the boundary
!===============================================================================
SUBROUTINE add_flow(self, mem, select)
  class(boundaries_t)     :: self
  class(mesh_t), pointer  :: m(:)
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: mem
  integer                 :: ub, uo
  real                    :: t0, v0, flow, t
  integer, optional       :: select
  !-----------------------------------------------------------------------------
  m => self%mesh
  t = timer%time
  t0 = self%t0
  v0 = self%v0

  if (t<t0) then
    flow = v0/2.*(1.-cos(math%pi*t/t0))
  else
    flow = v0
  end if

  if (self%is_set(bits%xu)) then
    uo = m(1)%uo
    ub = m(1)%ub
    mem(uo:ub,:,:) = flow
  else if (self%is_set(bits%yu)) then
    uo = m(2)%uo
    ub = m(2)%ub
    mem(:,uo:ub,:) = flow
  else if (self%is_set(bits%zu)) then
    uo = m(3)%uo
    ub = m(3)%ub
    mem(:,:,uo:ub) = flow
  end if
END SUBROUTINE add_flow

END MODULE boundaries_mod
