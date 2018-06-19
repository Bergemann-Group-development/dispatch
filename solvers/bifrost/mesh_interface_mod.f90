!> =============================================================================
!> Interface module, with a data type that can store all required BIFROST mesh
!> details as part of the solver_t data type, for get/put operations.
!> =============================================================================
MODULE mesh_interface_mod
  ! Dispatch modules
  USE trace_mod
  ! Bifrost modules
  USE params_mod
  USE parallel_mod
  USE stagger_mod
  implicit none
  private
  type, public:: bifrost_mesh_t
    real,dimension(:),   allocatable :: xm,xmdn,dxidxup,dxidxdn
    real,dimension(:),   allocatable :: dx1d
    real,dimension(:),   allocatable :: ym,ymdn,dyidyup,dyidydn
    real,dimension(:),   allocatable :: dy1d
    real,dimension(:),   allocatable :: zm,zmdn,dzidzup,dzidzdn
    real,dimension(:),   allocatable :: dz1d
    real,dimension(:,:), allocatable :: zlincoeff,dzlincoeff
    real,dimension(:,:), allocatable :: zuincoeff,dzuincoeff
    real :: dx,dy,dz
    integer(4):: xs,xe,xsb,xeb
    integer(4):: ys,ye,ysb,yeb
    integer(4):: zs,ze,zsb,zeb
    character(len=132) :: meshfile
  contains
    procedure:: init
    procedure:: put
    procedure:: get
  end type
  character(len=80), save:: id = &
  '$Id$'
CONTAINS

!> =============================================================================
!> Allocate patch space for Bifrost mesh data
!> =============================================================================
SUBROUTINE init (self, mb)
  class(bifrost_mesh_t):: self
  integer, optional:: mb(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('bifrost_mesh_t%init')
  allocate (self%xm      (size(xm     )))
  allocate (self%xmdn    (size(xmdn   )))
  allocate (self%dxidxup (size(dxidxup)))
  allocate (self%dxidxdn (size(dxidxdn)))
  allocate (self%dx1d    (size(dx1d   )))
  allocate (self%ym      (size(ym     )))
  allocate (self%ymdn    (size(ymdn   )))
  allocate (self%dyidyup (size(dyidyup)))
  allocate (self%dyidydn (size(dyidydn)))
  allocate (self%dy1d    (size(dy1d   )))
  allocate (self%zm      (size(zm     )))
  allocate (self%zmdn    (size(zmdn   )))
  allocate (self%dzidzup (size(dzidzup)))
  allocate (self%dzidzdn (size(dzidzdn)))
  allocate (self%dz1d    (size(dz1d   )))
  allocate (self%zlincoeff  (size(zlincoeff ,1), size(zlincoeff ,2)))
  allocate (self%dzlincoeff (size(dzlincoeff,1), size(dzlincoeff,2)))
  print *, 'bifrost_mesh_t%init:', size(xm), size(ym), size(zm)
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Put the patch values back into the Bifrost variables
!> =============================================================================
SUBROUTINE put (self)
  class(bifrost_mesh_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('bifrost_mesh_t%put')
  meshfile   = self%meshfile
  dx         = self%dx
  xm         = self%xm
  xmdn       = self%xmdn
  dxidxup    = self%dxidxup
  dxidxdn    = self%dxidxdn
  dx1d       = self%dx1d
  dy         = self%dy
  ym         = self%ym
  ymdn       = self%ymdn
  dyidyup    = self%dyidyup
  dyidydn    = self%dyidydn
  dy1d       = self%dy1d
  dz         = self%dz
  zm         = self%zm
  zmdn       = self%zmdn
  dzidzup    = self%dzidzup
  dzidzdn    = self%dzidzdn
  dz1d       = self%dz1d
  zlincoeff  = self%zlincoeff
  dzlincoeff = self%dzlincoeff
  zuincoeff  = self%zuincoeff
  dzuincoeff = self%dzuincoeff
  call trace%end
END SUBROUTINE put

!> =============================================================================
!> Copy the Bifrost mesh data to patch memory
!> =============================================================================
SUBROUTINE get (self)
  class(bifrost_mesh_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('bifrost_mesh_t%get')
  self%meshfile   = meshfile
  self%dx         = dx
  self%xm         = xm
  self%xmdn       = xmdn
  self%dxidxup    = dxidxup
  self%dxidxdn    = dxidxdn
  self%dx1d       = dx1d
  self%dy         = dy
  self%ym         = ym
  self%ymdn       = ymdn
  self%dyidyup    = dyidyup
  self%dyidydn    = dyidydn
  self%dy1d       = dy1d
  self%dz         = dz
  self%zm         = zm
  self%zmdn       = zmdn
  self%dzidzup    = dzidzup
  self%dzidzdn    = dzidzdn
  self%dz1d       = dz1d
  self%zlincoeff  = zlincoeff
  self%dzlincoeff = dzlincoeff
  self%zuincoeff  = zuincoeff
  self%dzuincoeff = dzuincoeff
  self%xs         = xs
  self%xe         = xe
  self%xsb        = xsb
  self%xeb        = xeb
  self%ys         = ys
  self%ye         = ye
  self%ysb        = ysb
  self%yeb        = yeb
  self%zs         = zs
  self%ze         = ze
  self%zsb        = zsb
  self%zeb        = zeb
  call trace%end
END SUBROUTINE get

END MODULE mesh_interface_mod
