!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Define a patch of a larger cube, which is split into n = [n1,n2,n2] patches.
!> Patches have meshes, and methods for computing derivatives
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE patch_mod
  USE io_mod
  USE mpi_mod
  USE trace_mod
  USE mesh_mod
  USE derivatives_mod
  USE upwind_mod
  USE central_mod
  implicit none
  private
  type, public:: patch_t
    integer:: id = 0
    real,    dimension(3):: size = 1.0
    integer, dimension(3):: n = [1,1,1]
    type(mesh_t), pointer:: mesh
    class(derivatives_t), pointer:: d
    logical:: do_upwind
  contains
    procedure:: init
    !procedure:: update
  end type
  integer:: n(3) = [4,4,1]
  real:: size(3) = 1.0, factor = 1.5
  logical:: do_upwind = .true., first_time = .true.
  integer:: id = 0
CONTAINS
!> =============================================================================
SUBROUTINE init (self, size, n, i, j, k)
  implicit none
  class(patch_t):: self
  real:: size(3)
  integer:: i, j, k, n(3)
  type(upwind_t), pointer:: upwind
  type(central_t), pointer:: central
  namelist /patch_params/ do_upwind, factor
  !-----------------------------------------------------------------------------
  ! Read patch parameters and copy to self
  !-----------------------------------------------------------------------------
  call trace_begin ('patch_t%init')
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, patch_params)
    if (io%master) write (*, patch_params)
  end if
  call mpi%barrier(delay=5)
  self%n = n
  self%size = size/n
  !-----------------------------------------------------------------------------
  ! Initialize mesh and place us in the unit cube
  !-----------------------------------------------------------------------------
  allocate (self%mesh)
  call self%mesh%init (self%size, self%n)
  call self%mesh%place_in_box (ii=[i,j,k])
  !-----------------------------------------------------------------------------
  ! Initialize derivative operators
  !-----------------------------------------------------------------------------
  if (do_upwind) then
    allocate (upwind)
    self%d => upwind
    self%d%factor = factor
  else
    allocate (central)
    self%d => central
  end if
  call self%d%init (self%mesh)
  call trace_end
END SUBROUTINE init
!> =============================================================================
END MODULE patch_mod
