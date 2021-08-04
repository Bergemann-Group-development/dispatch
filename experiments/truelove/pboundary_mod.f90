!===============================================================================
!> Boundary conditions for a 3-D box self-gravitating box
!> Assume `phi` = 0 in the boundaries, all other variables outflow.
!===============================================================================
MODULE pboundary_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE gpatch_mod
  USE bits_mod
  implicit none
  private

  type, public:: pboundary_t
  contains
    procedure:: init
    procedure:: conditions
  end type
  type(pboundary_t), public:: pboundaries

  logical, save:: fixedbcs = .false.

CONTAINS

!===============================================================================
SUBROUTINE init (self, patch)
  class(pboundary_t):: self
  class(gpatch_t):: patch
  logical, save:: first_time=.true.
  integer:: iostat
  namelist /boundary_params/ fixedbcs
  !-----------------------------------------------------------------------------
  !$omp critical (boundary_cr)
  if (first_time) then
    rewind (io%input)
    read (io%input, boundary_params, iostat=iostat)
    call io%namelist_warning ('boundary_params', iostat)
    write (io%namelist, boundary_params)
    first_time = .false.
  end if
  !$omp end critical (boundary_cr)
  if (fixedbcs) then
    patch%periodic = .false.
    if (patch%is_set(bits%root_grid)) then
      call patch%boundaries%set(bits%xl)
      call patch%boundaries%set(bits%xu)
      call patch%boundaries%set(bits%yl)
      call patch%boundaries%set(bits%yu)
      call patch%boundaries%set(bits%zl)
      call patch%boundaries%set(bits%zu)
    end if
  end if
END SUBROUTINE init

!===============================================================================
!> Impose boundary conditions.
!> All conditions are assumed "outflow" (i.e., zero gradient).
!===============================================================================
SUBROUTINE conditions (self, patch, it)
  class(pboundary_t) :: self
  class(gpatch_t), target :: patch
  class(mesh_t), pointer :: mx, my, mz
  integer :: ix, iy, iz, it
  !-----------------------------------------------------------------------------
  if (.not. fixedbcs) return
  call trace%begin('pboundary_t%conditions')

  mx => patch%mesh(1)
  my => patch%mesh(2)
  mz => patch%mesh(3)

  ! x-boundaries
  if (patch%n(1) > 1) then
    ! lower x
    if (patch%mesh(1)%lower_boundary) then
      do ix=mx%lb,mx%lo
        patch%mem(ix,my%lb:my%ub,mz%lb:mz%ub,:,it,1) = patch%mem(mx%li+1,my%lb:my%ub,mz%lb:mz%ub,:,it,1)
      end do
      if (patch%idx%phi > 0) then
        do ix=mx%lb,mx%lo
          patch%mem(ix,my%lb:my%ub,mz%lb:mz%ub,patch%idx%phi,it,1) = 0.0
        end do
      end if
      if (associated(patch%filled)) &
        patch%filled(1:mx%lo,my%lb:my%ub,mz%lb:mz%ub) = patch%level
    end if
  ! upper x
    if (patch%mesh(1)%upper_boundary) then
      do ix=mx%uo,mx%ub
        patch%mem(ix,my%lb:my%ub,mz%lb:mz%ub,:,it,1) = patch%mem(mx%ui-1,my%lb:my%ub,mz%lb:mz%ub,:,it,1)
      end do
      if (patch%idx%phi > 0) then
        do ix=mx%uo,mx%ub
          patch%mem(ix,my%lb:my%ub,mz%lb:mz%ub,patch%idx%phi,it,1) = 0.0
        end do
      end if
      if (associated(patch%filled)) &
        patch%filled(mx%uo:mx%ub,my%lb:my%ub,mz%lb:mz%ub) = patch%level
    end if
  end if

  ! y-boundaries
  if (patch%n(2) > 1) then
    ! lower y
    if (patch%mesh(2)%lower_boundary) then
      do iy=my%lb,my%lo
        patch%mem(mx%lb:mx%ub,iy,mz%lb:mz%ub,:,it,1) = patch%mem(mx%lb:mx%ub,my%li+1,mz%lb:mz%ub,:,it,1)
      end do
      if (patch%idx%phi > 0) then
        do iy=my%lb,my%lo
          patch%mem(mx%lb:mx%ub,iy,mz%lb:mz%ub,patch%idx%phi,it,1) = 0.0
        end do
      end if
      if (associated(patch%filled)) &
        patch%filled(mx%lb:mx%ub,1:my%lo,mz%lb:mz%ub) = patch%level
    end if
  ! upper y
    if (patch%mesh(2)%upper_boundary) then
      do iy=my%uo,my%ub
        patch%mem(mx%lb:mx%ub,iy,mz%lb:mz%ub,:,it,1) = patch%mem(mx%lb:mx%ub,my%ui-1,mz%lb:mz%ub,:,it,1)
      end do
      if (patch%idx%phi > 0) then
        do iy=my%uo,my%ub
          patch%mem(mx%lb:mx%ub,iy,mz%lb:mz%ub,patch%idx%phi,it,1) = 0.0
        end do
      end if
      if (associated(patch%filled)) &
        patch%filled(mx%lb:mx%ub,my%uo:my%ub,mz%lb:mz%ub) = patch%level
    end if
  end if

  ! z-boundaries
  if (patch%n(3) > 1) then
    ! lower z
    if (patch%mesh(3)%lower_boundary) then
      do iz=mz%lb,mz%lo
        patch%mem(mx%lb:mx%ub,my%lb:my%ub,iz,:,it,1) = patch%mem(mx%lb:mx%ub,my%lb:my%ub,mz%li+1,:,it,1)
      end do
      if (patch%idx%phi > 0) then
        do iz=mz%lb,mz%lo
          patch%mem(mx%lb:mx%ub,my%lb:my%ub,iz,patch%idx%phi,it,1) = 0.0
        end do
      end if
      if (associated(patch%filled)) &
        patch%filled(mx%lb:mx%ub,my%lb:my%ub,1:mz%lo) = patch%level
    end if
    ! upper z
    if (patch%mesh(3)%upper_boundary) then
      do iz=mz%uo,mz%ub
        patch%mem(mx%lb:mx%ub,my%lb:my%ub,iz,:,it,1) = patch%mem(mx%lb:mx%ub,my%lb:my%ub,mz%ui-1,:,it,1)
      end do
      if (patch%idx%phi > 0) then
        do iz=mz%uo,mz%ub
          patch%mem(mx%lb:mx%ub,my%lb:my%ub,iz,patch%idx%phi,it,1) = 0.0
        end do
      end if
      if (associated(patch%filled)) &
        patch%filled(mx%lb:mx%ub,my%lb:my%ub,mz%uo:mz%ub) = patch%level
    end if
  end if
  call trace%end ()
END SUBROUTINE conditions

END MODULE pboundary_mod
