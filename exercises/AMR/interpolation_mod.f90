!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE interpolation_mod
  USE io_mod
  USE trace_mod
  USE patches_mod
  implicit none
  private
  type, public :: interpolation_t
    real   :: factor=1.5
    logical:: do_upwind=.false.
  contains
    procedure:: init
    procedure:: prolongate
    procedure:: restrict
  end type
  logical:: first_time = .true.
CONTAINS
!> =============================================================================
SUBROUTINE init (self)
  class(interpolation_t):: self
  integer:: dims(3)
  real:: ds(3)
  logical:: do_upwind = .false.
  real:: factor = 1.5
  namelist /interpol_params/ do_upwind, factor
  !.............................................................................
  call trace_begin('interpolation_t%init', 1)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, interpol_params)
    if (io%master) write (*, interpol_params)
  end if
  self%do_upwind = do_upwind
  self%factor = factor
  call trace_end
END SUBROUTINE init
!> =============================================================================
!> Prolongate data from a coarse patch to a fine patch.
!> init is a switch to indicate we need to prolongate all the patch (initial refinement)
!> guard is a switch that signals updating of the guard cells
!> =============================================================================
SUBROUTINE prolongate (self, coarse, fine, init, guard)
  class(interpolation_t):: self
  class(patch_t):: coarse, fine
  logical :: init, guard
  !.............................................................................
  call trace_begin ('interpolation_t%prolongate', 3)
  if (first_time) then
    first_time = .false.
  end if
  call trace_end
END SUBROUTINE prolongate
!> =============================================================================
!> Restrict data from a fine patch to a coarse patch.
!> upload is a switch to indicate we need to restrict in interior of the fine patch (upload of data)
!> guard is a switch that signals updating of the guard cells
!> =============================================================================
SUBROUTINE restrict (self, coarse, fine, upload, guard)
  class(interpolation_t):: self
  class(patch_t):: coarse, fine
  logical :: upload, guard
  integer :: ix,iy,iz,jx,jy,jz,kx,ky,kz,lb(3),ub(3)
  real    :: inp
  !.............................................................................
  call trace_begin ('interpolation_t%prolongation', 3)
  if (first_time) then
    first_time = .false.
  end if
  if (upload) then
    lb = merge(coarse%li,coarse%li + coarse%n/2,   fine%position < coarse%position) ! lower index in coarse cell
    ub = merge(coarse%li + coarse%n/2-1,coarse%ui, fine%position < coarse%position) ! upper index in coarse cell
    where(coarse%n==1) ! redundant, but better guard against round-off errors somewhere
      lb=1
      ub=1
    endwhere
    inp = 1./product(merge(2.,1.,coarse%n > 1))                                 ! division factor
    ! This can be done more efficient using the seven different cases (1D, 2D, 3D in arbitrary configs)
    do iz=lb(3),ub(3)
      jz = iz-lb(3) + fine%li(3)
      kz = jz + merge(1,0,fine%n(3)>1)
    do iy=lb(2),ub(2)
      jy = iy-lb(2) + fine%li(2)
    do ix=lb(1),ub(1)
      jx = ix-lb(1) + fine%li(1)
      coarse%mem(ix,iy,iz,:,coarse%it,:) = sum(fine%mem(jx:kx,jy:ky,jz:kz,:,fine%it,:)) * inp
    enddo
    enddo
    enddo
  endif
  if (guard) then
    ! test for overlaps and check where data has been uploaded already
  endif
  call trace_end
END SUBROUTINE restrict
!
END MODULE interpolation_mod
