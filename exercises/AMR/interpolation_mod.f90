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
  integer :: offset(3),ix,iy,iz,cx,cy,cz,ox,oy,oz
  logical :: flat(3)
  !.............................................................................
  call trace_begin ('interpolation_t%prolongate', 3)
  if (first_time) then
    first_time = .false.
  end if

  if (guard) then
    ! test for overlaps and check where data has been uploaded already
  endif

  if (init) then
    offset = merge(coarse%li-1,coarse%li + coarse%n/2-1,fine%position < coarse%position)
    flat = coarse%n == 1
    if (flat(3)) then   ! z-dim is flat
      if (flat(2)) then ! y-dim is flat: 1D run in x-direction
        do ix=fine%li(1), fine%ui(1)
          cx = (ix-1)/2 + offset(1)
          ox = cx + merge(-1,1,modulo(ix-fine%li(1),2)==0)
          fine%mem(ix,1,1,:,fine%it,:) = 0.75 * coarse%mem(cx,1,1,:,coarse%it,:) + &
                                         0.25 * coarse%mem(ox,1,1,:,coarse%it,:)
        enddo
      else
        if (flat(1)) then ! x-dim is flat: 1D run in y-direction
          do iy=fine%li(2), fine%ui(2)
            cy = (iy-1)/2 + offset(2)
            oy = cy + merge(-1,1,modulo(iy-fine%li(2),2)==0)
            fine%mem(1,iy,1,:,fine%it,:) = 0.75 * coarse%mem(1,cy,1,:,coarse%it,:) + &
                                           0.25 * coarse%mem(1,oy,1,:,coarse%it,:)
          enddo
        else ! 2D run in x-y plane
          do iy=fine%li(2), fine%ui(2)
            cy = (iy-1)/2 + offset(2)
            oy = cy + merge(-1,1,modulo(iy-fine%li(2),2)==0)
          do ix=fine%li(1), fine%ui(1)
            cx = (ix-1)/2 + offset(1)
            ox = cx + merge(-1,1,modulo(ix-fine%li(1),2)==0)
            fine%mem(ix,iy,1,:,fine%it,:) = 0.75*0.75 * coarse%mem(cx,cy,1,:,coarse%it,:) + &
                                            0.25*0.75 * coarse%mem(ox,cy,1,:,coarse%it,:) + &
                                            0.75*0.25 * coarse%mem(cx,oy,1,:,coarse%it,:) + &
                                            0.25*0.25 * coarse%mem(ox,oy,1,:,coarse%it,:)
          enddo
          enddo
        endif
      endif
    else
      if (flat(2)) then
        if (flat(1)) then ! 1D run in z-dir
          do iz=fine%li(3), fine%ui(3)
            cz = (iz-1)/2 + offset(3)
            oz = cz + merge(-1,1,modulo(iz-fine%li(3),2)==0)
            fine%mem(1,1,iz,:,fine%it,:) = 0.75 * coarse%mem(1,1,cz,:,coarse%it,:) + &
                                           0.25 * coarse%mem(1,1,oz,:,coarse%it,:)
          enddo
        else ! 2D run in x-z plane
          do iz=fine%li(3), fine%ui(3)
            cz = (iz-1)/2 + offset(3)
            oz = cz + merge(-1,1,modulo(iz-fine%li(3),2)==0)
          do ix=fine%li(1), fine%ui(1)
            cx = (ix-1)/2 + offset(1)
            ox = cx + merge(-1,1,modulo(ix-fine%li(1),2)==0)
            fine%mem(ix,1,iz,:,fine%it,:) = 0.75*0.75 * coarse%mem(cx,1,cz,:,coarse%it,:) + &
                                            0.25*0.75 * coarse%mem(ox,1,cz,:,coarse%it,:) + &
                                            0.75*0.25 * coarse%mem(cx,1,oz,:,coarse%it,:) + &
                                            0.25*0.25 * coarse%mem(ox,1,oz,:,coarse%it,:)
          enddo
          enddo
        endif
      else
        if (flat(1)) then ! 2D run in y-z plane
          do iz=fine%li(3), fine%ui(3)
            cz = (iz-1)/2 + offset(3)
            oz = cz + merge(-1,1,modulo(iz-fine%li(3),2)==0)
          do iy=fine%li(2), fine%ui(2)
            cy = (iy-1)/2 + offset(2)
            oy = cy + merge(-1,1,modulo(iy-fine%li(2),2)==0)
            fine%mem(1,iy,iz,:,fine%it,:) = 0.75*0.75 * coarse%mem(1,cy,cz,:,coarse%it,:) + &
                                            0.25*0.75 * coarse%mem(1,oy,cz,:,coarse%it,:) + &
                                            0.75*0.25 * coarse%mem(1,cy,oz,:,coarse%it,:) + &
                                            0.25*0.25 * coarse%mem(1,oy,oz,:,coarse%it,:)
          enddo
          enddo
        else ! 3D run
          do iz=fine%li(3), fine%ui(3)
            cz = (iz-1)/2 + offset(3)
            oz = cz + merge(-1,1,modulo(iz-fine%li(3),2)==0)
          do iy=fine%li(2), fine%ui(2)
            cy = (iy-1)/2 + offset(2)
            oy = cy + merge(-1,1,modulo(iy-fine%li(2),2)==0)
          do ix=fine%li(1), fine%ui(1)
            cx = (ix-1)/2 + offset(1)
            ox = cx + merge(-1,1,modulo(ix-fine%li(1),2)==0)
            fine%mem(ix,iy,iz,:,fine%it,:) = 0.75*0.75*0.75 * coarse%mem(cx,cy,cz,:,coarse%it,:) + &
                                             0.25*0.75*0.75 * coarse%mem(ox,cy,cz,:,coarse%it,:) + &
                                             0.75*0.25*0.75 * coarse%mem(cx,oy,cz,:,coarse%it,:) + &
                                             0.25*0.25*0.75 * coarse%mem(ox,oy,cz,:,coarse%it,:) + &
                                             0.75*0.75*0.25 * coarse%mem(cx,cy,oz,:,coarse%it,:) + &
                                             0.25*0.75*0.25 * coarse%mem(ox,cy,oz,:,coarse%it,:) + &
                                             0.75*0.25*0.25 * coarse%mem(cx,oy,oz,:,coarse%it,:) + &
                                             0.25*0.25*0.25 * coarse%mem(ox,oy,oz,:,coarse%it,:)
          enddo
          enddo
          enddo
        endif
      endif
    endif
  endif

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
  call trace_begin ('interpolation_t%restrict', 3)
  if (first_time) then
    first_time = .false.
  end if
 
  if (guard) then
    ! test for overlaps and check where data has been uploaded already
  endif

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
  call trace_end
END SUBROUTINE restrict
!
END MODULE interpolation_mod
