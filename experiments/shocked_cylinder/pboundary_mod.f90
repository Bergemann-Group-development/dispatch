!===============================================================================
!> Boundary conditions for a 2-D box
!===============================================================================
MODULE pboundary_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE mhd_mod
  USE ic_mod
  USE kinds_mod
  implicit none
  private

  type, public:: pboundary_t
  contains
    procedure conditions
  end type
  type(pboundary_t), public:: pboundary

CONTAINS

!===============================================================================
!> Impose boundary conditions, taking staggering into account.
!> Only for use with Cartesian patches.
!===============================================================================
SUBROUTINE conditions (self, initial, patch, it)
  class(pboundary_t) :: self
  class(ic_t) :: initial
  class(mhd_t), target :: patch
  class(mesh_t), pointer :: mx, my, mz
  integer :: ix, iy, iz, offi, offo, it
  real(4) :: x, y, z
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: q
  !-----------------------------------------------------------------------------
  if (patch%mesh_type /= mesh_types%Cartesian) return
  call trace_begin('pboundary_t%conditions')

  mx => patch%mesh(1)
  my => patch%mesh(2)
  mz => patch%mesh(3)
  q  => patch%mem(:,:,:,:,it,1)
  ! x-boundaries
  if (patch%n(1) > 1) then
    do iz=mz%lb,mz%ub
      do iy=my%lb,my%ub
        ! inflow boundary (lower x)
        do ix=mx%lb,mx%lo
          x = mx%centre_nat + patch%mesh(1)%r(ix)
          if ( x <= (mx%p - 0.5*mx%s) )then
            q(ix,iy,iz,patch%idx%d ) = initial%rhol
            q(ix,iy,iz,patch%idx%e ) = initial%presl / (patch%gamma - 1.0)
            q(ix,iy,iz,patch%idx%px) = initial%vxl * initial%rhol
            q(ix,iy,iz,patch%idx%py) = initial%vyl * initial%rhol
            q(ix,iy,iz,patch%idx%pz) = initial%vzl * initial%rhol
          end if
        end do
        ! outflow boundary (upper x)
        do ix=mx%uo,mx%ub
          x = mx%centre_nat + patch%mesh(1)%r(ix)
          if ( x >= (mx%p + 0.5*mx%s) )then
            q(ix,iy,iz,patch%idx%d ) = q(mx%ui-1,iy,iz,patch%idx%d )
            q(ix,iy,iz,patch%idx%e ) = q(mx%ui-1,iy,iz,patch%idx%e )
            q(ix,iy,iz,patch%idx%px) = q(mx%ui-1,iy,iz,patch%idx%px)
            q(ix,iy,iz,patch%idx%py) = q(mx%ui-1,iy,iz,patch%idx%py)
            q(ix,iy,iz,patch%idx%pz) = q(mx%ui-1,iy,iz,patch%idx%pz)
          end if
        end do
      end do
    end do
  end if

  ! y-boundaries; reflection boundaries
  if (patch%n(2) > 1) then
    offi = 2*my%lo + 1
    offo = 2*my%uo - 1
    do iz=mz%lb,mz%ub
      do ix=mx%lb,mx%ub
        ! (lower y)
        do iy=my%lo,my%lb,-1
          q(ix,iy,iz,patch%idx%d ) =  q(ix,offi-iy,iz,patch%idx%d )
          q(ix,iy,iz,patch%idx%e ) =  q(ix,offi-iy,iz,patch%idx%e )
          q(ix,iy,iz,patch%idx%px) =  q(ix,offi-iy,iz,patch%idx%px)
          q(ix,iy,iz,patch%idx%py) = -q(ix,offi-iy,iz,patch%idx%py)
          q(ix,iy,iz,patch%idx%pz) =  q(ix,offi-iy,iz,patch%idx%pz)
        end do
        ! (upper y)
        do iy=my%uo,my%ub,+1
          q(ix,iy,iz,patch%idx%d ) =  q(ix,offo-iy,iz,patch%idx%d )
          q(ix,iy,iz,patch%idx%e ) =  q(ix,offo-iy,iz,patch%idx%e )
          q(ix,iy,iz,patch%idx%px) =  q(ix,offo-iy,iz,patch%idx%px)
          q(ix,iy,iz,patch%idx%py) = -q(ix,offo-iy,iz,patch%idx%py)
          q(ix,iy,iz,patch%idx%pz) =  q(ix,offo-iy,iz,patch%idx%pz)
        end do
      end do
    end do
  end if

  ! periodic in z (is this necessary?)

  nullify (mx, my, mz)
  call trace_end
END SUBROUTINE conditions

END MODULE pboundary_mod
