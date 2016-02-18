!> =============================================================================
!> Module to assist in handling boundaries in the Poisson solver.  If phi is
!> imposed at and outside some boundary, then the error should be ignored at
!> the very boundary, s
!> =============================================================================
MODULE boundary_mod
  USE trace_mod
  USE patches_mod
  implicit none
  private
  type, public:: boundary_t
    real:: r_BE
  contains
    procedure:: sphere
  end type
  type(boundary_t), public:: boundary
CONTAINS
!> =============================================================================
!> Enforce constant values outside a radius r=r_BE (not equalt to the BE-radius)
!> =============================================================================
SUBROUTINE sphere (self, patch, f, value, extra)
  class(boundary_t):: self
  class(patch_t):: patch
  integer:: i, j, k, ii(3), iv
  real:: x, y, z, r, ntot, nout, value, extra
  real(kind=8), dimension(patch%gn(1),patch%gn(2),patch%gn(3)) :: f
associate (gn=>patch%gn, ds=>patch%ds, pos=>patch%position)
  !.............................................................................
  call trace_begin('boundary_mod::boundary_t%bdry_phi', 4)
  do k=1,gn(3)
  do j=1,gn(2)
  do i=1,gn(1)
    x = (i-(gn(1)+1)*0.5)*ds(1) + pos(1)
    y = (j-(gn(2)+1)*0.5)*ds(2) + pos(2)
    z = (k-(gn(3)+1)*0.5)*ds(3) + pos(3)
    r = sqrt(x**2+y**2+z**2)
    if (r+extra > self%r_BE) then
      f(i,j,k) = value
    end if
  end do
  end do
  end do
  call trace_end
end associate
END SUBROUTINE
END MODULE boundary_mod
