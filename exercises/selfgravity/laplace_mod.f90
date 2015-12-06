!> =============================================================================
!> For now, a placeholder...   The operators are from the PP-code
!> =============================================================================
MODULE laplace_mod
  USE io_mod
  USE trace_mod
  USE patches_mod
  implicit none
  private
CONTAINS
!=======================================================================
! 7pt 6th order Laplace operator, for filtering high-k
!=======================================================================
FUNCTION laplace6 (self, f) result(h)
  class(patch_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)) :: f, h
  integer ix, iy, iz, lb(3), ub(3)
  real lpm0, lpm1, lpm2, lpm3, odx6, ody6, odz6
!.......................................................................
  lpm3 = + 1.
  lpm2 =  -6.
  lpm1 = +15.
  lpm0 = -20.
  odx6 = (1./self%ds(1))**6
  ody6 = (1./self%ds(2))**6
  odz6 = (1./self%ds(3))**6

  lb = self%li
  ub = self%ui
  do iz=lb(3),ub(3)
  do iy=lb(2),ub(2)
  do ix=lb(1),ub(1)
    h(ix,iy,iz) = ( lpm3*(f(ix-3,iy,iz)+f(ix+3,iy,iz)) &
                  + lpm2*(f(ix-2,iy,iz)+f(ix+2,iy,iz)) &
                  + lpm1*(f(ix-1,iy,iz)+f(ix+1,iy,iz)) &
                  + lpm0* f(ix  ,iy,iz) ) * odx6       &
                + ( lpm3*(f(ix,iy-3,iz)+f(ix,iy+3,iz)) &
                  + lpm2*(f(ix,iy-2,iz)+f(ix,iy+2,iz)) &
                  + lpm1*(f(ix,iy-1,iz)+f(ix,iy+1,iz)) &
                  + lpm0* f(ix,iy,  iz) ) * ody6       &
                + ( lpm3*(f(ix,iy,iz-3)+f(ix,iy,iz+3)) &
                  + lpm2*(f(ix,iy,iz-2)+f(ix,iy,iz+2)) &
                  + lpm1*(f(ix,iy,iz-1)+f(ix,iy,iz+1)) &
                  + lpm0* f(ix,iy,iz  ) ) * odz6
  enddo
  enddo
  enddo
END FUNCTION laplace6
!=======================================================================
! 7pt 2nd order Laplace operator, for filtering high-k
!=======================================================================
FUNCTION laplace_6th (self, f) result(h)
  class(patch_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)) :: f, h
  integer ix, iy, iz, lb(3), ub(3)
  real lpm0, lpm1, lpm2, lpm3, odxq, odyq, odzq
!.......................................................................
  lpm3 =  + 2./180.
  lpm2 =  -27./180.
  lpm1 = +270./180.
  lpm0 = -490./180.
  odxq = (1./self%ds(1))**2
  odyq = (1./self%ds(2))**2
  odzq = (1./self%ds(3))**2
  lb = self%li
  ub = self%ui
  do iz=lb(3),ub(3)
  do iy=lb(2),ub(2)
  do ix=lb(1),ub(1)
    h(ix,iy,iz) = ( lpm3*(f(ix-3,iy,iz)+f(ix+3,iy,iz)) &
                  + lpm2*(f(ix-2,iy,iz)+f(ix+2,iy,iz)) &
                  + lpm1*(f(ix-1,iy,iz)+f(ix+1,iy,iz)) &
                  + lpm0* f(ix  ,iy,iz) ) * odxq       &
                + ( lpm3*(f(ix,iy-3,iz)+f(ix,iy+3,iz)) &
                  + lpm2*(f(ix,iy-2,iz)+f(ix,iy+2,iz)) &
                  + lpm1*(f(ix,iy-1,iz)+f(ix,iy+1,iz)) &
                  + lpm0* f(ix,iy,  iz) ) * odyq       &
                + ( lpm3*(f(ix,iy,iz-3)+f(ix,iy,iz+3)) &
                  + lpm2*(f(ix,iy,iz-2)+f(ix,iy,iz+2)) &
                  + lpm1*(f(ix,iy,iz-1)+f(ix,iy,iz+1)) &
                  + lpm0* f(ix,iy,iz  ) ) * odzq
  enddo
  enddo
  enddo
END FUNCTION laplace_6th

!=======================================================================
! 5pt 2nd order Laplace operator, for filtering high-k
!=======================================================================
FUNCTION laplace_4th (self, f) result(h)
  class(patch_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)) :: f, h
  integer ix, iy, iz, lb(3), ub(3)
  real lpm0, lpm1, lpm2, odxq, odyq, odzq
!.......................................................................
  lpm2 =  -1./12.
  lpm1 = +16./12.
  lpm0 = -30./12.
  odxq = (1./self%ds(1))**2
  odyq = (1./self%ds(2))**2
  odzq = (1./self%ds(3))**2
  lb = self%li
  ub = self%ui
  do iz=lb(3),ub(3)
  do iy=lb(2),ub(2)
  do ix=lb(1),ub(1)
    h(ix,iy,iz) = ( lpm2*(f(ix-2,iy,iz)+f(ix+2,iy,iz)) &
                  + lpm1*(f(ix-1,iy,iz)+f(ix+1,iy,iz)) &
                  + lpm0* f(ix  ,iy,iz) ) * odxq       &
                + ( lpm2*(f(ix,iy-2,iz)+f(ix,iy+2,iz)) &
                  + lpm1*(f(ix,iy-1,iz)+f(ix,iy+1,iz)) &
                  + lpm0* f(ix,iy,  iz) ) * odyq       &
                + ( lpm2*(f(ix,iy,iz-2)+f(ix,iy,iz+2)) &
                  + lpm1*(f(ix,iy,iz-1)+f(ix,iy,iz+1)) &
                  + lpm0* f(ix,iy,iz  ) ) * odzq
  enddo
  enddo
  enddo
END FUNCTION laplace_4th

!=======================================================================
! 3pt 2nd order Laplace operator, for filtering high-k
!=======================================================================
FUNCTION laplace_2nd (self, f) result(h)
  class(patch_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)) :: f, h
  integer ix, iy, iz, lb(3), ub(3)
  real lpm0, lpm1, lpm2, odxq, odyq, odzq
!.......................................................................
  lpm1 = +1.0
  lpm0 = -2.0
  odxq = (1./self%ds(1))**2
  odyq = (1./self%ds(2))**2
  odzq = (1./self%ds(3))**2
  lb = self%li
  ub = self%ui
  do iz=lb(3),ub(3)
  do iy=lb(2),ub(2)
  do ix=lb(1),ub(1)
    h(ix,iy,iz) = ( lpm1*(f(ix-1,iy,iz)+f(ix+1,iy,iz)) &
                  + lpm0* f(ix  ,iy,iz) ) * odxq       &
                + ( lpm1*(f(ix,iy-1,iz)+f(ix,iy+1,iz)) &
                  + lpm0* f(ix,iy,  iz) ) * odyq       &
                + ( lpm1*(f(ix,iy,iz-1)+f(ix,iy,iz+1)) &
                  + lpm0* f(ix,iy,iz  ) ) * odzq
  enddo
  enddo
  enddo
END FUNCTION laplace_2nd
!=======================================================================
END MODULE laplace_mod
