!> =============================================================================
!> Direct force of gravity
!> =============================================================================
MODULE force_mod
  USE io_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  USE mesh_mod
  USE planet_mod
  USE star_mod
  USE stagger_mod
  implicit none
  private
  type, public:: force_t
    real:: omega(3)
    procedure(star_and_planet), pointer:: selected=>null()
  contains
    procedure init
  end type
  type(force_t), public:: force
CONTAINS

!> =============================================================================
!> Select type of force -- could be more than one planet, for example
!> =============================================================================
SUBROUTINE init (self)
  class(force_t):: self
  self%selected => star_and_planet
END SUBROUTINE

!> =============================================================================
!> Force of gravity from a star and a planet
!> =============================================================================
FUNCTION star_and_planet (self, time, Ux, Uy, Uz, m) RESULT (ff)
  class(force_t):: self
  real(8):: time
  type(mesh_t):: m(3)
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3):: ff
  real, dimension(:,:,:):: Ux, Uy, Uz
  !.............................................................................
  integer:: ix, iy, iz, iv
  real(8):: kx, ky, kz, rs(3), rp(3), rs1, rp1
  real(8):: x, y, z, xs, ys, zs, xp, yp, zp, ds2(3)
  !-----------------------------------------------------------------------------
  call trace_begin('force_t:single_solenoidal')
  ds2 = 0.5*m%d
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + m(3)%r(iz)
    zp = planet%position(3) - z
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + m(2)%r(iy)
      yp = planet%position(2) - y
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + m(1)%r(ix)
        xp = planet%position(1) - x
        rp1 = sqrt((xp+ds2(1))**2+yp**2+zp**2)
        ff(ix,iy,iz,1) = planet%mass*scaling%grav*(xp+ds2(1))/rp1**3
        rp1 = sqrt(xp**2+(yp+ds2(2))**2+zp**2)
        ff(ix,iy,iz,2) = planet%mass*scaling%grav*(yp+ds2(2))/rp1**3
        rp1 = sqrt(xp**2+yp**2+(zp-ds2(3))**2)
        ff(ix,iy,iz,3) = planet%mass*scaling%grav*(zp+ds2(3))/rp1**3
        !-----------------------------------------------------------------------
        ! Add the vertical and radial rotating coordinate frame accelerations
        !-----------------------------------------------------------------------
        ff(ix,iy,iz,3) = ff(ix,iy,iz,3) - self%omega(3)**2*(zp-m(3)%d)
        ff(ix,iy,iz,2) = ff(ix,iy,iz,2) + self%omega(3)**2*(yp-m(2)%d)
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Add the classical Coriolis force terms, correctly centered (not important)
  !-----------------------------------------------------------------------------
  ff(:,:,:,1) = ff(:,:,:,1) + 2.*self%omega(3)*xdn(yup(Uy))
  ff(:,:,:,2) = ff(:,:,:,2) - 2.*self%omega(3)*ydn(xup(Ux))
  call trace_end
END FUNCTION star_and_planet

END MODULE
