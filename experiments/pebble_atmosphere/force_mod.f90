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
!> Force of gravity from a star and a planet.  To interprete the local forces,
!> we choose a position such that the star is in the positive y-direction, and
!> hence a positive Omega_z corresponds to positive rotation speed in the x-
!> direction
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
        !-----------------------------------------------------------------------
        ! The force from the planet is included explicitly.  Since the [xyz]p
        ! vector is twoards the planet, the staggered mesh adjustment is positive
        !-----------------------------------------------------------------------
        rp1 = sqrt((xp+ds2(1))**2+yp**2+zp**2)
        ff(ix,iy,iz,1) = planet%mass*scaling%grav*(xp+ds2(1))/rp1**3
        rp1 = sqrt(xp**2+(yp+ds2(2))**2+zp**2)
        ff(ix,iy,iz,2) = planet%mass*scaling%grav*(yp+ds2(2))/rp1**3
        rp1 = sqrt(xp**2+yp**2+(zp-ds2(3))**2)
        ff(ix,iy,iz,3) = planet%mass*scaling%grav*(zp+ds2(3))/rp1**3
        !-----------------------------------------------------------------------
        ! Add the vertical and radial rotating coordinate frame accelerations.
        ! The force in the radial direction increases with increasing y, 
        ! giving an equilbrium velocity in the x-direction that is positive
        ! (faster rotation) for a positive y-displacement
        !-----------------------------------------------------------------------
        ff(ix,iy,iz,3) = ff(ix,iy,iz,3) - self%omega(3)**2*(zp+ds2(3))
        ff(ix,iy,iz,2) = ff(ix,iy,iz,2) + self%omega(3)**2*(yp+ds2(2))
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Add the classical Coriolis force terms, correctly centered (not important),
  ! assuming (as above) that the motion is in the x-y plane.  A positive omega(3)
  ! is motion in the negative x-direction, if y points in the radial direction.
  ! Hence, a positive velocity in the x-direction reduces the centripetal force,
  ! and leads to acceleration in the negative y-direction -- vice versa for 
  ! motion in the y-direction
  !-----------------------------------------------------------------------------
  ff(:,:,:,1) = ff(:,:,:,1) + 2.*self%omega(3)*xdn(yup(Uy))
  ff(:,:,:,2) = ff(:,:,:,2) - 2.*self%omega(3)*ydn(xup(Ux))
  call trace_end
END FUNCTION star_and_planet

END MODULE
