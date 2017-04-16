!> =============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!> =============================================================================
MODULE force_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE planet_mod
  USE scaling_mod
  implicit none
  private
  type, public:: force_t
    procedure(star_and_planet), pointer:: selected=>null()
  contains
    procedure :: init
  end type
 logical, save:: first_time=.true.
CONTAINS

!> =============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!> =============================================================================
SUBROUTINE init (self, default)
   class(force_t):: self
   character(len=*), optional:: default
   integer:: iostat
   character(len=20), save:: type='planet_atmosphere'
   namelist /force_params/ type
   !----------------------------------------------------------------------------
   call trace_begin('force_t%init')
   !$omp critical (input_cr)
   if (first_time) then
     first_time = .false.
     if (present(default)) type = default
     rewind (io%input)
     read (io%input, force_params, iostat=iostat)
     if (io%master) write (*, force_params)
   end if
   !$omp end critical (input_cr)
   select case (trim(type))
   case ('planet_atmosphere')
     self%selected => star_and_planet
   case ('disk')
     self%selected => disk
   case default
     nullify(self%selected)
   end select
   call trace_end
END SUBROUTINE init

!> =============================================================================
!> Force of gravity from a star and a planet.  To interprete the local forces,
!> we choose a position such that the star is in the positive y-direction, and
!> hence a positive Omega_z corresponds to positive rotation speed in the x-
!> direction
!> =============================================================================
FUNCTION star_and_planet (self, time, d, p, m) RESULT (ff)
  class(force_t):: self
  real(8):: time
  real, dimension(:,:,:), pointer:: d
  real, dimension(:,:,:,:), pointer:: p
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3):: ff
  !.............................................................................
  integer:: ix, iy, iz
  real(8):: xp, yp, zp, rp
  real:: ux, uy, newton
  !-----------------------------------------------------------------------------
  call trace_begin('force_t%star_and_planet')
  ff = 0.0
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  newton = planet%mass*scaling%grav
  !-----------------------------------------------------------------------------
  ! To allow vectorization there should be no ifs in this loop
  !-----------------------------------------------------------------------------
  do iz=m(3)%lb,m(3)%ub
    zp = m(3)%p + r3(iz) - planet%position(3)
    do iy=m(2)%lb,m(2)%ub
      yp = m(2)%p + r2(iy) - planet%position(2)
      do ix=m(1)%lb,m(1)%ub
        xp = m(1)%p + r1(ix) - planet%position(1)
        !-----------------------------------------------------------------------
        ! The force from the planet is included explicitly, and is oriented
        ! away from the planet, which is at the center of the local coordinates.
        ! The L1 and L2 Lagrange points are the result of net centrifugal
        ! forces that point in the opposite direction to the force of gravity,
        ! cancelling at a distance where $\Omega^2 r$ is equal to $GM/r^2$
        !-----------------------------------------------------------------------
        rp = sqrt(xp**2+yp**2+zp**2)
        rp = max(rp,planet%radius*0.1)
        ff(ix,iy,iz,1) = -newton*xp/rp**3
        ff(ix,iy,iz,2) = -newton*yp/rp**3
        ff(ix,iy,iz,3) = -newton*zp/rp**3
        !-----------------------------------------------------------------------
        ! Add the vertical and radial rotating coordinate frame accelerations.
        ! The net force in the radial direction increases with increasing x,
        ! giving that x now points in the radial direction.  It has two
        ! contributions:  The centripetal force \Omega^2 x, and the reduction
        ! of the force of gravity from the central star, which contributs a
        ! factor of 2 more, for a total factor of 3.  Vertically, there
        ! is a net force that comes from the projection of the force of
        ! gravity onto the z-axis, and has the same scaling constant as the
        ! radial force, but the opposite sign.
        !-----------------------------------------------------------------------
        ff(ix,iy,iz,1) = ff(ix,iy,iz,1) + 3.0*planet%omega(3)**2*xp
        ff(ix,iy,iz,3) = ff(ix,iy,iz,3) -     planet%omega(3)**2*zp
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Add the classical Coriolis force term: -2*(\Omega x U). A positive omega(3)
  ! corresponds to rotation in the y-direction, if x points radially outwards.
  ! Hence, a positive velocity in the y-direction increases the centripetal force,
  ! and leads to acceleration in the positive x-direction -- vice versa for
  ! motion in the x-direction; a motion inwards (-x) leads to acceleration in
  ! the forward (+y) direction
  !-----------------------------------------------------------------------------
  ff(:,:,:,1) = d*ff(:,:,:,1) - 2.*(planet%omega(2)*p(:,:,:,3)-planet%omega(3)*p(:,:,:,2))
  ff(:,:,:,2) = d*ff(:,:,:,2) - 2.*(planet%omega(3)*p(:,:,:,1)-planet%omega(1)*p(:,:,:,3))
  ff(:,:,:,3) = d*ff(:,:,:,3) - 2.*(planet%omega(1)*p(:,:,:,2)-planet%omega(2)*p(:,:,:,1))
    end associate
  call trace_end
END FUNCTION star_and_planet

!> =============================================================================
!> Force of gravity from a star, in a rotating disk coordinate system
!> =============================================================================
FUNCTION disk (self, time, d, p, m) RESULT (ff)
  class(force_t):: self
  real(8):: time
  real, dimension(:,:,:), pointer:: d
  real, dimension(:,:,:,:), pointer:: p
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3):: ff
  !.............................................................................
  integer:: ix, iy, iz
  real(8):: xp, yp, zp
  real:: ux, uy
  !-----------------------------------------------------------------------------
  call trace_begin('force_t%disk')
  ff = 0.0
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  !-----------------------------------------------------------------------------
  ! To allow vectorization there should be no ifs in this loop
  !-----------------------------------------------------------------------------
  do iz=m(3)%lb,m(3)%ub
    zp = m(3)%p + r3(iz) - planet%position(3)
    do iy=m(2)%lb,m(2)%ub
      yp = m(2)%p + r2(iy) - planet%position(2)
      do ix=m(1)%lb,m(1)%ub
        xp = m(1)%p + r1(ix) - planet%position(1)
        !-----------------------------------------------------------------------
        ! Add the vertical and radial rotating coordinate frame accelerations.
        ! The net force in the radial direction increases with increasing x,
        ! giving that x now points in the radial direction.  It has two
        ! contributions:  The centripetal force \Omega^2 x, and the reduction
        ! of the force of gravity from the central star, which contributs a
        ! factor of 2 more, for a total factor of 3.  Vertically, there
        ! is a net force that comes from the projection of the force of
        ! gravity onto the z-axis, and has the same scaling constant as the
        ! radial force, but the opposite sign.
        !-----------------------------------------------------------------------
        ff(ix,iy,iz,1) = ff(ix,iy,iz,1) + 3.0*planet%omega(3)**2*xp
        ff(ix,iy,iz,3) = ff(ix,iy,iz,3) -     planet%omega(3)**2*zp
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Add the classical Coriolis force term: -2*(\Omega x U). A positive omega(3)
  ! is motion in the positive x-direction, if y points towards the system center.
  ! Hence, a positive velocity in the x-direction increases the centripetal force,
  ! and leads to acceleration in the negative y-direction -- vice versa for
  ! motion in the y-direction; a motion inwards (+y) leads to acceleration in
  ! the forward (+x) direction
  !-----------------------------------------------------------------------------
  ff(:,:,:,1) = d*ff(:,:,:,1) - 2.*(planet%omega(2)*p(:,:,:,3)-planet%omega(3)*p(:,:,:,2))
  ff(:,:,:,2) = d*ff(:,:,:,2) - 2.*(planet%omega(3)*p(:,:,:,1)-planet%omega(1)*p(:,:,:,3))
  ff(:,:,:,3) = d*ff(:,:,:,3) - 2.*(planet%omega(1)*p(:,:,:,2)-planet%omega(2)*p(:,:,:,1))
  !-----------------------------------------------------------------------------
  end associate
  call trace_end
END FUNCTION disk

END MODULE
