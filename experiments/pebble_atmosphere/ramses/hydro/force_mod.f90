!> =============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!> =============================================================================
MODULE force_mod
  USE io_mod
  USE trace_mod
  !USE vector_mod
  USE mesh_mod
  USE random_mod
  USE planet_mod
  USE star_mod
  USE scaling_mod
  implicit none
  private
  real(8), parameter:: pi2=8.0*atan(1.0)
  type, public:: force_t
    real:: omega(3)
    real:: frc = 1.0
    procedure(star_and_planet), pointer:: selected=>null()
  contains
    procedure :: init
    procedure :: friction
    procedure :: friction2
  end type
 logical, save:: first_time=.true.
CONTAINS

!> =============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!> =============================================================================
SUBROUTINE init (self)
   class(force_t):: self
   integer:: iostat
   real   :: artificial_friction_time
   character(len=32), save:: type='star_and_planet'
   namelist /force_params/ type
   !----------------------------------------------------------------------------
   call trace_begin('force_t%init')
   !$omp critical (input_cr)
   if (first_time) then
     first_time = .false.
     rewind (io%input)
     read (io%input, force_params, iostat=iostat)
     if (io%master) write (*, force_params)
   end if
   !$omp end critical (input_cr)
   select case (type)
   case ('single_solenoidal')
     self%selected => single_solenoidal
   case ('star_and_planet')
     self%selected => star_and_planet
   case default
     nullify(self%selected)
   end select
   call trace_end
END SUBROUTINE init

!> =====================================================================
SUBROUTINE friction(self, time, m, ff)
  class(force_t) :: self
  real(8) :: time
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3):: ff
  ! --------------------------------------------------------------------
  if (time < planet%frc_time*2.) then
    ff = ff - ff*(exp(-(time/planet%frc_time)**2))
  end if
END SUBROUTINE friction

SUBROUTINE friction2(self, time)
  class(force_t) :: self
  real(8) :: time
  ! --------------------------------------------------------------------
  if (time < (planet%frc_time * 2)) then
    self%frc = 1.-(exp(-(time/planet%frc_time)**2))
  else
    self%frc = 1.0
  end if
END SUBROUTINE friction2

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
  integer:: ix, iy, iz, iv
  real(8):: kx, ky, kz, rs(3), rp(3), rs1, rp1
  real(8):: x, y, z, xp, yp, zp, ds2(3), R_Hill
  !-----------------------------------------------------------------------------
  call trace_begin('force_t%star_and_planet')
  call self%friction2(time)
  ds2 = 0.5*m%d
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + r3(iz)
    zp = z - planet%position(3)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + r2(iy)
      yp = y - planet%position(2)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + r1(ix)
        xp = x - planet%position(1)
        !-----------------------------------------------------------------------
        ! The force from the planet is included explicitly, and is oriented
        ! away from the planet, which is at the center of the local coordinates.
        ! The L1 and L2 Lagrange points are the result of net centrifugal
        ! forces that point in the opposite direction to the force of gravity,
        ! cancelling at a distance where $\Omega^2 r$ is equal to $GM/r^2$
        !-----------------------------------------------------------------------
        rp1 = sqrt(xp**2+yp**2+zp**2)
        ff(ix,iy,iz,1) = -planet%mass*scaling%grav*xp/rp1**3
        ff(ix,iy,iz,2) = -planet%mass*scaling%grav*yp/rp1**3 
        ff(ix,iy,iz,3) = -planet%mass*scaling%grav*zp/rp1**3
        !-----------------------------------------------------------------------
        ! Add the vertical and radial rotating coordinate frame accelerations.
        ! The force in the radial direction increases with increasing y,
        ! giving an equilbrium velocity in the x-direction that is positive
        ! (faster rotation) for a positive y-displacement.  Vertically, there
        ! is a net force that comes from the projection of the force of
        ! gravity onto the z-axis, and has the same scaling constant as the
        ! radial force, but the opposite sign.
        !-----------------------------------------------------------------------
        ff(ix,iy,iz,2) = ff(ix,iy,iz,2) + 3.0*planet%omega(3)**2*yp  *self%frc
        ff(ix,iy,iz,3) = ff(ix,iy,iz,3) -     planet%omega(3)**2*zp  *self%frc
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Radial Hill distance = L2
  !-----------------------------------------------------------------------------
  if (io%verbose > 1) then
    R_Hill = (planet%mass*scaling%grav/(3.*planet%omega(3)**2))**0.33333333
    if (abs(m(1)%p)          < 0.5*m(1)%d .and. &
        abs(m(2)%p - R_Hill) < 0.5*m(2)%s .and. &
        abs(m(3)%p)          < 0.5*m(3)%d) then
      ix = m(1)%o
      iz = m(3)%o
      iy = m(2)%o + (R_Hill-m(2)%p)/m(2)%d
      x = m(1)%p + r1(ix) - planet%position(1)
      y = m(2)%p + r2(iy) - planet%position(2)
      z = m(3)%p + r3(iz) - planet%position(3)
      rp1 = sqrt(x**2+y**2+z**2)
      print *, 'radial   R_Hill, force:', ix,iy,iz, R_Hill, rp1, ff(ix,iy,iz,2)
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Vertical Hill distance
  !-----------------------------------------------------------------------------
  if (io%verbose > 1) then
    if (abs(m(1)%p)          < 0.5*m(1)%d .and. &
        abs(m(2)%p)          < 0.5*m(2)%d .and. &
        abs(m(3)%p - R_Hill) < 0.5*m(3)%s) then
      ix = m(1)%o
      iy = m(2)%o
      iz = m(3)%o + (R_Hill-m(3)%p)/m(3)%d
      x = m(1)%p + r1(ix) - planet%position(1)
      y = m(2)%p + r2(iy) - planet%position(2)
      z = m(3)%p + r3(iz) - planet%position(3)
      rp1 = sqrt(x**2+y**2+z**2)
      print *, 'vertical R_Hill, force:', ix,iy,iz, R_Hill, rp1, ff(ix,iy,iz,3)
    end if
  end if
  end associate
  !-----------------------------------------------------------------------------
  ! Add the classical Coriolis force term: -2*(\Omega x U). A positive omega(3)
  ! is motion in the positive x-direction, if y points towards the system center.
  ! Hence, a positive velocity in the x-direction increases the centripetal force,
  ! and leads to acceleration in the negative y-direction -- vice versa for
  ! motion in the y-direction; a motion inwards (+y) leads to acceleration in
  ! the forward (+x) direction
  !-----------------------------------------------------------------------------
  ff(:,:,:,1) = d*ff(:,:,:,1) - 2.*(planet%omega(2)*p(:,:,:,3)-planet%omega(3)*p(:,:,:,2)) *self%frc
  ff(:,:,:,2) = d*ff(:,:,:,2) - 2.*(planet%omega(3)*p(:,:,:,1)-planet%omega(1)*p(:,:,:,3)) *self%frc
  ff(:,:,:,3) = d*ff(:,:,:,3) - 2.*(planet%omega(1)*p(:,:,:,2)-planet%omega(2)*p(:,:,:,1)) *self%frc
 ! call self%friction(time, m, ff)
  call trace_end
END FUNCTION star_and_planet


!> =============================================================================
!> Forcing on a single 3-D wavenumber, changing amplitude and phase periodically
!> =============================================================================
FUNCTION single_solenoidal (self, time, d, p, m) RESULT (ff)
  class(force_t)        :: self
  real(8)               :: time
  real, dimension(:,:,:), pointer:: d
  real, dimension(:,:,:,:), pointer:: p
  class(mesh_t), pointer, dimension(:) :: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3) :: ff
  !.............................................................................
  type(random_t)        :: random
  integer               :: ix, iy, iz, iostat
  real(8)               :: kx, ky, kz
  !.............................................................................
  integer, save:: seed=-11
  real, save:: a(3)=0.0
  real, save:: t_turn=0.3, t_save=0.0, a0(3)=0.1, k(3)=1.0, phi(3)=0.0
  namelist /solenoidal_params/ k, a0, t_turn
  !-----------------------------------------------------------------------------
  call trace_begin('force_t%single_solenoidal')
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, solenoidal_params, iostat=iostat)
    if (io%master) write (*, solenoidal_params)
  end if
  !$omp end critical (input_cr)
  if (time >= t_save) then
    a = a0*random%ran3()
    t_save = t_save + t_turn
!     if (io%verbose>3) print *,'forcing', a, phi
  end if
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    kz = pi2*k(3)*(r3(iz)+m(3)%p-0.5d0*m(3)%s)/m(3)%b + phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky = pi2*k(2)*(r2(iy)+m(2)%p-0.5d0*m(2)%s)/m(2)%b + phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx = pi2*k(1)*(r1(ix)+m(1)%p-0.5d0*m(1)%s)/m(1)%b + phi(1)
        ff(ix,iy,iz,1) = a(1)*(sin(ky)+sin(kz))
        ff(ix,iy,iz,2) = a(2)*(sin(kz)+sin(kx))
        ff(ix,iy,iz,3) = a(3)*(sin(kx)+sin(ky))
      end do
    end do
  end do
  end associate
  call trace_end
END FUNCTION single_solenoidal

END MODULE
