!===============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!===============================================================================
MODULE force_mod
  USE io_mod
  USE trace_mod
  USE vector_mod
  USE mesh_mod
  USE random_mod
  USE planet_mod
  USE star_mod
  USE scaling_mod
  USE stagger_mod
  implicit none
  private
  real(8), parameter:: pi2=8.0*atan(1.0)
  type, public:: force_t
    real:: omega(3)
    procedure(star_and_planet), pointer:: selected=>null()
  contains
    procedure init
  end type
 logical, save:: first_time=.true.
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self)
   class(force_t):: self
   integer:: iostat
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
     if (first_time) print*,'UNKNOWN FORCING TYPE'
   end select
   call trace_end
END SUBROUTINE init

!===============================================================================
!> Force of gravity from a star and a planet.  To interprete the local forces,
!> we choose a position such that the star is in the positive y-direction, and
!> hence a positive Omega_z corresponds to positive rotation speed in the x-
!> direction
!===============================================================================
FUNCTION star_and_planet (self, time, d, p, Ux, Uy, Uz, m) RESULT (ff)
  class(force_t):: self
  real(8):: time
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3):: ff
  real, dimension(:,:,:):: Ux, Uy, Uz, d
  real, dimension(:,:,:,:):: p
  !.............................................................................
  integer:: ix, iy, iz, iv
  real(8):: kx, ky, kz, rs(3), rp(3), rs1, rp1
  real(8):: x, y, z, xs, ys, zs, xp, yp, zp, ds2(3)
  !-----------------------------------------------------------------------------
  call trace_begin('force_t:single_solenoidal')
  ds2 = 0.5*m%d
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + r3(iz)
    zp = planet%position(3) - z
    zs = zp + ds2(3)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + r2(iy)
      yp = planet%position(2) - y
      ys = yp + ds2(2)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + r1(ix)
        xp = planet%position(1) - x
        xs = xp + ds2(1)
        !-----------------------------------------------------------------------
        ! The force from the planet is included explicitly.  Since the [xyz]p
        ! vector is twoards the planet, the staggered mesh adjustment is positive
        !-----------------------------------------------------------------------
        rp1 = sqrt(xs**2+yp**2+zp**2)
        ff(ix,iy,iz,1) = planet%mass*scaling%grav*xs/rp1**3
        rp1 = sqrt(xp**2+ys**2+zp**2)
        ff(ix,iy,iz,2) = planet%mass*scaling%grav*ys/rp1**3
        rp1 = sqrt(xp**2+yp**2+zs**2)
        ff(ix,iy,iz,3) = planet%mass*scaling%grav*zs/rp1**3
        !-----------------------------------------------------------------------
        ! Add the vertical and radial rotating coordinate frame accelerations.
        ! The force in the radial direction increases with increasing y,
        ! giving an equilbrium velocity in the x-direction that is positive
        ! (faster rotation) for a positive y-displacement
        !-----------------------------------------------------------------------
        ff(ix,iy,iz,3) = ff(ix,iy,iz,3) -     planet%omega(3)**2*zs
        ff(ix,iy,iz,2) = ff(ix,iy,iz,2) + 3.0*planet%omega(3)**2*ys
      end do
    end do
  end do
  !print *,'planet_force: mass, grav =', planet%mass, scaling%grav
  !print *,'planet_force: omega      =', planet%omega
  !print *,'planet_force: position   =', planet%position
  end associate
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

!===============================================================================
!> Forcing on a single 3-D wavenumber, changing amplitude and phase periodically
!===============================================================================
FUNCTION single_solenoidal (self, time, d, p, Ux, Uy, Uz, m) RESULT (ff)
  class(force_t)        :: self
  real(8)               :: time
  class(mesh_t), pointer, dimension(:) :: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3) :: ff
  real, dimension(:,:,:):: Ux, Uy, Uz, d
  real, dimension(:,:,:,:):: p
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
  call trace_begin('force_t:single_solenoidal')
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
