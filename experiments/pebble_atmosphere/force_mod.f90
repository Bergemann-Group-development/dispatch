!===============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!===============================================================================
MODULE force_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE random_mod
  USE planet_mod
  USE scaling_mod
  implicit none
  private
  real(8), parameter:: pi2=8.0*atan(1.0)
  type, public:: force_t
    integer:: id
    character(len=64):: solver
    procedure(star_and_planet), pointer:: selected=>null()
  contains
    procedure init
  end type
 logical, save:: first_time=.true.
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self, solver, id, mesh)
   class(force_t):: self
   integer:: id
   class(mesh_t), pointer:: mesh(:)
   character(len=64):: solver
   integer:: iostat
   character(len=20), save:: type='planet_atmosphere'
   namelist /force_params/ type
   !----------------------------------------------------------------------------
   call trace_begin('force_t%init')
   call io%lock%set
   if (first_time) then
     first_time = .false.
     rewind (io%input)
     read (io%input, force_params, iostat=iostat)
     if (io%master) write (*, force_params)
   end if
   call io%lock%unset
   select case (trim(type))
   case ('planet_atmosphere')
     self%selected => star_and_planet
   case ('disk')
     self%selected => disk
   case default
     nullify(self%selected)
   end select
   self%solver = trim(solver)
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
  real, dimension(:,:,:), pointer:: Ux, Uy, Uz, d
  real, dimension(:,:,:,:), pointer:: p
  !.............................................................................
  ff = disk_perhaps_with_planet (self, time, d, p, Ux, Uy, Uz, m, with_planet=.true.)
END FUNCTION star_and_planet

!===============================================================================
!> Force of gravity from a star and a planet.  To interprete the local forces,
!> we choose a position such that the star is in the positive y-direction, and
!> hence a positive Omega_z corresponds to positive rotation speed in the x-
!> direction
!===============================================================================
FUNCTION disk (self, time, d, p, Ux, Uy, Uz, m) RESULT (ff)
  class(force_t):: self
  real(8):: time
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3):: ff
  real, dimension(:,:,:), pointer:: Ux, Uy, Uz, d
  real, dimension(:,:,:,:), pointer:: p
  !.............................................................................
  ff = disk_perhaps_with_planet (self, time, d, p, Ux, Uy, Uz, m)
END FUNCTION disk

!===============================================================================
!> Force of gravity from a star and a planet.  To interprete the local forces,
!> we choose a position such that the star is in the positive y-direction, and
!> hence a positive Omega_z corresponds to positive rotation speed in the x-
!> direction
!===============================================================================
FUNCTION disk_perhaps_with_planet (self, time, d, p, Ux, Uy, Uz, m, with_planet) RESULT (ff)
  class(force_t):: self
  real(8):: time
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3):: ff
  real, dimension(:,:,:), pointer:: Ux, Uy, Uz, d
  real, dimension(:,:,:,:), pointer:: p
  logical, optional:: with_planet
  !.............................................................................
  integer:: ix, iy, iz, iv
  real(8):: rp, newton
  real(8):: x, y, z, xs, ys, zs, xp, yp, zp, ds2(3), uxx, uyy
  !-----------------------------------------------------------------------------
  call trace_begin('force_t%disk_perhaps_with_planet')
  if (self%solver(1:7)=='ramses_') then
    ds2 = 0.
  else if (self%solver(1:7)=='stagger') then
    ds2 = 0.5*m%d
  else if (self%solver(1:6)=='zeus3d') then
    ds2 = 0.5*m%d
  else
    error stop 'force'
  end if
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)

  if (present(with_planet)) then
    newton = planet%mass*scaling%grav
  else
    newton = 0.0
  end if
  ff = 0.0
  if (planet%r_smooth <=0.0) then
    do iz=m(3)%lb,m(3)%ub
      z = m(3)%p + r3(iz)
      zp = z  - planet%position(3)
      zs = zp - ds2(3)
      do iy=m(2)%lb,m(2)%ub
        y = m(2)%p + r2(iy)
        yp = y  - planet%position(2)
        ys = yp - ds2(2)
        do ix=m(1)%lb,m(1)%ub
          x = m(1)%p + r1(ix)
          xp = x  - planet%position(1)
          xs = xp - ds2(1)
          !---------------------------------------------------------------
          ! The force from the planet is included explicitly, and is oriented
          ! away from the planet, which is at the center of the local coordinates.
          ! The L1 and L2 Lagrange points are the result of net centrifugal
          ! forces that point in the opposite direction to the force of gravity,
          ! cancelling at a distance where $\Omega^2 r$ is equal to $GM/r^2$
          !---------------------------------------------------------------
          rp = sqrt(xs**2+yp**2+zp**2)
          rp = max(rp,planet%radius*0.1)
          ff(ix,iy,iz,1) = -newton*xs/rp**3
          rp = sqrt(xp**2+ys**2+zp**2)
          rp = max(rp,planet%radius*0.1)
          ff(ix,iy,iz,2) = -newton*ys/rp**3
          rp = sqrt(xp**2+yp**2+zs**2)
          rp = max(rp,planet%radius*0.1)
          ff(ix,iy,iz,3) = -newton*zs/rp**3
          !---------------------------------------------------------------
          ! Add the vertical and radial rotating coordinate frame accelerations.
          ! The net force in the radial direction increases with increasing x,
          ! giving that x now points in the radial direction.  It has two
          ! contributions:  The centripetal force \Omega^2 x, and the reduction
          ! of the force of gravity from the central star, which contributes a
          ! factor of 2 more, for a total factor of 3.  Vertically, there
          ! is a net force that comes from the projection of the force of
          ! gravity onto the z-axis, and has the same scaling constant as the
          ! radial force, but the opposite sign.
          !---------------------------------------------------------------
          ff(ix,iy,iz,1) = ff(ix,iy,iz,1) + 3.0*planet%omega(3)**2*xs
          ff(ix,iy,iz,3) = ff(ix,iy,iz,3) -     planet%omega(3)**2*zs
        end do
      end do
    end do
  else
    do iz=m(3)%lb,m(3)%ub
      z = m(3)%p + r3(iz)
      zp = z  - planet%position(3)
      zs = zp - ds2(3)
      do iy=m(2)%lb,m(2)%ub
        y = m(2)%p + r2(iy)
        yp = y  - planet%position(2)
        ys = yp - ds2(2)
        do ix=m(1)%lb,m(1)%ub
          x = m(1)%p + r1(ix)
          xp = x  - planet%position(1)
          xs = xp - ds2(1)
          !---------------------------------------------------------------
          ! The force from the planet is included explicitly, and is oriented
          ! away from the planet, which is at the center of the local coordinates.
          ! The L1 and L2 Lagrange points are the result of net centrifugal
          ! forces that point in the opposite direction to the force of gravity,
          ! cancelling at a distance where $\Omega^2 r$ is equal to $GM/r^2$
          !---------------------------------------------------------------
          rp = sqrt(xs**2+yp**2+zp**2)
          ff(ix,iy,iz,1) = planet%newton(rp,xs)
          rp = sqrt(xp**2+ys**2+zp**2)
          ff(ix,iy,iz,2) = planet%newton(rp,ys)
          rp = sqrt(xp**2+yp**2+zs**2)
          ff(ix,iy,iz,3) = planet%newton(rp,zs)
          !---------------------------------------------------------------
          ! Add the vertical and radial rotating coordinate frame accelerations.
          ! The net force in the radial direction increases with increasing x,
          ! giving that x now points in the radial direction.  It has two
          ! contributions:  The centripetal force \Omega^2 x, and the reduction
          ! of the force of gravity from the central star, which contributes a
          ! factor of 2 more, for a total factor of 3.  Vertically, there
          ! is a net force that comes from the projection of the force of
          ! gravity onto the z-axis, and has the same scaling constant as the
          ! radial force, but the opposite sign.
          !---------------------------------------------------------------
          ff(ix,iy,iz,1) = ff(ix,iy,iz,1) + 3.0*planet%omega(3)**2*xs
          ff(ix,iy,iz,3) = ff(ix,iy,iz,3) -     planet%omega(3)**2*zs
        end do
      end do
    end do
  end if

  !---------------------------------------------------------------------
  ! Friction term, to damp out initial sloshing motions, relative to the
  ! average shear.
  !---------------------------------------------------------------------
  if (planet%friction > 0.0) then
    if (self%solver(1:6)=='ramses') then
      do iz=m(3)%lb,m(3)%ub
        z = m(3)%p + r3(iz)
        zp = z  - planet%position(3)
        do iy=m(2)%lb,m(2)%ub
          y = m(2)%p + r2(iy)
          yp = y  - planet%position(2)
          do ix=m(1)%lb,m(1)%ub
            x = m(1)%p + r1(ix)
            xp = x  - planet%position(1)
            uxx = -planet%v_accrete
            uyy = -(1.5+planet%trap)*planet%omega(3)*xp - planet%v_head
            ff(ix,iy,iz,1) = ff(ix,iy,iz,1) - planet%friction*(p(ix,iy,iz,1)/d(ix,iy,iz)-uxx)
            ff(ix,iy,iz,2) = ff(ix,iy,iz,2) - planet%friction*(p(ix,iy,iz,2)/d(ix,iy,iz)-uyy)
            ff(ix,iy,iz,3) = ff(ix,iy,iz,3) - planet%friction*(p(ix,iy,iz,3)/d(ix,iy,iz))
          end do
        end do
      end do
    else if (self%solver(1:6)=='stagge'.or.self%solver(1:6)=='zeus3d') then
      do iz=m(3)%lb,m(3)%ub
        z = m(3)%p + r3(iz)
        zp = z  - planet%position(3)
        zs = zp - ds2(3)
        do iy=m(2)%lb,m(2)%ub
          y = m(2)%p + r2(iy)
          yp = y  - planet%position(2)
          ys = yp - ds2(2)
          do ix=m(1)%lb,m(1)%ub
            x = m(1)%p + r1(ix)
            xp = x  - planet%position(1)
            xs = xp - ds2(1)
            uxx = -planet%v_accrete
            uyy = -(1.5+planet%trap)*planet%omega(3)*xs - planet%v_head
            ff(ix,iy,iz,1) = ff(ix,iy,iz,1) - planet%friction*(Ux(ix,iy,iz)-uxx)
            ff(ix,iy,iz,2) = ff(ix,iy,iz,2) - planet%friction*(Uy(ix,iy,iz)-uyy)
            ff(ix,iy,iz,3) = ff(ix,iy,iz,3) - planet%friction*(Uz(ix,iy,iz))
          end do
        end do
      end do
    else
      print *, trim(self%solver)
      error stop "force_mod(1) :: unknown solver"
    end if
  end if
  !---------------------------------------------------------------------
  ! Add the classical Coriolis force term: -2*(\Omega x U). A positive omega(3)
  ! corresponds to rotation in the y-direction, if x points radially outwards.
  ! Hence, a positive velocity in the y-direction increases the centripetal force,
  ! and leads to acceleration in the positive x-direction -- vice versa for
  ! motion in the x-direction; a motion inwards (-x) leads to acceleration in
  ! the forward (+y) direction
  !---------------------------------------------------------------------
  if (self%solver(1:6)=='ramses') then
    ff(:,:,:,1) = ff(:,:,:,1) - 2.*(planet%omega(2)*p(:,:,:,3)-planet%omega(3)*p(:,:,:,2))/d
    ff(:,:,:,2) = ff(:,:,:,2) - 2.*(planet%omega(3)*p(:,:,:,1)-planet%omega(1)*p(:,:,:,3))/d
    ff(:,:,:,3) = ff(:,:,:,3) - 2.*(planet%omega(1)*p(:,:,:,2)-planet%omega(2)*p(:,:,:,1))/d
  else if (self%solver(1:6)=='stagge'.or.self%solver(1:6)=='zeus3d') then
    ff(:,:,:,1) = ff(:,:,:,1) - 2.*(planet%omega(2)*st_dn(st_up(Uz,3),1) - planet%omega(3)*st_dn(st_up(Uy,2),1))
    ff(:,:,:,2) = ff(:,:,:,2) - 2.*(planet%omega(3)*st_dn(st_up(Ux,1),2) - planet%omega(1)*st_dn(st_up(Uz,3),2))
    ff(:,:,:,3) = ff(:,:,:,3) - 2.*(planet%omega(1)*st_dn(st_up(Uy,2),3) - planet%omega(2)*st_dn(st_up(Ux,1),3))
  else
    print *, trim(self%solver)
    error stop "force_mod(2) :: unknown solver"
  end if
  !---------------------------------------------------------------------
  ! The force of gravity from the planet, -G m/r^2, balances the
  ! tripled centripetal force 3\Omega*r at r**3 = (G m)/(3 \Omega^2)
  ! and since \Omega^2 = G M/R^3, r = R (m/(3M))^(1/3)
  !---------------------------------------------------------------------
  if (self%id==io%id_debug) then
    if (maxval(abs(m(2:3)%p)/m(2:3)%d) < 0.5) then
      iz = (m(3)%li+m(3)%ui)/2
      iy = (m(2)%li+m(2)%ui)/2
      do ix=m(1)%li,m(1)%ui
        xs = m(1)%p + r1(ix) - planet%position(1) - ds2(1)
        if (self%solver(1:6)=='stagge'.or.self%solver(1:6)=='zeus3d') then
          print '("FORCE:",i6,i4,f12.2,1p,4e14.3)', self%id, ix, m(1)%p + r1(ix), &
            ff(ix,iy,iz,1), Uy(ix,iy,iz), &
            2.*planet%omega(3)*Uy(ix,iy,iz), 3.*planet%omega(3)**2*xs
        else
          print '("FORCE:",i6,i4,f12.2,1p,4e14.3)', self%id, ix, m(1)%p + r1(ix), &
            ff(ix,iy,iz,1), p(ix,iy,iz,2)/d(ix,iy,iz), &
            2.*planet%omega(3)*p(ix,iy,iz,2)/d(ix,iy,iz), 3.*planet%omega(3)**2*xs
        end if
      end do
    end if
  end if
  end associate
  call trace_end
END FUNCTION disk_perhaps_with_planet

!===============================================================================
!> Forcing on a single 3-D wavenumber, changing amplitude and phase periodically
!===============================================================================
FUNCTION single_solenoidal (self, time, Ux, Uy, Uz, m) RESULT (ff)
  class(force_t)        :: self
  real(8)               :: time
  class(mesh_t), pointer, dimension(:) :: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3) :: ff
  real, dimension(:,:,:):: Ux, Uy, Uz
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
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, solenoidal_params, iostat=iostat)
    if (io%master) write (*, solenoidal_params)
  end if
  call io%lock%unset
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

FUNCTION st_up (a, i) result (b)
  real, dimension(:,:,:) :: a
  real, dimension(size(a,1), size(a,2), size(a,3)):: b
  real    :: c
  integer :: ix, iy, iz, i, n(3)
!...............................................................................
  c = 0.5
  n = shape(a)
  select case(i)
  case(1)
    b(n(1),:,:) = a(n(1),:,:)
    do iz=1,n(3)
      do iy=1,n(2)
        do ix=1,n(1)-1
          b(ix,iy,iz) = c*(a(ix+1,iy,iz) + a(ix,iy,iz))
        end do
      end do
    end do
  case(2)
    b(:,n(2),:) = a(:,n(2),:)
    do iz=1,n(3)
      do iy=1,n(2)-1
        do ix=1,n(1)
          b(ix,iy,iz) = c*(a(ix,iy+1,iz) + a(ix,iy,iz))
        end do
      end do
    end do
  case(3)
    b(:,:,n(3)) = a(:,:,n(3))
    do iz=1,n(3)-1
      do iy=1,n(2)
        do ix=1,n(1)
          b(ix,iy,iz) = c*(a(ix,iy,iz+1) + a(ix,iy,iz))
        end do
      end do
    end do
  end select
END FUNCTION st_up

FUNCTION st_dn (a, i) result (b)
  real, dimension(:,:,:) :: a
  real, dimension(size(a,1), size(a,2), size(a,3)):: b
  real    :: c
  integer :: ix, iy, iz, i, n(3)
!...............................................................................
  c = 0.5
  n = shape(a)
  select case(i)
  case(1)
    b(1,:,:) = a(1,:,:)
    do iz=1,n(3)
      do iy=1,n(2)
        do ix=2,n(1)
          b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix-1,iy,iz))
        end do
      end do
    end do
  case(2)
    b(:,1,:) = a(:,1,:)
    do iz=1,n(3)
      do iy=2,n(2)
        do ix=1,n(1)
          b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy-1,iz))
        end do
      end do
    end do
  case(3)
    b(:,:,1) = a(:,:,1)
    do iz=2,n(3)
      do iy=1,n(2)
        do ix=1,n(1)
          b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy,iz-1))
        end do
      end do
    end do
  end select
END FUNCTION st_dn

END MODULE
