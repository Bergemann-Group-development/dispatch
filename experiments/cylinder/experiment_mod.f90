!===============================================================================
!> Experiment initialization and update
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE cpatch_mod
  USE bits_mod
  USE boundaries_mod
  USE mesh_mod
  USE scaling_mod
  USE download_mod
  USE scratch_mod
  implicit none
  private
  type, public, extends(cpatch_t):: experiment_t
    real:: mass=1.0
    real:: radius=1.0
  contains
    procedure:: init
    procedure:: update
    procedure:: dnload
  end type
  real    :: c_sound     =  0.01
  real    :: v_kepler    =  1.
  real    :: d_mid       =  1.0
  real    :: p_mid       = -1.0
  logical :: with_planet = .false.
CONTAINS

!===============================================================================
!> Read parameters and set up initial state
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !.............................................................................
  logical, save:: first_time=.true.
  namelist /experiment_params/ c_sound, v_kepler, d_mid, p_mid, with_planet
  !.............................................................................
  integer:: l(3), u(3)
  real(8):: v_k
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init')
  call self%cpatch_t%init
  !----------------------------------------------------------------------------
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io_unit%input)
    read (io_unit%input, experiment_params)
    write (io_unit%output, experiment_params)
  end if
  call io%lock%unset
  !----------------------------------------------------------------------------
  v_k = v_kepler/sqrt(self%r)                   ! mid point Kepler speed
  self%omega = [0d0, 0d0, v_k/self%r]           ! angular velocity
  self%gamma = 1.0                              ! isothermal
  self%csound = c_sound                         ! solver sound speed
  !----------------------------------------------------------------------------
  l = [self%mesh(3)%lb, self%mesh(2)%lb, self%mesh(1)%lb]
  u = [self%mesh(3)%ub, self%mesh(2)%ub, self%mesh(1)%ub]
  call initial_condition (self, l, u)
  l = [self%mesh(3)%li, self%mesh(2)%li, self%mesh(1)%li]
  u = [self%mesh(3)%ui, self%mesh(2)%ui, self%mesh(1)%ui]
  d  => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%d ,self%it,1)
  e  => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e ,self%it,1)
  px => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%px,self%it,1)
  py => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%py,self%it,1)
  pz => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%pz,self%it,1)
  print '(i6,1p,4g12.3)', self%id, minval(d), maxval(d), maxval(abs(px/d)), maxval(abs(py/d))
  call trace_end
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call boundary (self)
  call force (self)
  call self%cpatch_t%update
  call trace_end
END SUBROUTINE update

!===============================================================================
!===============================================================================
SUBROUTINE dnload (self, only)
  class(experiment_t):: self
  integer, optional:: only
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%dnload')
  if (self%use_data_hub) then
    call self%data_hub%update (self%link, only=only)
  else
    call download%download_link (self%link, only=only)
  end if
  call trace%end()
END SUBROUTINE dnload

!===============================================================================
!> Power law isothermal disk conditions
!===============================================================================
SUBROUTINE initial_condition (self, l, u)
  class(experiment_t):: self
  integer:: l(3), u(3)
  !.............................................................................
  integer:: ix, iy, iz
  real(8):: x, xs, y, ys, z, r, phi, v_k, h, ds2
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  !----------------------------------------------------------------------------
  ds2 = 0d0
  if (self%kind(1:7)=='stagger' .or. &
      self%kind(1:6)=='zeus3d') then
    ds2 = -self%ds(1)/2d0
  end if
  d  => self%mem(:,:,:,self%idx%d ,self%it,1)
  e  => self%mem(:,:,:,self%idx%e ,self%it,1)
  px => self%mem(:,:,:,self%idx%px,self%it,1)
  py => self%mem(:,:,:,self%idx%py,self%it,1)
  pz => self%mem(:,:,:,self%idx%pz,self%it,1)
  do iz=l(3),u(3)
  do iy=l(2),u(2)
  do ix=l(1),u(1)
    x = self%r + self%mesh(1)%r(ix)
    y =        + self%mesh(2)%r(iy)
    z = self%z + self%mesh(3)%r(iz)
    xs = x + ds2                                        ! staggered
    ys = y + ds2                                        ! staggered
    !---------------------------------------------------------------------------
    phi = atan2(xs,y)                                   ! local angle
    r = sqrt(xs**2 + y**2)                              ! local radius
    v_k = v_kepler/sqrt(r)                              ! Kepler speed
    h = r*c_sound/v_k                                   ! scale height
    v_k = v_k - self%r*self%omega(3)                    ! net Kepler speed
    d (ix,iy,iz) = d_mid*r**p_mid*exp(-(z/h)**2)        ! Gaussian profile
    px(ix,iy,iz) = - v_k*sin(phi)*d(ix,iy,iz)           ! radial momentum
    !---------------------------------------------------------------------------
    phi = atan2(x,ys)                                   ! local angle
    r = sqrt(x**2 + ys**2)                              ! local radius
    v_k = v_kepler/sqrt(r)                              ! Kepler speed
    h = r*c_sound/v_k                                   ! scale height
    v_k = v_k - self%r*self%omega(3)                    ! net Kepler speed
    d (ix,iy,iz) = d_mid*r**p_mid*exp(-(z/h)**2)        ! Gaussian profile
    py(ix,iy,iz) = + v_k*cos(phi)*d(ix,iy,iz)           ! azimuthal momentum
    pz(ix,iy,iz) = 0.0                                  ! vertical momentum
    !---------------------------------------------------------------------------
    phi  = atan2(x ,y)                                  ! local angle
    r  = sqrt(x**2 + y**2)                              ! local radius
    v_k = v_kepler/sqrt(r)                              ! Kepler speed
    h = r*c_sound/v_k                                   ! scale height
    v_k = v_k - self%r*self%omega(3)                    ! net Kepler speed
    d (ix,iy,iz) = d_mid*r**p_mid*exp(-(z/h)**2)        ! Gaussian profile
    e (ix,iy,iz) = d(ix,iy,iz)*c_sound**2               ! isothermal
  end do
  end do
  end do
  if (self%kind(1:7)=='ramses_' .and. &
      self%gamma > 1.0) then
    e = d*c_sound**2/(self%gamma*(self%gamma-1.0)) &
      + 0.5*(px**2+py**2+pz**2)/d
  end if
END SUBROUTINE initial_condition

!===============================================================================
!===============================================================================
SUBROUTINE boundary (self)
  class(experiment_t):: self
  integer:: l(3), u(3)
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%boundary')
  if (self%boundaries%is_set (bits%xl)) then
    l = [self%mesh(3)%lb, self%mesh(2)%lb, self%mesh(1)%lb]
    u = [self%mesh(3)%lo, self%mesh(2)%lb, self%mesh(1)%lb]
    call initial_condition (self, l, u)
  end if
  if (self%boundaries%is_set (bits%xu)) then
    l = [self%mesh(3)%uo, self%mesh(2)%lb, self%mesh(1)%lb]
    u = [self%mesh(3)%ub, self%mesh(2)%lb, self%mesh(1)%lb]
    call initial_condition (self, l, u)
  end if
  if (self%boundaries%is_set (bits%zl)) then
    l = [self%mesh(3)%lb, self%mesh(2)%lb, self%mesh(1)%lb]
    u = [self%mesh(3)%ub, self%mesh(2)%lb, self%mesh(1)%lo]
    call initial_condition (self, l, u)
  end if
  if (self%boundaries%is_set (bits%zu)) then
    l = [self%mesh(3)%lb, self%mesh(2)%lb, self%mesh(1)%uo]
    u = [self%mesh(3)%ub, self%mesh(2)%lb, self%mesh(1)%ub]
    call initial_condition (self, l, u)
  end if
  call trace_end
END SUBROUTINE boundary

!===============================================================================
!> Classical Coriolis force and similar force terms
!===============================================================================
SUBROUTINE force (self)
  class(experiment_t):: self
  !.............................................................................
  real(8):: time
  real, dimension(:,:,:), pointer:: Ux, Uy, Uz, d
  real, dimension(:,:,:,:), pointer:: p, ff, UU
  integer:: ix, iy, iz, iv
  real(8):: rp, newton
  real(8):: x, y, z, xs, ys, zs, xp, yp, zp, ds2(3), uxx, uyy
  class(mesh_t), pointer, dimension(:):: m
  integer:: l(3), u(3)
  !-----------------------------------------------------------------------------
  call trace_begin('experiment_t%force')
  m => self%mesh
  if (self%kind(1:7)=='ramses_') then
    ds2 = 0.
  else if (self%kind(1:7)=='stagger') then
    ds2 = 0.5*m(1)%d
  else if (self%kind(1:6)=='zeus3d') then
    ds2 = 0.5*m(1)%d
  else
    error stop 'force'
  end if
  allocate (scratch%force_per_unit_mass(m(1)%gn,m(2)%gn,m(3)%gn,3))
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r, ff => scratch%force_per_unit_mass)
  if (with_planet) then
    newton = self%mass*scaling%grav
  else
    newton = 0.0
  end if
  ff = 0.0
  do iz=m(3)%lb,m(3)%ub
    z = r3(iz)
    zs = z - ds2(3)
    do iy=m(2)%lb,m(2)%ub
      y = r2(iy)
      ys = y - ds2(2)
      do ix=m(1)%lb,m(1)%ub
        x = r1(ix)
        xs = x - ds2(1)
        !---------------------------------------------------------------
        ! The force from the planet is included explicitly, and is oriented
        ! away from the planet, which is at the center of the local coordinates.
        ! The L1 and L2 Lagrange points are the result of net centrifugal
        ! forces that point in the opposite direction to the force of gravity,
        ! cancelling at a distance where $\Omega^2 r$ is equal to $GM/r^2$
        !---------------------------------------------------------------
        !rp = sqrt(xs**2+yp**2+zp**2)
        !rp = max(rp,self%radius*0.1)
        !ff(ix,iy,iz,1) = -newton*xs/rp**3
        !rp = sqrt(xp**2+ys**2+zp**2)
        !rp = max(rp,self%radius*0.1)
        !ff(ix,iy,iz,2) = -newton*ys/rp**3
        !rp = sqrt(xp**2+yp**2+zs**2)
        !rp = max(rp,self%radius*0.1)
        !ff(ix,iy,iz,3) = -newton*zs/rp**3
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
        ff(ix,iy,iz,1) = ff(ix,iy,iz,1) + 3.0*self%omega(3)**2*xs
        ff(ix,iy,iz,3) = ff(ix,iy,iz,3) -     self%omega(3)**2*zs
      end do
    end do
  end do
  !---------------------------------------------------------------------
  ! Add the classical Coriolis force term: -2*(\Omega x U). A positive omega(3)
  ! corresponds to rotation in the y-direction, if x points radially outwards.
  ! Hence, a positive velocity in the y-direction increases the centripetal force,
  ! and leads to acceleration in the positive x-direction -- vice versa for
  ! motion in the x-direction; a motion inwards (-x) leads to acceleration in
  ! the forward (+y) direction
  !---------------------------------------------------------------------
  if (self%kind(1:6)=='ramses') then
    d => self%mem(:,:,:,self%idx%d,self%it,1)
    p => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
    ff(:,:,:,1) = ff(:,:,:,1) - 2.*(self%omega(2)*p(:,:,:,3)-self%omega(3)*p(:,:,:,2))/d
    ff(:,:,:,2) = ff(:,:,:,2) - 2.*(self%omega(3)*p(:,:,:,1)-self%omega(1)*p(:,:,:,3))/d
    ff(:,:,:,3) = ff(:,:,:,3) - 2.*(self%omega(1)*p(:,:,:,2)-self%omega(2)*p(:,:,:,1))/d
  else if (self%kind(1:6)=='stagge'.or.self%kind(1:6)=='zeus3d') then
    allocate (UU(self%gn(1),self%gn(2),self%gn(3),3))
    Ux => UU(:,:,:,1)
    Uy => UU(:,:,:,2)
    Uz => UU(:,:,:,3)
    call self%p2U (UU, self%it)
    ff(:,:,:,1) = ff(:,:,:,1) - 2.*(self%omega(2)*st_dn(st_up(Uz,3),1) - self%omega(3)*st_dn(st_up(Uy,2),1))
    ff(:,:,:,2) = ff(:,:,:,2) - 2.*(self%omega(3)*st_dn(st_up(Ux,1),2) - self%omega(1)*st_dn(st_up(Uz,3),2))
    ff(:,:,:,3) = ff(:,:,:,3) - 2.*(self%omega(1)*st_dn(st_up(Uy,2),3) - self%omega(2)*st_dn(st_up(Ux,1),3))
    deallocate (UU)
  else
    print *, trim(self%kind)
    error stop "force_mod(2) :: unknown solver"
  end if
  l = [self%mesh(3)%li, self%mesh(2)%li, self%mesh(1)%li]
  u = [self%mesh(3)%ui, self%mesh(2)%ui, self%mesh(1)%ui]
  !print '(i6,1p,4g13.3)', self%id, &
    !minval(abs(ff(l(1):u(1),l(2):u(2),l(3):u(3),1))), &
    !maxval(abs(ff(l(1):u(1),l(2):u(2),l(3):u(3),1))), &
    !minval(abs(ff(l(1):u(1),l(2):u(2),l(3):u(3),2))), &
    !maxval(abs(ff(l(1):u(1),l(2):u(2),l(3):u(3),2)))
  end associate
  call trace_end
END SUBROUTINE force

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

END MODULE experiment_mod
