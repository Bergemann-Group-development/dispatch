!===============================================================================
!> Module to set the initial conditions on each patch.
!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE mesh_mod
  USE index_mod
  USE math_mod
  USE scaling_mod
  USE units_mod
  implicit none
  private

  real(8), save:: aspect_ratio = 0.05     ! aspect ratio of the disk (H/R), where H is the scale height
  real(8), save:: opening_angle = 2.0     ! given as a multiple of `aspectratio`
  real(8), save:: disc_r0 = 0.02          ! [AU] inner disc radius/scale radius
  real(8), save:: disc_rmax = 10.0        ! [AU] size of the disc
  real(8), save:: disc_mass = 0.1         ! [M_Sun] total disc mass
  real(8), save:: star_mass = 1.0         ! [M_Sun] mass of central star
  real(8), save:: rpindex = 1.5           ! (negative) radial power index of the disc density
  real(8), save:: pbeta0 = 0.922e3        ! initial plasma-beta of the field at `r0_disc`
  real(8), save:: mbend = 0.5             ! initial magnetic field bending parameter
  real(8), save:: gamma                   ! ratio of specific heats; should be overridden; this is here
                                          ! for convenience!!
  real(8), save:: density_floor = 1.0e-24 ! density floor value
  logical, save:: rotating_disk = .true.  ! switch for turning disk rotation on/off

  type, public:: initial_t
    logical:: mhd
    real(8):: aspect_ratio
    real(8):: opening_angle
    real(8):: disc_r0
    real(8):: disc_rmax
    real(8):: disc_mass
    real(8):: rpindex
    real(8):: pbeta0
    real(8):: mbend
    real(8):: gamma
    real(8):: star_mass
    real(8):: density_floor
    logical:: with_density_floor
    logical:: rotating_disk
    character(len=64):: solver
  contains
    procedure:: init
    procedure:: condition
  end type

CONTAINS
!===============================================================================
!> Setup the data structure, read in any namelist variables.
!===============================================================================
SUBROUTINE init (self, solver, gamma)
  class(initial_t):: self
  character(len=64):: solver
  real, optional:: gamma
  logical, save:: first_time=.true.
  logical, save:: with_density_floor=.true.
  integer:: iostat
  namelist /IC_params/ aspect_ratio, opening_angle, disc_r0, rpindex, disc_mass, &
                            pbeta0, mbend, disc_rmax, star_mass, with_density_floor, &
                            density_floor, rotating_disk
  !----------------------------------------------------------------------------
  call trace%begin('initial_t%init')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, IC_params, iostat=iostat)
    if (iostat/=0) call io%namelist_warning ('IC_params')
    if (io%master) write (*, IC_params)
  end if
  call io%lock%unset

  self%aspect_ratio = aspect_ratio
  self%opening_angle = opening_angle
  self%disc_r0 = disc_r0 * cgs%au ! convert to cgs
  self%disc_rmax = disc_rmax * cgs%au ! convert to cgs
  self%disc_mass = disc_mass * cgs%m_sun
  self%rpindex = rpindex
  self%pbeta0 = pbeta0
  self%mbend = mbend
  self%gamma = gamma
  self%star_mass = star_mass * cgs%m_sun
  self%density_floor = density_floor
  self%with_density_floor = with_density_floor
  self%rotating_disk = rotating_disk
  self%solver = trim(solver)

  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Assign the initial conditions.
!> The static gravitational potential is set in `potential_mod`.
!===============================================================================
SUBROUTINE condition (self, m, ff, idx)
  class(initial_t):: self
  type(index_t):: idx
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: ff
  class(mesh_t), pointer, dimension(:):: m

  class(mesh_t), pointer:: m1, m2, m3
  integer:: i1, i2, i3, ione, jone
  real(8):: gm_s, vkep, rho_mid, disk_angle, rho0, pgas0, q0, r_sph, r_cyl, bpol0, &
            cot_theta, z, theta, x1, x2, switch_rotating, vkeptilde2
  real(8), dimension(:,:,:), allocatable:: psi ! magnetic flux function
  logical, save:: first_time = .true.
  !----------------------------------------------------------------------------
  call trace%begin ('initial_t%condition')

  m1 => m(1)
  m2 => m(2)
  m3 => m(3)

  ione = merge(1, 0, m1%n > 1)
  jone = merge(1, 0, m2%n > 1)

  gm_s = scaling%grav * self%star_mass
  disk_angle = atan(self%opening_angle * self%aspect_ratio)
  ! using the total disc mass, calculate the midplane density at r0;
  ! note that the "error function" is not part of the Fortran standard.
  q0 = 3.0 - self%rpindex
  rho0 = (2.0 * math%pi)**1.5 * self%disc_r0**self%rpindex / q0 &
       * (self%disc_rmax**q0 - self%disc_r0**q0) * self%aspect_ratio &
       * erf(1.0/sqrt(2.0)/self%aspect_ratio)
  rho0 = self%disc_mass / rho0
  switch_rotating = 1.0
  if (.not. self%rotating_disk) switch_rotating = 0.0

  pgas0 = rho0 * (gm_s / self%disc_r0) / (1.0 + self%rpindex)
  ! initial magnetic field strength is set primarily by `pbeta0`.
  if (self%pbeta0 == 0) then
    bpol0 = 0.0
  else
    bpol0 = sqrt(2.0 * pgas0 / self%pbeta0)
  end if

  if (first_time .and. io%master) then
    print*,"rho_0 =",rho0
    print*,"Corona/disk density ratio:",self%aspect_ratio**2 * (1.0 + self%rpindex)
    print*,"1 Keplerian orbit at 'disc_r0' (in s):",2.0 * math%pi / sqrt(scaling%grav * self%star_mass / self%disc_r0**3)
    print*,"1 Keplerian orbit at 'disc_rmax' (in s):",2.0 * math%pi / sqrt(scaling%grav * self%star_mass / self%disc_rmax**3)
    print*,"Initial field strength at the midplane and r = r_0:",bpol0
    first_time = .false.
  end if

  ! initialise magnetic flux function
  allocate (psi(m1%gn,m2%gn,m3%gn))
  if (self%pbeta0 == 0.0) then
    psi = 0.0
  else
  ! NOTE!: the flux function should be *edge-centred*!!!
    do i3=m3%lb,m3%gn
      do i2=m2%lb,m2%gn
        ! face-centred
        if (i2 == m2%gn) then
          x2 = m2%centre_nat + m2%r(m2%ub) + 0.5*m2%d
        else
          x2 = m2%centre_nat + m2%r(i2) - 0.5*m2%d
        end if
        do i1=m1%lb,m1%gn
          ! face-centred
          if (i1 == m1%gn) then
            x1 = m1%centre_nat + m1%r(m1%ub) + 0.5*m1%d
          else
            x1 = m1%centre_nat + m1%r(i1) - 0.5*m1%d
          end if
          if (m1%type == 'spherical_mesh') then
            theta = x2
            r_sph = abs(x1)
            z = abs(r_sph * cos(theta))
            r_cyl = abs(r_sph * sin(theta))
          else if (m1%type == 'cylindrical_mesh') then
            z = x1
            r_cyl = x2
            r_sph = sqrt(z**2 + r_cyl**2)
            theta = acos(z / r_cyl)
          end if
          cot_theta = 1.0 / (tan(theta) + tiny(0.0))
          ! the introduction of the `sign` function is to ensure B_theta changes sign
          ! across a symmetry axis.
          psi(i1,i2,i3) = 4.0 * self%disc_r0**2 * bpol0 * (r_cyl / self%disc_r0)**0.75 &
                        * self%mbend**1.25 / 3.0 / (self%mbend**2 + cot_theta**2)**0.625 &
                        * sign(1.0d0,sin(theta))
        end do
      end do
    end do
  end if

  ! intialise field variables
  do i3=m3%lb,m3%ub
    do i2=m2%lb,m2%ub
      x2 = m2%centre_nat + m2%r(i2) ! cell-centred
      do i1=m1%lb,m1%ub
        x1 = m1%centre_nat + m1%r(i1) ! cell-centred
        if (m1%type == 'spherical_mesh') then
          theta = x2
          r_sph = abs(x1)
          z = abs(r_sph * cos(theta))
          r_cyl = abs(r_sph * sin(theta))
        else if (m1%type == 'cylindrical_mesh') then
          z = x1
          r_cyl = x2
          r_sph = sqrt(z**2 + r_cyl**2)
          theta = acos(z / r_cyl)
        end if

        vkep = sqrt(gm_s / r_cyl)
        vkeptilde2 = gm_s / r_sph
        rho_mid = rho0 * (r_sph / self%disc_r0)**(-self%rpindex) 
        if (theta < (0.5*math%pi - disk_angle) .or. theta > (0.5*math%pi + disk_angle)) then
          ! in the corona
          ff(i1,i2,i3,idx%d ) = rho_mid * self%aspect_ratio**2 * (1.0 + self%rpindex) &
                             * exp(-0.5 * cos(0.5*math%pi - disk_angle)**2 / self%aspect_ratio**2)
          pgas0 = ff(i1,i2,i3,idx%d) * vkeptilde2 / (1.0 + self%rpindex)
          ff(i1,i2,i3,idx%e ) = pgas0 / (self%gamma - 1.0)
          ff(i1,i2,i3,idx%pz) = 0.0
        else
          ! in the disk
          ff(i1,i2,i3,idx%d ) = rho_mid * exp(-0.5 * cos(theta)**2 / self%aspect_ratio**2)
          pgas0 = ff(i1,i2,i3,idx%d) * vkeptilde2 * self%aspect_ratio**2
          ff(i1,i2,i3,idx%e ) = pgas0 / (self%gamma - 1.0)
          ff(i1,i2,i3,idx%pz) = switch_rotating * vkep * sqrt(1.0 - (1.0 + self%rpindex) * self%aspect_ratio**2) &
                              * ff(i1,i2,i3,idx%d ) * m1%h31c(i1) * m2%h32c(i2)
        end if

        ! apply density floor (if desired)
        if (self%with_density_floor) then
          ff(i1,i2,i3,idx%d) = max(ff(i1,i2,i3,idx%d), density_floor)
        end if

        ff(i1,i2,i3,idx%px) = 0.0
        ff(i1,i2,i3,idx%py) = 0.0
        if (self%mhd) then
          ! initially purely poloidal field
          ! the form below should be safe for both spherical and cylindrical coords.
          ff(i1,i2,i3,idx%bx) = (psi(i1,i2+jone,i3) - psi(i1,i2,i3)) &
                              / (m2%d * m1%h2f(i1) * m1%h31f(i1) * m2%h32c(i2))
          ff(i1,i2,i3,idx%by) = -(psi(i1+ione,i2,i3) - psi(i1,i2,i3)) &
                              / (m1%d * m1%h31c(i1) * m2%h32f(i2) + tiny(0.0))
          ff(i1,i2,i3,idx%bz) = 0.0
        end if
      end do
      ! recall that the ZEUS solver has one extra zone (ub+1) for staggered quantities.
      ! it is particularly important to account for this when the boundary values of the
      ! magnetic fields are not being updated directly (but rather using EMFs).
      if (self%solver(1:4) == 'zeus' .and. m1%n > 1) then
        ff(m1%gn,i2,i3,idx%bx) = (psi(m1%gn,i2+jone,i3) - psi(m1%gn,i2,i3)) &
                              / (m2%d * m1%h2f(m1%gn) * m1%h31f(m1%gn) * m2%h32c(i2))
      end if
    end do
    if (self%solver(1:4) == 'zeus' .and. m2%n > 1) then
      do i1=m1%lb,m1%ub
        ff(i1,m2%gn,i3,idx%by) = -(psi(i1+ione,m2%gn,i3) - psi(i1,m2%gn,i3)) &
                               / (m1%d * m1%h31c(i1) * m2%h32f(m2%gn))
      end do
    end if
  end do
  if (self%solver(1:4) == 'zeus' .and. m3%n > 1) then
    ff(m1%lb:m1%ub,m2%lb:m2%ub,m3%gn,idx%bz) = 0.0
  end if

  nullify (m1, m2, m3)
  deallocate (psi)
  call trace%end()
END SUBROUTINE condition

END MODULE initial_mod
