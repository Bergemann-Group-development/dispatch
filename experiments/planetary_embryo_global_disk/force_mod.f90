!> =============================================================================
!> Direct force of gravity
!> =============================================================================
MODULE force_mod
  USE io_mod
  USE trace_mod
  USE scaling_mod
  USE mesh_mod
  USE mhd_mod
  USE parameters_mod
  implicit none
  private

  type, public:: force_t
    real(8), pointer :: omega
    procedure(star_only), pointer:: selected=>null()
  contains
    procedure init
  end type
  type(force_t), public:: force

CONTAINS
!> =============================================================================
!> Select type of force -- could be more than one planet, for example
!> =============================================================================
SUBROUTINE init (self, patch)
  class(force_t) :: self
  class(mhd_t), target :: patch
  !-----------------------------------------------------------------------------
  call trace_begin('force_t%init')

  ! select implementation for the force source term
  self%selected => star_and_planet

  self%omega => patch%grid_vel(3)

  call trace_end
END SUBROUTINE

!> =============================================================================
!> Returns the combined gravitational potential of a star and planet in
!  cylindrical coords.
!> =============================================================================
FUNCTION star_only (self, m, time, dt) RESULT (gp)
  class(force_t):: self
  class(mesh_t), dimension(:), allocatable, target :: m
  class(mesh_t), pointer :: m1, m2, m3
  real(8), intent(in) :: time, dt
  real(8), dimension(m(1)%gn,m(2)%gn,m(3)%gn):: gp
  !.............................................................................
  integer:: i1, i2, i3
  real(8):: z, r, phi, z_s, r_s, gm_s, rvec_s
  !-----------------------------------------------------------------------------
  call trace_begin('force_mod::star_only')

  m1 => m(1)
  m2 => m(2)
  m3 => m(3)

  gm_s = scaling%grav * experiment_params%mstar
  do i3=m3%lb,m3%ub
    phi = m3%centre_nat + m3%r(i3)
    do i2=m2%lb,m2%ub
      r = m2%centre_nat + m2%r(i2)
      r_s = r - 0.0
      do i1=m1%lb,m1%ub
        z = m1%centre_nat + m1%r(i1)
        z_s = z - 0.0

        rvec_s = sqrt(z_s**2 + r_s**2)
        gp(i1,i2,i3) = -gm_s / rvec_s
      end do
    end do
  end do

  nullify (m1, m2, m3)
  call trace_end
END FUNCTION star_only

!> =============================================================================
!> Returns the combined gravitational potential of a star and planet in
!  cylindrical coords.
!> =============================================================================
FUNCTION star_and_planet (self, m, time, dt) RESULT (gp)
  class(force_t):: self
  class(mesh_t), dimension(:), allocatable, target :: m
  class(mesh_t), pointer :: m1, m2, m3
  real(8), intent(in) :: time, dt
  real(8), dimension(3) :: planet_pos
  real(8), dimension(m(1)%gn,m(2)%gn,m(3)%gn):: gp, gp2, rvec_p, gp3, gp4
  !.............................................................................
  integer:: i1, i2, i3
  real(8):: z, r, phi, z_s, r_s, gm_s, gm_p, rvec_s
  !-----------------------------------------------------------------------------
  call trace_begin('force_mod::star_and_planet')

  m1 => m(1)
  m2 => m(2)
  m3 => m(3)

  gm_s = scaling%grav * experiment_params%mstar
  gm_p = scaling%grav * experiment_params%planet%mass
  planet_pos = experiment_params%planet%cyl_pos
  do i3=m3%lb,m3%ub
    phi = m3%centre_nat + m3%r(i3)
    do i2=m2%lb,m2%ub
      r = m2%centre_nat + m2%r(i2)
      r_s = r - 0.0
      do i1=m1%lb,m1%ub
        z = m1%centre_nat + m1%r(i1)
        z_s = z - 0.0

        rvec_s = sqrt(z_s**2 + r_s**2)
        rvec_p(i1,i2,i3) = planet_pos(2)**2 + r**2 &
                 - 2.0d0 * planet_pos(2) * r * cos(abs(phi - planet_pos(3)))
        gp(i1,i2,i3) = -gm_s / rvec_s - gm_p / sqrt(rvec_p(i1,i2,i3))
      end do
    end do
  end do

  nullify (m1, m2, m3)
  call trace_end
END FUNCTION star_and_planet

END MODULE
