!===============================================================================
!> Direct force of gravity
!===============================================================================
MODULE potential_mod
  USE io_mod
  USE trace_mod
  USE scaling_mod
  USE mesh_mod
  USE parameters_mod
  USE kinds_mod
  implicit none
  private

  type, public:: potential_t
    real(8):: omega
    procedure(star_and_planet), pointer:: selected=>null()
  contains
    procedure init
  end type
  type(potential_t), public:: potential

CONTAINS
!===============================================================================
!> Select type of force -- could be more than one planet, for example
!===============================================================================
SUBROUTINE init (self)
  class(potential_t) :: self
  !-----------------------------------------------------------------------------
  call trace_begin('potential_t%init')

  ! select implementation for the potential source term
  self%selected => star_and_planet

  call trace_end
END SUBROUTINE

!===============================================================================
!> Returns the combined gravitational potential of a star and planet in
!> cylindrical coords.
!===============================================================================
FUNCTION star_and_planet (self, m, time, dt) RESULT (gp)
  class(potential_t):: self
  class(mesh_t), dimension(:), pointer :: m
  class(mesh_t), pointer :: m1, m2, m3
  real(8), intent(in) :: time, dt
  real(kind=KindScalarVar), dimension(m(1)%gn,m(2)%gn,m(3)%gn):: gp
  !.............................................................................
  integer:: i1, i2, i3
  real(8):: z, r, phi, z_s, r_s, gm_s, gm_p, rvec_s, rvec_p, soft, radius_p_i, r_dot_rp
  !-----------------------------------------------------------------------------
  call trace_begin('potential_mod::star_and_planet')

  m1 => m(1)
  m2 => m(2)
  m3 => m(3)

  ! update the planet mass
  call planet%grow(time)
  ! udpate the planet position
  call planet%move(dt, self%omega)

  gm_s = scaling%grav * experiment_params%mstar
  gm_p = scaling%grav * planet%mass
  soft = experiment_params%gpp_soften * experiment_params%aspect_rat * experiment_params%aplanet
  radius_p_i = 1.0d0 / (planet%position(1)**2 + planet%position(2)**2)
  ! for safety
  gp = 0.0
  !
  do i3=m3%lb,m3%ub
    phi = m3%centre_nat + m3%r(i3)
    do i2=m2%lb,m2%ub
      r = m2%centre_nat + m2%r(i2)
      r_s = r - 0.0
      do i1=m1%lb,m1%ub
        z = m1%centre_nat + m1%r(i1)
        z_s = z - 0.0

        rvec_s = sqrt(z_s**2 + r_s**2)
        ! apply the law of cosines to find the distance between the current cell and the planet.
        rvec_p = planet%position(2)**2 + r**2 &
               - 2.0d0 * planet%position(2) * r * cos(abs(phi - planet%position(3)))
        ! dot product between current cell position and the planet's position.
        r_dot_rp = z * planet%position(1) + r * planet%position(2) * cos(phi - planet%position(3))
        ! 1st: stellar potential; 2nd: planetary point potential;
        ! 3rd: indirect term due to the centre-of-mass not coinciding exactly with the position of the star.
        gp(i1,i2,i3) = -gm_s / rvec_s - gm_p / sqrt(rvec_p + soft**2) &
                     + gm_p * radius_p_i * r_dot_rp
      end do
    end do
  end do

  nullify (m1, m2, m3)
  call trace_end
END FUNCTION star_and_planet

END MODULE
