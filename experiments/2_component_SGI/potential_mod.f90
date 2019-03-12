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
    procedure(star), pointer:: selected=>null()
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
  self%selected => star

  call trace_end
END SUBROUTINE

!===============================================================================
!> Returns the combined gravitational potential of a star and planet in
!> cylindrical coords.
!===============================================================================
FUNCTION star (self, m) RESULT (gp)
  class(potential_t):: self
  class(mesh_t), dimension(:), pointer :: m
  class(mesh_t), pointer :: m1, m2, m3
  real(kind=KindScalarVar), dimension(m(1)%gn,m(2)%gn,m(3)%gn):: gp
  !.............................................................................
  integer:: i1, i2, i3
  real(8):: z, r, phi, z_s, r_s, gm_s, rvec_s
  !-----------------------------------------------------------------------------
  call trace%begin('potential_mod::star')

  m1 => m(1)
  m2 => m(2)
  m3 => m(3)

  gm_s = scaling%grav * experiment_params%mstar
  ! for safety
  !allocate(gp(m1%gn,m2%gn,m3%gn))
  gp = 0.0
  !
  do i3=m3%lb,m3%ub
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
  call trace%end
END FUNCTION star

END MODULE
