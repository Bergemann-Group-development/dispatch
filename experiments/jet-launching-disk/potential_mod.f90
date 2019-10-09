!===============================================================================
!> Direct force of gravity
!===============================================================================
MODULE potential_mod
  USE io_mod
  USE trace_mod
  USE scaling_mod
  USE mesh_mod
  USE mpi_mod
  USE index_mod
  USE kinds_mod
  implicit none
  private

  type, public:: potential_t
  contains
    procedure init
  end type
  type(potential_t), public:: potential

CONTAINS
!===============================================================================
!> Select type of force -- could be more than one planet, for example
!===============================================================================
SUBROUTINE init (self, m, idx, star_mass, ff)
  class(potential_t) :: self
  type(index_t):: idx
  class(mesh_t), dimension(:), pointer :: m
  real(8):: star_mass
  real(kind=KindScalarVar), dimension(:,:,:,:):: ff
  !-----------------------------------------------------------------------------
  call trace%begin('potential_t%init')

  ! select implementation for the potential source term
  call central_star_potential(m, idx, star_mass, ff)

  call trace%end()
END SUBROUTINE

!===============================================================================
!> Returns the combined gravitational potential of a star and planet in
!> cylindrical coords.
!===============================================================================
SUBROUTINE central_star_potential(m, idx, star_mass, ff)
  class(mesh_t), dimension(:), pointer :: m
  type(index_t):: idx
  real(8):: star_mass
  real(kind=KindScalarVar), dimension(:,:,:,:):: ff
  !.............................................................................
  class(mesh_t), pointer :: m1, m2, m3
  integer:: i1, i2, i3
  real(8):: r_sph, gm_s, x1, x2
  !-----------------------------------------------------------------------------
  call trace%begin('potential_mod::central_star_potential')

  ! fool-proofing
  if (idx%phi <= 0) then
    call mpi%abort('There is no slot in mem available for phi! Abort!')
  end if

  m1 => m(1)
  m2 => m(2)
  m3 => m(3)

  gm_s = scaling%grav * star_mass
  do i3=m3%lb,m3%ub
    do i2=m2%lb,m2%ub
      x2 = m2%centre_nat + m2%r(i2)
      do i1=m1%lb,m1%ub
        x1 = m1%centre_nat + m1%r(i1)
        if (m1%type == 'spherical_mesh') then
          r_sph = abs(x1)
        else if (m1%type == 'cylindrical_mesh') then
          r_sph = sqrt(x1**2 + x2**2)
        end if
        ff(i1,i2,i3,idx%phi) = -gm_s / r_sph
      end do
    end do
  end do

  nullify (m1, m2, m3)
  call trace%end()
END SUBROUTINE central_star_potential

END MODULE
