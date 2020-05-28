!===============================================================================
!> Direct force of gravity
!===============================================================================
MODULE force_mod
  USE io_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  !USE vector_mod
  USE mesh_mod
  USE planet_mod
  USE star_mod
  implicit none
  private
  type, public:: force_t
    procedure(star_and_planet), pointer:: selected=>null()
  contains
    procedure init
  end type
  type(force_t), public:: force
CONTAINS

!===============================================================================
!> Select type of force -- could be more than one planet, for example
!===============================================================================
SUBROUTINE init (self, solver)
  class(force_t):: self
  character(len=*):: solver
  self%selected => star_and_planet
END SUBROUTINE

!===============================================================================
!> Force of gravity from a star and a planet
!===============================================================================
FUNCTION star_and_planet (self, time, d, p, Ux, Uy, Uz, m) RESULT (ff)
  class(force_t)                                :: self
  real(8)                                       :: time
  real, dimension(:,:,:)  , pointer             :: Ux, Uy, Uz
  real, dimension(:,:,:)  , pointer, contiguous :: d
  real, dimension(:,:,:,:), pointer             :: p
  class(mesh_t), dimension(:), pointer          :: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3)    :: ff
  !.............................................................................
  integer:: ix, iy, iz, iv
  real(8):: kx, ky, kz, rs(3), rp(3), rs1, rp1
  real(8):: x, y, z, xs, ys, zs, xp, yp, zp
  !-----------------------------------------------------------------------------
  call trace_begin('force_t:single_solenoidal')
  do iv=1,3
    do iz=m(3)%lb,m(3)%ub
      z = m(3)%p + m(3)%r(iz) + m(3)%h(iv)*m(3)%d
      zs = z - star%position(3)
      zp = z - planet%position(3)
      rs(3) = zs
      rp(3) = zp
      do iy=m(2)%lb,m(2)%ub
        y = m(2)%p + m(2)%r(iy) + m(2)%h(iv)*m(2)%d
        yp = y - planet%position(2)
        ys = y - star%position(2)
        rs(2) = ys
        rp(2) = yp
        do ix=m(1)%lb,m(1)%ub
          x = m(1)%p + m(1)%r(ix) + m(1)%h(iv)*m(1)%d
          xp = x - planet%position(1)
          xs = x - star%position(1)
          rs(1) = ys
          rp(1) = yp
          rp1 = sqrt(xp**2+yp**2+zp**2)
          rs1 = sqrt(xs**2+ys**2+zs**2)
          ff(ix,iy,iz,iv) = &
            -star%mass  *scaling%grav*rs(iv)/rs1**3 &
            -planet%mass*scaling%grav*rp(iv)/rp1**3
        end do
      end do
    end do
  end do
  call trace_end
END FUNCTION star_and_planet

END MODULE
