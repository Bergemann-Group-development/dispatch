!> =============================================================================
!> Radiative transfer shadow test.
!> http://adsabs.harvard.edu/abs/2015A&A...574A..81R
!> =============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  USE mesh_mod
  USE mhd_mod
  implicit none
  private

  real(8), parameter:: pi2=8d0*atan(1d0)

  type, public:: initial_t
    real:: clump_pos(3), clump_a, clump_b, d0, t0, clump_dfact
  contains
    procedure :: init
    procedure :: initial_conditions
  end type

  logical, save:: first_time=.true.

CONTAINS
!> =============================================================================
!> Read namelist
!> =============================================================================
SUBROUTINE init (self)
  class(initial_t):: self
  real, save:: clump_pos(3)=0.0, clump_a=0.1, clump_b=0.06, d0=1.0, t0=290.0, clump_dfact=10.0
  namelist /initial_params/ clump_pos, clump_a, clump_b, d0, t0, clump_dfact
  !----------------------------------------------------------------------------
  call trace_begin('initial_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, initial_params)
    if (io%master) write (*, initial_params)
    call scaling%init()
  end if
  !$omp end critical (input_cr)

  self%clump_pos = clump_pos
  self%clump_a = clump_a
  self%clump_b = clump_b
  self%d0 = d0/scaling%d
  self%t0 = t0/scaling%temp
  self%clump_dfact = clump_dfact

  call trace_end
END SUBROUTINE init

SUBROUTINE initial_conditions (self, patch)
!> =============================================================================
!> Set the initial conditions
!> =============================================================================
  class(initial_t) :: self
  class(mhd_t), target :: patch
  class(mesh_t), pointer :: mx, my, mz
  real, dimension(:,:,:,:), pointer :: f
  integer :: ix, iy, iz
  real :: p0, q, x, y, z
  !----------------------------------------------------------------------------
  call trace_begin('initial_t::initial_conditions')

  mx => patch%mesh(1)
  my => patch%mesh(2)
  mz => patch%mesh(3)
  f => patch%mem(:,:,:,:,1,1)

  ! note: I've assumed a unit system where k_B / mu * m_pro is 1.
  p0 = self%d0 * self%t0
  do iz=mz%lb,mz%ub
    z = mz%p + mz%r(iz)
    do iy=my%lb,my%ub
      y = my%p  + my%r(iy)
      do ix=mx%lb,mx%ub
        x = mx%p  + mx%r(ix)

        ! currently, the y and z coordinates are treated symmetrically; this
        ! leads to a cylindrically symmetric problem in 3-D.
        q = ((x - self%clump_pos(1)) / self%clump_a)**2 &
          + ((y - self%clump_pos(2)) / self%clump_b)**2 &
          + ((z - self%clump_pos(3)) / self%clump_b)**2

        ! the if-then-else is to prevent overflow when (x,y,z) is very from zero.
        if (q > 88.0) then
          f(ix,iy,iz,1) = self%d0
        else
          f(ix,iy,iz,1) = self%d0 + (self%clump_dfact * self%d0 - self%d0) &
                                  / (1.0 + exp(q))
        end if

        ! the clump is initially in pressure equilibrium with the ambient medium
        f(ix,iy,iz,2) = (log(p0) - log(f(ix,iy,iz,1)) * patch%gamma) &
                      * f(ix,iy,iz,1) / (patch%gamma - 1.0)


        ! momentum components
        f(ix,iy,iz,3) = 0.0
        f(ix,iy,iz,4) = 0.0
        f(ix,iy,iz,5) = 0.0
        ! magnetic field components
        f(ix,iy,iz,6) = 0.0
        f(ix,iy,iz,7) = 0.0
        f(ix,iy,iz,8) = 0.0
      end do
    end do
  end do

  nullify (mx, my, mz, f)
  call trace_end
END SUBROUTINE initial_conditions

END MODULE
