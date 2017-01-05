!> =============================================================================
!> Initializiation for the MHD blast experiment
!> The set up draws heavily from the "blast" problem in AZEuS/ZEUS-3D.
!> =============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE mesh_mod
  USE mhd_mod
  implicit none
  private

  type, public :: initial_t
    real :: r0, centre(3), d0(2), p0(2), vx0(2), vy0(2), vz0(2), bx0(2), by0(2), bz0(2)
  contains
    procedure :: init
  end type

  logical, save :: first_time = .true.

CONTAINS

!> =============================================================================
!> Initialise a spherical region with a specific radius wherein the MHD
!> variables differ from the remainder of the domain.
!> This experiment should be able to generate the initial conditions for Sedov
!> and MHD blast (Londrillo & Del Zanna 2000; Gardiner & Stone 2010; Ramsey
!> et al. 2012, etc.) type benchmarks.
!>
!> Namelist parameters:
!>  r0: the radius of the initial spherical "discontinuity"
!>  centre: the coordinates of the centre of the initial discontinuity
!>  d0: the gas density in the ambient [d0(1)] and central regions [d0(2)]
!>  p0: the gas pressure "  "     "    [p0(1)]  "     "       "    [p0(2)]
!>  vx0: the x-velocity in the ambient [vx0(1)] and central regions [vx0(2)]
!>  vy0: the y-velocity  "  "     "    [vy0(1)]  "     "       "    [vy0(2)]
!>  vz0: the z-velocity  "  "     "    [vz0(1)]  "     "       "    [vz0(2)]
!>  bx0: the x-mag. field " "     "    [bx0(1)]  "     "       "    [bx0(2)]
!>  by0: the y-mag. field " "     "    [by0(1)]  "     "       "    [by0(2)]
!>  bz0: the z-mag. field " "     "    [bz0(1)]  "     "       "    [bz0(2)]
!>
!> =============================================================================
SUBROUTINE init (self, patch)
  class(initial_t):: self
  class(mhd_t), target:: patch
  class(mesh_t), pointer:: m1, m2, m3
  integer:: i1, i2, i3
  real(4):: speedsq, bsq
  real(4):: r0=0.125, centre(3)=0.0, d0(2)=1.0, p0(2)=0.6, vx0(2)=0.0, vy0(2)=0.0, &
             vz0(2)=0.0, bx0(2)=0.0, by0(2)=0.0, bz0(2)=0.0
  real(8):: x, y, z, r2
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: f => null()
  namelist / initial_params / r0, centre, d0, p0, vx0, vy0, vz0, bx0, by0, bz0
  !----------------------------------------------------------------------------
  call trace_begin('initial_t%init')

  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, initial_params)
    if (io%master) write (*, initial_params)
  end if
  !$omp end critical (input_cr)

  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)

  f => patch%mem(:,:,:,:,1,1)

  r0 = max(r0, m1%d, m2%d, m3%d)
  self%r0 = r0
  self%centre = centre
  self%d0 = d0
  self%p0 = p0
  self%vx0 = vx0
  self%vy0 = vy0
  self%vz0 = vz0
  self%bx0 = bx0
  self%by0 = by0
  self%bz0 = bz0

  speedsq = sqrt(vx0(1)**2 + vy0(1)**2 + vz0(1)**2)
  bsq = sqrt(bx0(1)**2 + by0(1)**2 + bz0(1)**2)
  ! set up atmosphere/ambient.
  do i3=m3%lb,m3%ub
    z = m3%centre_nat + patch%mesh(3)%r(i3)
    do i2=m2%lb,m2%ub
      y = m2%centre_nat + patch%mesh(2)%r(i2)
      do i1=m1%lb,m1%ub
        x = m1%centre_nat + patch%mesh(1)%r(i1)
        f(i1,i2,i3,didx) = d0(1)

        f(i1,i2,i3,sidxs(1)) = d0(1) * vx0(1)
        f(i1,i2,i3,sidxs(2)) = d0(1) * vy0(1)
        f(i1,i2,i3,sidxs(3)) = d0(1) * vz0(1)
        f(i1,i2,i3,Bidxs(1)) = bx0(1)
        f(i1,i2,i3,Bidxs(2)) = by0(1)
        f(i1,i2,i3,Bidxs(3)) = bz0(1)

        if (patch%kind == 'stagger_mhd_patch') then
          f(i1,i2,i3,Sidx)  = (log(p0(1)) - log(d0(1)) * patch%gamma) &
                            * d0(1) / (patch%gamma - 1.0)
        else if (patch%kind == 'ramses_mhd_patch') then
          f(i1,i2,i3,etidx) = p0(1) + 0.5 * d0(1) * speedsq + 0.5 * bsq
        else
          f(i1,i2,i3,e1idx) = p0(1)
        end if

        if (patch%kind == 'zeus_mhd_patch') then
          ! the ZEUS module stores the angular momentum rather than the linear momentum
          f(i1,i2,i3,sidxs(2)) = f(i1,i2,i3,sidxs(2)) * m1%h2c(i1)
          f(i1,i2,i3,sidxs(3)) = f(i1,i2,i3,sidxs(3)) * m1%h31c(i1) * m2%h32c(i2)
        end if
      end do
    end do
  end do

  ! set up central region.
  speedsq = sqrt(vx0(2)**2 + vy0(2)**2 + vz0(2)**2)
  bsq = sqrt(bx0(2)**2 + by0(2)**2 + bz0(2)**2)
  do i3=m3%lb,m3%ub
    z = m3%centre_nat + patch%mesh(3)%r(i3)
    do i2=m2%lb,m2%ub
      y = m2%centre_nat + patch%mesh(2)%r(i2)
      do i1=m1%lb,m1%ub
        x = m1%centre_nat + patch%mesh(1)%r(i1)
        r2 = (x - centre(1))**2 + (y - centre(2))**2 + (z - centre(3))**2
        if (r2 .le. real(r0,kind=8)**2) then
          f(i1,i2,i3,didx) = d0(2)

          f(i1,i2,i3,sidxs(1)) = d0(2) * vx0(2)
          f(i1,i2,i3,sidxs(2)) = d0(2) * vy0(2)
          f(i1,i2,i3,sidxs(3)) = d0(2) * vz0(2)
          f(i1,i2,i3,Bidxs(1)) = bx0(2)
          f(i1,i2,i3,Bidxs(2)) = by0(2)
          f(i1,i2,i3,Bidxs(3)) = bz0(2)

          if (patch%kind == 'stagger_mhd_patch') then
            f(i1,i2,i3,Sidx)  = (log(p0(2)) - log(d0(2)) * patch%gamma) &
                              * d0(2) / (patch%gamma - 1.0)
          else if (patch%kind == 'ramses_mhd_patch') then
            f(i1,i2,i3,etidx) = p0(2) + 0.5 * d0(2) * speedsq + 0.5 * bsq
          else
            f(i1,i2,i3,e1idx) = p0(2)
          end if

          if (patch%kind == 'zeus_mhd_patch') then
            ! the ZEUS module stores the angular momentum rather than the linear momentum
            f(i1,i2,i3,sidxs(2)) = f(i1,i2,i3,sidxs(2)) * m1%h2c(i1)
            f(i1,i2,i3,sidxs(3)) = f(i1,i2,i3,sidxs(3)) * m1%h31c(i1) * m2%h32c(i2)
          end if
        end if
      end do
    end do
  end do

  nullify (m1, m2, m3, f)
  call trace_end
END SUBROUTINE init

END MODULE initial_mod
