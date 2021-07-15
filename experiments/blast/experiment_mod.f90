!===============================================================================
!> Initializiation and update of each ZEUS patch
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE mesh_mod
  USE solver_mod
  USE eos_mod
  implicit none
  private

  type, public, extends(solver_t):: experiment_t
  contains
    procedure :: init => init_experiment
    procedure :: update => update_experiment
  end type

  logical, save :: first_time = .true.

CONTAINS

!===============================================================================
!> Initialise this experiment.
!===============================================================================
SUBROUTINE init_experiment (self)
  class(experiment_t) :: self
  !----------------------------------------------------------------------------
  if (self%mem_allocated) return
  call trace_begin('experiment_t%init')
  self%mhd=.true.
  call self%solver_t%init()
  call initial_conditions(self)
  call self%output
  call trace_end
END SUBROUTINE init_experiment

!===============================================================================
!> Update procedure for this experiment.
!===============================================================================
SUBROUTINE update_experiment (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')

  call self%solver_t%update()

  call trace_end
END SUBROUTINE update_experiment

!===============================================================================
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
!===============================================================================
SUBROUTINE initial_conditions (patch)
  class(experiment_t), target:: patch
  class(mesh_t), pointer:: m1, m2, m3
  integer:: i1, i2, i3
  real(4), save:: speedsq, bsq, dxmax
  real(4), save:: r0=0.125, centre(3)=0.0, d0(2)=1.0, p0(2)=0.6, vx0(2)=0.0, vy0(2)=0.0, &
             vz0(2)=0.0, bx0(2)=0.0, by0(2)=0.0, bz0(2)=0.0
  real(8):: x, y, z, r2
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: f
  integer:: iostat
  integer, save:: case
  namelist / experiment_params / r0, centre, d0, p0, vx0, vy0, vz0, bx0, by0, bz0,case
  !----------------------------------------------------------------------------
  call trace_begin('initial_t%init')

  call io%lock%set('initial_t%init')
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    call io%namelist_warning ('experiment_params', iostat)
    if (io%master) write (io%output, experiment_params)
  end if
  call io%lock%unset('initial_t%init')

  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)

  f => patch%mem(:,:,:,:,1,1)
  dxmax = maxval(patch%ds,mask=(patch%n>1))
  r0 = max(r0, dxmax)
  select case(case)
  case(2)
    vz0=0.0
    bz0=0.0
  case(3)
    vy0=0.0
    by0=0.0
  case(4)
    vx0=0.0
    bx0=0.0
  case default
  end select

  speedsq = vx0(1)**2 + vy0(1)**2 + vz0(1)**2
  if (patch%idx%bx > 0) then
    bsq = bx0(1)**2 + by0(1)**2 + bz0(1)**2
  else
    bsq = 0.0
  end if
  ! set up atmosphere/ambient.
  do i3=m3%lb,m3%ub
    z = m3%centre_nat + patch%mesh(3)%r(i3)
    do i2=m2%lb,m2%ub
      y = m2%centre_nat + patch%mesh(2)%r(i2)
      do i1=m1%lb,m1%ub
        x = m1%centre_nat + patch%mesh(1)%r(i1)

        f(i1,i2,i3,patch%idx%d)  = d0(1)
        f(i1,i2,i3,patch%idx%px) = d0(1) * vx0(1)
        f(i1,i2,i3,patch%idx%py) = d0(1) * vy0(1)
        f(i1,i2,i3,patch%idx%pz) = d0(1) * vz0(1)
        if (patch%idx%bx > 0) then
          f(i1,i2,i3,patch%idx%bx) = bx0(1)
          f(i1,i2,i3,patch%idx%by) = by0(1)
          f(i1,i2,i3,patch%idx%bz) = bz0(1)
        end if

        if (patch%idx%s>0) then
          f(i1,i2,i3,patch%idx%s)  = (eos%s0+log(p0(1)*d0(1)**-eos%gamma)) * d0(1) / (eos%gamma - 1.0)
        else if (patch%idx%e>0) then
          f(i1,i2,i3,patch%idx%e) = p0(1)/(eos%gamma-1)
          if (patch%kind(1:6) == 'ramses') then
            f(i1,i2,i3,patch%idx%e) = f(i1,i2,i3,patch%idx%e) + 0.5 *d0(1)  * speedsq + 0.5 * bsq
          end if
        end if
        if (trim(patch%kind) == 'zeus_mhd_patch') then
          ! the ZEUS module stores the angular momentum rather than the linear momentum
          f(i1,i2,i3,patch%idx%py) = f(i1,i2,i3,patch%idx%py) * m1%h2c(i1)
          f(i1,i2,i3,patch%idx%pz) = f(i1,i2,i3,patch%idx%pz) * m1%h31c(i1) * m2%h32c(i2)
        end if
      end do
    end do
  end do

  ! set up central region.
  speedsq = vx0(2)**2 + vy0(2)**2 + vz0(2)**2
  if (patch%idx%bx > 0) then
    bsq = bx0(2)**2 + by0(2)**2 + bz0(2)**2
  else
    bsq = 0.0
  end if
  do i3=m3%lb,m3%ub
    z = m3%centre_nat + patch%mesh(3)%r(i3)
    do i2=m2%lb,m2%ub
      y = m2%centre_nat + patch%mesh(2)%r(i2)
      do i1=m1%lb,m1%ub
        x = m1%centre_nat + patch%mesh(1)%r(i1)
        select case(case)
        case (2)
          r2 = (x - centre(1))**2 + (y - centre(2))**2
        case (3)
          r2 = (x - centre(1))**2 + (z - centre(3))**2
        case (4)
          r2 = (y - centre(2))**2 + (z - centre(3))**2
        case default
          r2 = (x - centre(1))**2 + (y - centre(2))**2 + (z - centre(3))**2
        end select
        if (r2 .le. real(r0,kind=8)**2) then
          f(i1,i2,i3,patch%idx%d) = d0(2)
          f(i1,i2,i3,patch%idx%px) = d0(2) * vx0(2)
          f(i1,i2,i3,patch%idx%py) = d0(2) * vy0(2)
          f(i1,i2,i3,patch%idx%pz) = d0(2) * vz0(2)
          if (patch%idx%bx > 0) then
            f(i1,i2,i3,patch%idx%bx) = bx0(2)
            f(i1,i2,i3,patch%idx%by) = by0(2)
            f(i1,i2,i3,patch%idx%bz) = bz0(2)
          end if

          if (patch%idx%s>0) then
            f(i1,i2,i3,patch%idx%s)  = (eos%s0+log(p0(2)*d0(2)**-eos%gamma)) * d0(2) / (eos%gamma - 1.0)
          else if (patch%idx%e>0) then
            f(i1,i2,i3,patch%idx%e) = p0(2)/(eos%gamma-1)
            if (patch%kind(1:6) == 'ramses') then
              f(i1,i2,i3,patch%idx%e) = f(i1,i2,i3,patch%idx%e) + 0.5 * d0(2) * speedsq + 0.5 * bsq
            end if
          end if
          if (trim(patch%kind) == 'zeus_mhd_patch') then
            ! the ZEUS module stores the angular momentum rather than the linear momentum
            f(i1,i2,i3,patch%idx%py) = f(i1,i2,i3,patch%idx%py) * m1%h2c(i1)
            f(i1,i2,i3,patch%idx%pz) = f(i1,i2,i3,patch%idx%pz) * m1%h31c(i1) * m2%h32c(i2)
          end if
        end if
      end do
    end do
  end do

  nullify (m1, m2, m3, f)
  call trace_end
END SUBROUTINE initial_conditions

END MODULE experiment_mod
