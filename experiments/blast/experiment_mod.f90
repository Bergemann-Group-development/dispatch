!===============================================================================
!> Initializiation and update of each ZEUS patch
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE mesh_mod
  USE pboundary_mod
  USE solver_mod
  USE pboundary_mod
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
  call self%solver_t%init()
  call initial_conditions(self)
  call pboundaries%conditions(self, self%it)
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

  call pboundaries%conditions(self, self%it)
  call self%solver_t%update()
  call pboundaries%conditions(self, self%new)

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
!>  e0: the gas thermal en. "     "    [e0(1)]  "     "       "    [e0(2)]
!>  vx0: the x-velocity in the ambient [vx0(1)] and central regions [vx0(2)]
!>  vy0: the y-velocity  "  "     "    [vy0(1)]  "     "       "    [vy0(2)]
!>  vz0: the z-velocity  "  "     "    [vz0(1)]  "     "       "    [vz0(2)]
!>  bx0: the x-mag. field " "     "    [bx0(1)]  "     "       "    [bx0(2)]
!>  by0: the y-mag. field " "     "    [by0(1)]  "     "       "    [by0(2)]
!>  bz0: the z-mag. field " "     "    [bz0(1)]  "     "       "    [bz0(2)]
!>
!===============================================================================
SUBROUTINE initial_conditions (patch)
  class(experiment_t):: patch
  class(mesh_t), pointer:: m1, m2, m3
  integer:: i1, i2, i3, iostat
  real(4):: speedsq, bsq, dxmax
  real(4):: r0=0.125, centre(3)=0.0, d0(2)=1.0, e0(2)=0.6, vx0(2)=0.0, vy0(2)=0.0, &
             vz0(2)=0.0, bx0(2)=0.0, by0(2)=0.0, bz0(2)=0.0
  real(8):: x, y, z, r2, e2p
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: f
  namelist / experiment_params / r0, centre, d0, e0, vx0, vy0, vz0, bx0, by0, bz0
  !----------------------------------------------------------------------------
  call trace_begin('initial_t%init')

  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    call io%namelist_warning ('experiment_params', iostat)
    if (io%master) write (io%output, experiment_params)
  end if
  call io%lock%unset

  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)

  f => patch%mem(:,:,:,:,1,1)

  dxmax = maxval(patch%ds,mask=(patch%n>1))
  r0 = max(r0, dxmax)

  speedsq = vx0(1)**2 + vy0(1)**2 + vz0(1)**2
  bsq = bx0(1)**2 + by0(1)**2 + bz0(1)**2
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
        if (patch%nv > 5) then
          f(i1,i2,i3,patch%idx%bx) = bx0(1)
          f(i1,i2,i3,patch%idx%by) = by0(1)
          f(i1,i2,i3,patch%idx%bz) = bz0(1)
        end if

        if (trim(patch%kind) == 'stagger2_mhd_patch') then
          f(i1,i2,i3,patch%idx%s)  = (log(e0(1)*e2p) - log(d0(1)) * patch%gamma) &
                            * d0(1) / (patch%gamma - 1.0)
        else if (trim(patch%kind(1:9)) == 'stagger2e') then
          f(i1,i2,i3,patch%idx%e) = e0(1)
        else if (trim(patch%kind(1:6)) == 'ramses') then
          f(i1,i2,i3,patch%idx%e) = e0(1) + 0.5 * d0(1) * speedsq + 0.5 * bsq
        else if (trim(patch%kind) == 'zeus_mhd_patch') then
          f(i1,i2,i3,patch%idx%e) = e0(1)

          ! the ZEUS module stores the angular momentum rather than the linear momentum
          f(i1,i2,i3,patch%idx%py) = f(i1,i2,i3,patch%idx%py) * m1%h2c(i1)
          f(i1,i2,i3,patch%idx%pz) = f(i1,i2,i3,patch%idx%pz) * m1%h31c(i1) * m2%h32c(i2)
        end if
      end do
    end do
  end do

  ! set up central region.
  speedsq = vx0(2)**2 + vy0(2)**2 + vz0(2)**2
  bsq = bx0(2)**2 + by0(2)**2 + bz0(2)**2
  e2p = (1d0 - patch%gamma)
  do i3=m3%lb,m3%ub
    z = m3%centre_nat + patch%mesh(3)%r(i3)
    do i2=m2%lb,m2%ub
      y = m2%centre_nat + patch%mesh(2)%r(i2)
      do i1=m1%lb,m1%ub
        x = m1%centre_nat + patch%mesh(1)%r(i1)
        r2 = (x - centre(1))**2 + (y - centre(2))**2 + (z - centre(3))**2
        if (r2 .le. real(r0,kind=8)**2) then
          f(i1,i2,i3,patch%idx%d) = d0(2)
          f(i1,i2,i3,patch%idx%px) = d0(2) * vx0(2)
          f(i1,i2,i3,patch%idx%py) = d0(2) * vy0(2)
          f(i1,i2,i3,patch%idx%pz) = d0(2) * vz0(2)
          if (patch%nv > 5) then
            f(i1,i2,i3,patch%idx%bx) = bx0(2)
            f(i1,i2,i3,patch%idx%by) = by0(2)
            f(i1,i2,i3,patch%idx%bz) = bz0(2)
          end if

          if (trim(patch%kind) == 'stagger2_mhd_patch') then
            f(i1,i2,i3,patch%idx%s)  = (log(e0(2)) - log(d0(2)) * patch%gamma) &
                              * d0(2) / (patch%gamma - 1.0)
          else if (trim(patch%kind) == 'stagger2e_mhd_patch') then
            f(i1,i2,i3,patch%idx%e) = e0(2)
          else if (trim(patch%kind(1:6)) == 'ramses') then
            f(i1,i2,i3,patch%idx%e) = e0(2) + 0.5 * d0(2) * speedsq + 0.5 * bsq
          else if (trim(patch%kind) == 'zeus_mhd_patch') then
            f(i1,i2,i3,patch%idx%e) = e0(2)

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
