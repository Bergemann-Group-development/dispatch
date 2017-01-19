!> =============================================================================
!> Initializiation and update of each ZEUS patch
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE scaling_mod
  USE mesh_mod
  USE timestep_mod
  USE mhd_mod
  USE force_mod
  USE boundary_mod
  USE parameters_mod
  implicit none
  private

  type, public, extends(mhd_t):: experiment_t
  contains
    procedure :: init => init_experiment
    procedure :: update => update_experiment
  end type

  real(8), parameter :: pi = asin(1.0d0) * 2.0d0

CONTAINS

!> =============================================================================
!> Initialise this experiment.
!> =============================================================================
SUBROUTINE init_experiment (self)
  class(experiment_t) :: self
  !----------------------------------------------------------------------------
  if (self%mem_allocated) return
  call trace_begin('experiment_t%init')

  self%mhd = .false.
  self%staggered = .true.

  call self%mhd_t%init()
  call initial_conditions(self)
  call boundary%conditions(self)
  call force%init(self)

  call trace_end
END SUBROUTINE init_experiment

!> =============================================================================
!> Update procedure for this experiment.
!> =============================================================================
SUBROUTINE update_experiment (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')

  call boundary%conditions(self)
  self%gp%q = force%selected(self%mesh, self%time, self%dtime)
  call self%mhd_t%update()
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call boundary%conditions(self)

  call trace_end
END SUBROUTINE update_experiment

!> =============================================================================
!> Initialise the MHD variables
!> The current implementation implements MMSN (Minimum Mass Solar Nebula) conditions.
!> =============================================================================
SUBROUTINE initial_conditions (patch)
  class(mhd_t), target :: patch
  class(mesh_t), pointer :: m1, m2, m3
  integer :: i1, i2, i3
  real(8) :: x1, x2, x3, cylr, cosphi, vr, vphi, vkep, tgas0
  real(8) :: pos_cyl(3)
  !----------------------------------------------------------------------------
  call trace_begin('initial_mod::initial_conditions')
  
  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)
  
  pos_cyl = CartesianToCylindrical(patch%position(1), patch%position(2), patch%position(3))

  do i3=m3%lb,m3%ub
    x3 = m3%centre_nat + m3%r(i3)
    do i2=m2%lb,m2%ub
      x2 = m2%centre_nat + m2%r(i2)
      do i1=m1%lb,m1%ub
        x1 = m1%centre_nat + m1%r(i1)

        select type (m1)
        type is (Cartesian_mesh)
          cylr = sqrt(x1**2 + x2**2)
        type is (cylindrical_mesh)
          cylr = sqrt(pos_cyl(2)**2 + x2**2 - 2.0d0 * x2 * pos_cyl(2) * cos(pi - x3))
        type is (spherical_mesh)
          cylr = sqrt(pos_cyl(2)**2 + x1**2 - 2.0d0 * x1 * pos_cyl(2) * cos(pi - x3))
        end select

        patch%d%q (i1,i2,i3)  = MMSN_surface_density(cylr)
        ! the ideal gas constant, k / mu / m_p, is assumed incorporated into `scaling%temp`.
        tgas0 = MMSN_temperature(cylr)
        patch%e1%q(i1,i2,i3)  = tgas0 * patch%d%q(i1,i2,i3) &
                              / (patch%gamma - 1.0)

        ! the ZEUS module stores the angular momentum rather than the linear momentum
        ! note: to obtain consistent results across mesh types, this step is done in the
        ! same cylindrical coord. system as used to obtain the physical variables (rather
        ! than the current patch coord. system).
        vkep = Keplerian_speed(cylr)

        select type (m1)
        type is (cylindrical_mesh)
          patch%s%x%q(i1,i2,i3) = 0.0
          patch%s%y%q(i1,i2,i3) = 0.0
          patch%s%z%q(i1,i2,i3) = vkep * patch%d%q(i1,i2,i3) * m1%h31c(i1) * m2%h32c(i2)
        type is (spherical_mesh)
          vr   = -vkep * sin(x3)
          vphi = -vkep * cos(x3)
          patch%s%x%q(i1,i2,i3) = vr * patch%d%q(i1,i2,i3)
          patch%s%y%q(i1,i2,i3) = 0.0d0
          patch%s%z%q(i1,i2,i3) = vphi * patch%d%q(i1,i2,i3) * m1%h31c(i1) * m2%h32c(i2)
        end select

        ! no magnetic fields in this set up
        patch%b%x%q(i1,i2,i3) = 0.0
        patch%b%y%q(i1,i2,i3) = 0.0
        patch%b%z%q(i1,i2,i3) = 0.0
        ! the gravitational potential is assigned in `force_mod`.
      end do
    end do
  end do

  nullify (m1, m2, m3)
  call trace_end
END SUBROUTINE initial_conditions

END MODULE experiment_mod
