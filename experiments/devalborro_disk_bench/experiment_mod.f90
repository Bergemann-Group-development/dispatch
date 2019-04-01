!===============================================================================
!> Initializiation and update of each ZEUS patch
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE solver_mod
  USE pboundary_mod
  USE potential_mod
  USE parameters_mod
  USE scaling_mod
  implicit none
  private

  type, public, extends(solver_t):: experiment_t
  contains
    procedure :: init => init_experiment
    procedure :: update => update_experiment
  end type

  real(8), parameter :: pi = asin(1.0d0) * 2.0d0

CONTAINS

!===============================================================================
!> Initialise this experiment.
!===============================================================================
SUBROUTINE init_experiment (self)
  class(experiment_t) :: self
  !----------------------------------------------------------------------------
  if (self%mem_allocated) return
  call trace_begin('experiment_t%init')

  self%periodic(2) = .false.
  call self%solver_t%init()
  call initial_conditions(self)
  call potential%init()
  call pboundaries%conditions(self, self%it)

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
  ! update the potential (for a growing or moving planet)
  potential%omega = self%grid_vel(3)
  self%mem(:,:,:,self%idx%phi,self%it,1) = potential%selected(self%mesh, self%time, self%dtime)
  call self%solver_t%update()
  call pboundaries%conditions(self, self%new)

  call trace_end
END SUBROUTINE update_experiment

!===============================================================================
!> Initialise the MHD variables
!===============================================================================
SUBROUTINE initial_conditions (patch)
  class(solver_t), target :: patch
  class(mesh_t), pointer :: m1, m2, m3
  integer :: i1, i2, i3
  real(4) :: z, r, phi, zeta_offset
  real(8) :: sigma0, eint
  real(8), dimension(:), allocatable :: vkep
  !----------------------------------------------------------------------------
  call trace%begin('experiment_mod::initial_conditions')
  
  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)
  associate(q=>patch%mem(:,:,:,:,patch%it,1),idx=>patch%idx)

  ! calculate the Keplerian velocity
  allocate(vkep(m2%gn))
  vkep(:m2%gn-1) = sqrt(experiment_params%mstar * scaling%grav / (m2%centre_nat + patch%mesh(2)%r(:)))
  vkep(m2%gn) = vkep(m2%gn-1)
  sigma0 = 2.0d-3 * experiment_params%mstar / pi / experiment_params%aplanet**2

  select type (m3)
  type is (cylindrical_mesh)
    zeta_offset = m3%zeta
  end select
  do i3=m3%lb,m3%ub
    phi = m3%centre_nat + patch%mesh(3)%r(i3)
    do i2=m2%lb,m2%ub
      r = m2%centre_nat + patch%mesh(2)%r(i2)
      do i1=m1%lb,m1%ub
        z = m1%centre_nat + patch%mesh(1)%r(i1)
        q(i1,i2,i3,idx%d) = sigma0
        eint = (experiment_params%aspect_rat * vkep(i2))**2 * sigma0 / (patch%gamma - 1.0)
        if (patch%kind(1:14) == 'zeus_mhd_patch' .or. patch%kind(1:9) == 'stagger2e') then
          q(i1,i2,i3,idx%e) = eint
        end if
        q(i1,i2,i3,idx%px) = 0.0
        q(i1,i2,i3,idx%py) = 0.0
        q(i1,i2,i3,idx%pz) = vkep(i2) * sqrt(1.0 - experiment_params%aspect_rat**2) * sigma0
        if (size(patch%mem,4) > idx%bx) then
          q(i1,i2,i3,idx%bx) = 0.0
          q(i1,i2,i3,idx%by) = 0.0
          q(i1,i2,i3,idx%bz) = 0.0
        end if
        ! initialise total energy
        if (patch%kind(1:14) == 'zeus_mhd_patch' .or. patch%kind(1:6) == 'ramses') then
          q(i1,i2,i3,idx%et) = eint + 0.5 * 0.25 * sigma0 &
                                    * (q(i1,i2,i3,idx%pz)+q(i1,i2+1,i3,idx%px))**2
        end if
        ! the ZEUS module stores the angular momentum rather than the linear momentum
        if (patch%kind(1:14) == 'zeus_mhd_patch') then
          q(i1,i2,i3,idx%py) = q(i1,i2,i3,idx%py) * m1%h2c(i1)
          q(i1,i2,i3,idx%pz) = q(i1,i2,i3,idx%pz) * m1%h31c(i1) * m2%h32c(i2)
        end if
      end do
    end do
  end do

  deallocate (vkep)
  end associate
  nullify (m1, m2, m3)
  call trace%end
END SUBROUTINE initial_conditions

END MODULE experiment_mod
