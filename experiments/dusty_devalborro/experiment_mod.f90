!> =============================================================================
!> Initializiation and update of each ZEUS patch
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE solver_mod
  USE pboundary_mod
  USE potential_mod
  USE parameters_mod
  USE scaling_mod
  USE index_mod
  USE download_mod
  implicit none
  private

  type, public, extends(solver_t):: experiment_t
  contains
    procedure :: init => init_experiment
    procedure :: update => update_experiment
    procedure :: dnload
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
  self%nv = experiment_params%nv                        ! Number of total variables 
  self%nc = experiment_params%nc                        ! Number of fluid components
  call self%mhd_t%init()
  call initial_conditions(self)
  call potential%init()
  call pboundaries%conditions(self, self%it)

  call trace_end
END SUBROUTINE init_experiment

!> =============================================================================
!> Update procedure for this experiment.
!> =============================================================================
SUBROUTINE update_experiment (self)
  class(experiment_t):: self

  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call pboundaries%conditions(self, self%it)
  ! update the potential (for a growing or moving planet)
  potential%omega = self%grid_vel(3)
  call self%mhd_t%update()
  call pboundaries%conditions(self, self%new)
  call trace_end
END SUBROUTINE update_experiment

!> =============================================================================
!> Initialise the MHD variables
!> =============================================================================
SUBROUTINE initial_conditions (patch)
  class(solver_t), target :: patch
  class(mesh_t), pointer :: m1, m2, m3
  real, dimension(:,:,:,:), pointer:: p, pd
  real, dimension(:,:,:),   pointer:: d, px , py , pz, dd, pxd, pyd, pzd
  integer :: i1, i2, i3
  real(4) :: z, r, phi, zeta_offset
  real(8) :: sigma0, P0
  real(8), dimension(:), allocatable :: vkep
  !----------------------------------------------------------------------------
  call trace_begin('experiment_mod::initial_conditions')

  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)

  d  => patch%mem(:,:,:,patch%idx%d        ,patch%it,1)
  p  => patch%mem(:,:,:,patch%idx%px:patch%idx%pz,patch%it,1)
  px => p(:,:,:,1); py => p(:,:,:,2); pz => p(:,:,:,3);

  if (patch%nc==2) then
    patch%no_pressure = 2
    dd  => patch%mem(:,:,:,patch%idx%d+5,patch%it,1)
    pd  => patch%mem(:,:,:,patch%idx%px+5:patch%idx%pz+5,patch%it,1)
    pxd => pd(:,:,:,1); pyd => pd(:,:,:,2); pzd => pd(:,:,:,3);
  end if

  ! calculate the Keplerian velocity
  allocate(vkep(m2%gn))
  vkep(:) = sqrt(experiment_params%mstar * scaling%grav / (m2%centre_nat + patch%mesh(2)%r(:)))
  sigma0 = 2.0d-3 * experiment_params%mstar / pi / experiment_params%aplanet**2
  select type (m3)
  type is (cylindrical_mesh)
    zeta_offset = m3%zeta
  end select
  do i3=m3%lb,m3%ub
    phi = m3%centre_nat + m3%r(i3)
    do i2=m2%lb,m2%ub
      r = m2%centre_nat + m2%r(i2)
      do i1=m1%lb,m1%ub
        z = m1%centre_nat + m1%r(i1)
        d (i1,i2,i3) = sigma0
        P0           = (experiment_params%aspect_rat * vkep(i2))**2 * sigma0
        px(i1,i2,i3) = 0.0
        py(i1,i2,i3) = 0.0
        pz(i1,i2,i3) = vkep(i2) * sigma0!* (1.0 - experiment_params%aspect_rat**2)
        ! the stagger2/multifluid module stores the angular momentum rather than the linear momentum
        py(i1,i2,i3) = py(i1,i2,i3) * m1%h2c(i1)
        pz(i1,i2,i3) = pz(i1,i2,i3) * m1%h31c(i1) * m2%h32c(i2)

        if      (trim(patch%kind) == 'stagger2_hd_patch' ) then
          patch%mem(i1,i2,i3,patch%idx%s,patch%it,1) = (log(P0) - log(sigma0) * patch%gamma) &
                                                        * sigma0 / (patch%gamma - 1.0)
        else if (trim(patch%kind) == 'stagger2e_hd_patch') then
          patch%mem(i1,i2,i3,patch%idx%e,patch%it,1) = P0 / (patch%gamma - 1.0)
        else if (trim(patch%kind) == 'ramses_hd_patch'   ) then
          patch%mem(i1,i2,i3,patch%idx%e,patch%it,1) = P0 / (patch%gamma - 1.0) &
                                                        + 0.5*sigma0*vkep(i2)**2
        else if (trim(patch%kind) == 'zeus_hd_patch'     ) then
          patch%mem(i1,i2,i3,patch%idx%e,patch%it,1) = P0 / (patch%gamma - 1.0)
        end if
        if (patch%nc == 2) then
          dd (i1,i2,i3) = sigma0 * experiment_params%dust2gas
          patch%mem(i1,i2,i3,patch%idx%s+5,patch%it,1) = 0.0
          pxd(i1,i2,i3) = 0.0
          pyd(i1,i2,i3) = 0.0
          pzd(i1,i2,i3) = vkep(i2) * sigma0 * experiment_params%dust2gas!* (1.0 - experiment_params%aspect_rat**2)
          pyd(i1,i2,i3) = pyd(i1,i2,i3) * m1%h2c(i1)
          pzd(i1,i2,i3) = pzd(i1,i2,i3) * m1%h31c(i1) * m2%h32c(i2)
        end if
      end do
    end do
  end do
  deallocate (vkep)
  nullify (m1, m2, m3)
  call trace_end
END SUBROUTINE initial_conditions

!> =============================================================================
!> Download guard zone values
!> =============================================================================
SUBROUTINE dnload (self, only)
  class(experiment_t):: self
  integer, optional:: only
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%dnload')
  call download%download_link (self%link, .true., only=only)
  call trace%end
END SUBROUTINE dnload

END MODULE experiment_mod
