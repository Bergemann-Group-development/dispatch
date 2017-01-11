!> =============================================================================
!> Initializiation for the shocked cylinder experiment
!> =============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE scaling_mod
  USE mesh_mod
  USE mhd_mod
  USE parameters_mod
  implicit none
  private

  real(8), parameter :: pi = 3.14159265358979323846d0

  public :: initial_conditions
CONTAINS

!> =============================================================================
!> Initialise the MHD variables
!> =============================================================================
SUBROUTINE initial_conditions (patch)
  class(mhd_t), target :: patch
  class(mesh_t), pointer :: m1, m2, m3
  integer :: i1, i2, i3
  real(4) :: z, r, phi, zeta_offset
  real(8) :: sigma0
  real(8), dimension(:), allocatable :: vkep
  !----------------------------------------------------------------------------
  call trace_begin('initial_mod::initial_conditions')
  
  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)

  ! calculate the Keplerian velocity
  allocate(vkep(m2%gn-1))
  vkep(:) = sqrt(experiment_params%mstar * scaling%grav / (m2%centre_nat + patch%mesh(2)%r(:)))
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
        patch%d%q(i1,i2,i3) = sigma0
        patch%e1%q(i1,i2,i3) = (experiment_params%aspect_rat * vkep(i2))**2 &
                             * sigma0 / (patch%gamma - 1.0)
        patch%s%x%q(i1,i2,i3) = 0.0
        patch%s%y%q(i1,i2,i3) = 0.0
        patch%s%z%q(i1,i2,i3) = vkep(i2) * sigma0
        patch%b%x%q(i1,i2,i3) = 0.0
        patch%b%y%q(i1,i2,i3) = 0.0
        patch%b%z%q(i1,i2,i3) = 0.0
        ! the ZEUS module stores the angular momentum rather than the linear momentum
        patch%s%y%q(i1,i2,i3) = patch%s%y%q(i1,i2,i3) * m1%h2c(i1)
        patch%s%z%q(i1,i2,i3) = patch%s%z%q(i1,i2,i3) * m1%h31c(i1) * m2%h32c(i2)
      end do
    end do
  end do

  deallocate (vkep)
  nullify (m1, m2, m3)
  call trace_end
END SUBROUTINE initial_conditions

END MODULE initial_mod
