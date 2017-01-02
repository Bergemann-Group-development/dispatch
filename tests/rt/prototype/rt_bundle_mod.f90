!> =====================================================================
!> The input to this module are the opacities and source functions at
!> the points belonging to the bundle, and the net radiative intensity
!> self%qin at the upstream face of the bundle.
!> 
!> The module task (solve) is to compute the corresponding net cooling
!> self%q throughout the volume, and also to compute the net radiative
!> intensity on the downstream bundle face.
!>
!> This module knows nothing about context, including boundaries.  If
!> the patch is a (virtual or boundary) patch, the calling procedure
!> has picked up the incoming values from the boundary conditions, and
!> otherwise it has them from the upstream patch.
!>
!> Rather than copy the input values to scratch arrays, the calling
!> procedure just points to the values, which are part of its data.
!> =====================================================================
MODULE rt_bundle_mod
  USE io_mod
  USE trace_mod
  USE rt_solver_mod
  implicit none
  private
  ! --------------------------------------------------------------------
  type, public:: rt_bundle_t
    integer                         :: dim(3)                           ! patch dimensions
    real, dimension(:,:,:), pointer :: q, s, dtau
    real, dimension(:,:),   pointer :: qin
    real, dimension(:,:),   pointer :: qout
  contains
    procedure :: init
    procedure :: solve
    procedure :: forward
    procedure :: reverse
  end type
  type(rt_bundle_t), public:: rt_bundle
! ----------------------------------------------------------------------
CONTAINS

!> =====================================================================
!> Initialize the solver, which chooses the solver type
!> =====================================================================
SUBROUTINE init (self)
  class (rt_bundle_t) :: self
  ! --------------------------------------------------------------------
  call trace%begin ('rt_bundle_t%init')
  call rt_solver%init
  call trace%end
END SUBROUTINE init

!> =====================================================================
!> Solve the radiative transfer equation with rt_solver%transfer
!> =====================================================================
SUBROUTINE solve (self)
  class (rt_bundle_t) :: self
  real, dimension(:),   pointer :: qfw0, qfw1
  real, dimension(:,:), pointer :: q, s, dtau
  integer :: iz
  ! --------------------------------------------------------------------
  call trace%begin ('rt_bundle_t%update')
  do iz=1,size(self%s,3)
    qfw0 => self%qin (:,iz)             ! input surface values
    s    => self%s   (:,:,iz)           ! source function
    dtau => self%dtau(:,:,iz)           ! optical depth increments
    qfw1 => self%qout(:,iz)             ! net volume cooling
    q    => self%q   (:,:,iz)           ! output surface values
    call rt_solver%transfer (dtau, S, qfw0, Q, qfw1)
  end do
  call trace%end
END SUBROUTINE solve

!> =====================================================================
!> Solve the radiative transfer equation with rt_solver%transfer
!> =====================================================================
SUBROUTINE forward (self)
  class (rt_bundle_t) :: self
  real, dimension(:),   pointer :: qfw0, qfw1
  real, dimension(:,:), pointer :: q, s, dtau
  integer :: iz
  ! --------------------------------------------------------------------
  call trace%begin ('rt_bundle_t%forward')
  do iz=1,size(self%s,3)
    qfw0 => self%qin (:,iz)             ! input surface values
    s    => self%s   (:,:,iz)           ! source function
    dtau => self%dtau(:,:,iz)           ! optical depth increments
    qfw1 => self%qout(:,iz)             ! net volume cooling
    q    => self%q   (:,:,iz)           ! output surface values
    call rt_solver%forward (dtau, S, qfw0, Q, qfw1)
  end do
  call trace%end
END SUBROUTINE forward

!> =====================================================================
!> Solve the radiative transfer equation with rt_solver%transfer
!> =====================================================================
SUBROUTINE reverse (self)
  class (rt_bundle_t) :: self
  real, dimension(:),   pointer :: qfw0, qfw1
  real, dimension(:,:), pointer :: q, s, dtau
  integer :: iz
  ! --------------------------------------------------------------------
  call trace%begin ('rt_bundle_t%reverse')
  do iz=1,size(self%s,3)
    qfw0 => self%qin (:,iz)             ! input surface values
    s    => self%s   (:,:,iz)           ! source function
    dtau => self%dtau(:,:,iz)           ! optical depth increments
    qfw1 => self%qout(:,iz)             ! net volume cooling
    q    => self%q   (:,:,iz)           ! output surface values
    call rt_solver%reverse (dtau, S, qfw0, Q, qfw1)
  end do
  call trace%end
END SUBROUTINE reverse

END MODULE rt_bundle_mod
