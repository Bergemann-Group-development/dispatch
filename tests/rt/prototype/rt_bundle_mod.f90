!=======================================================================
!> The inputs to this module are the opacities and source functions at
!> the points belonging to the bundle, and the net radiative intensity
!> self%qin at the upstream face of the bundle.
!>
!> The module task is to compute the corresponding net cooling
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
!=======================================================================
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
    procedure, nopass :: forward
    procedure, nopass :: reverse
    procedure, nopass :: transfer
  end type
  type(rt_bundle_t), public:: rt_bundle
! ----------------------------------------------------------------------
CONTAINS

!=======================================================================
!> Initialize the solver, which chooses the solver type
!=======================================================================
SUBROUTINE init (self)
  implicit none
  class (rt_bundle_t) :: self
  ! --------------------------------------------------------------------
  call trace%begin ('rt_bundle_t%init')
  call rt_solver%init
  call trace%end
END SUBROUTINE init

!=======================================================================
!> Solve the radiative transfer equation with rt_solver%transfer
!=======================================================================
SUBROUTINE forward (dtau, exp0, s, qin, q, qout)
  implicit none
  real, dimension(:,:,:), pointer :: dtau, exp0, s, q
  real, dimension(:,:),   pointer :: qin, qout
  !.....................................................................
  real, dimension(:),   pointer   :: qfw0, qfw1
  real, dimension(:,:), pointer   :: q1, s1, dtau1, ex0
  integer :: iz
  ! --------------------------------------------------------------------
  call trace%begin ('rt_bundle_t%forward')
  do iz=1,size(s,3)
    qfw0  => qin (:,iz)                 ! input surface values
    s1    => s   (:,:,iz)               ! source function
    dtau1 => dtau(:,:,iz)               ! optical depth increments
    ex0   => exp0(:,:,iz)               ! exp(-dtau)
    qfw1  => qout(:,iz)                 ! net volume cooling
    q1    => q   (:,:,iz)               ! output surface values
    call rt_solver%forward (dtau1, ex0, s1, qfw0, q1, qfw1)
  end do
  call trace%end
END SUBROUTINE forward

!=======================================================================
!> Solve the radiative transfer equation with rt_solver%transfer
!=======================================================================
SUBROUTINE reverse (dtau, exp0, s, qin, q, qout)
  implicit none
  real, dimension(:,:,:), pointer :: dtau, exp0, s, q
  real, dimension(:,:),   pointer :: qin, qout
  !.....................................................................
  real, dimension(:),     pointer :: qrv0, qrv1
  real, dimension(:,:),   pointer :: q1, s1, dtau1, ex0
  integer :: iz
  ! --------------------------------------------------------------------
  call trace%begin ('rt_bundle_t%reverse')
  do iz=1,size(s,3)
    qrv1  => qin (:,iz)                 ! input surface values
    s1    => s   (:,:,iz)               ! source function
    dtau1 => dtau(:,:,iz)               ! optical depth increments
    ex0   => exp0(:,:,iz)               ! exp(-dtau)
    qrv0  => qout(:,iz)                 ! net volume cooling
    q1    => q   (:,:,iz)               ! output surface values
    call rt_solver%reverse (dtau1, ex0, s1, qrv1, q1, qrv0)
  end do
  call trace%end
END SUBROUTINE reverse

!=======================================================================
!> Solve the radiative transfer equation with rt_solver%transfer
!=======================================================================
SUBROUTINE transfer (dtau, s, Ifw0, Ifw1, q, Irv1, Irv0)
  implicit none
  real, dimension(:,:,:), pointer :: dtau, s, q
  real, dimension(:,:),   pointer :: Ifw0, Ifw1, Irv1, Irv0
  !.....................................................................
  real, dimension(:),     pointer :: xi
  real, dimension(:,:),   pointer :: q1, s1, dtau1, ex0,q12,q13,q14
  integer :: iz, a, b
  ! --------------------------------------------------------------------
  call trace%begin ('rt_bundle_t%reverse')
  allocate(xi(size(Ifw1,1)))

  a=size(s,1)
  b=size(s,2)
  do iz=1,size(s,3)
    s1    => s   (:,:,iz)               ! source function
    dtau1 => dtau(:,:,iz)               ! optical depth increments
    q1    => q   (:,:,iz)               ! output surface values
                 
    call rt_solver%feautrier_stag(dtau=dtau1, s=s1,q=q1, xi=xi, doxi=.false.)
           
  end do

  deallocate(xi)
  call trace%end
END SUBROUTINE transfer


END MODULE rt_bundle_mod
