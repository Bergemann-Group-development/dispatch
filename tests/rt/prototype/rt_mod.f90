!> =====================================================================
!> The input to this module is the source function and opacity values
!> at all grid points, which are computed by the calling procedure, in
!> the rmhd_mod module.  The current module therefore does not need
!> to know anything about the mhd_mod values, nor does it have to know
!> anything about tables or equations of state.
!>
!> The task of the module is to compute the source function and optical
!> depth increments along all ray bundles that pass through the patch,
!> and add them to the net cooling rate, which is the output of the 
!> module.  In doing so, the module needs to pick up upstream or boundary
!> values of the net intensity q=I-s, and it needs to map the Cartesian
!> input values to values along the rays.
!>
!> In this simple prototype we assume there are only two bundles; both
!> vertical, and in opposite directions.
!> 
!> The module is allowed to generate two tasks for each direction, so
!> the exponetial factor only need to be computed once, and can be ued
!> in both solutions.
!> =====================================================================
MODULE rt_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE rt_bundle_mod
  implicit none
  private
  ! --------------------------------------------------------------------
  type, public, extends(task_t):: rt_t
    integer                           :: dim(3)                         ! patch dimensions
    real, dimension(:),       pointer :: ds                             ! step size along rays
    integer                           :: n_lambda                       ! number of wavelengths
    integer                           :: n_bundle                       ! number of bundles
    real, dimension(:), pointer       :: w_lambda                       ! wavelength weights
    real, dimension(:), pointer       :: w_bundle                       ! space angle weights
    type(rt_bundle_t)                 :: bundle
    real, dimension(:,:,:),   pointer :: rk
    real, dimension(:,:,:),   pointer :: dtau
    real, dimension(:,:,:,:), pointer :: src
    real, dimension(:,:,:),   pointer :: dedt
  contains
    procedure :: init
    procedure :: solve
    procedure :: optical_depth
    procedure :: interpolate
  end type
! ----------------------------------------------------------------------
CONTAINS

!> =====================================================================
!> We may or may not need locally allocated arrays.  If the calling 
!> procedure holds the necessary data (which it does, since it computes
!> them and is picking up the output), it just needs to set the pointers
!> in this module to point to the existing arrays.  However, the optical
!> depth increments and the corresponding exp(-dtau) factors are 
!> computed here, so they need array allocation.
!> =====================================================================
SUBROUTINE init (self)
  class (rt_t) :: self
  ! --------------------------------------------------------------------
  call trace_begin ('rt_t%init')
  rt_bundle%dim = self%dim
  allocate (self%w_lambda(self%n_lambda))
  allocate (self%w_bundle(self%n_bundle))
  allocate (self%dedt (self%dim(1), self%dim(2), self%dim(3)))
  self%grace = 1.0                                                      ! use values for one dtime
  self%dtime = 0.001                            ! FIXME
  call rt_bundle%init
  call trace_end
END SUBROUTINE init

!> =====================================================================
!> This procedure should loop over all bundle directions, spawning two
!> tasks for each direction, with solvers running in the forward and
!> reverse directions.
!>
!> The procedure implies an implicit wait for the upstream net intensity
!> which arrives at different times for the two directions.  The first
!> call increments the counter from 0 to 1, and when the 2nd call makes
!> it 2, the procedure knows that temporary work arrays should be 
!> deallocated.
!>
!> In this prototype case we use only one (vertical = y) direction, but
!> in the general case we loop over the included ray directions, doing
!> a waittask after the loop, so all tasks are active at the same time.
!> Here we spawn one task per wavelength, just to demonstrate tasks.
!>
!> The tasks that are missing upstream values need to to a taskwait,
!> and try again later.  However, in the current special case we always
!> have both boundart values, so there are no waits.
!> =====================================================================
SUBROUTINE solve (self)
  class (rt_t) :: self
  !.....................................................................
  integer:: ib, iw
  real, dimension(:),     pointer :: ds
  real, dimension(:,:,:), pointer :: s, rk, dtau, exp1
  ! --------------------------------------------------------------------
  call trace_begin ('rt_t%solve')
  self%dedt = 0.0                                                       ! reset net cooling
  do ib=1,self%n_bundle                                                 ! loop over bundles
    do iw=1,self%n_lambda                                               ! loop over wavelength
      !$omp task private(iw,s,rk,dtau)                                  ! spawn forward task
      s  => self%src(:,:,:,iw)                                          ! pick source function
      rk => self%rk(:,:,:)                                              ! pick opacity
      allocate (dtau(self%dim(1), self%dim(2), self%dim(3)))            ! temporary array
      ds => self%ds                                                     ! along this direction
      call self%optical_depth (ds,rk,dtau)                              ! compute dtau
      exp1 = exp(-dtau)                                                 ! expopential factors
      call rt_bundle%forward                                            ! solve forward
      self%dedt = self%dedt + rt_bundle%q*self%w_lambda(iw)*self%w_bundle(ib)
      call rt_bundle%reverse                                            ! solve reverse
      self%dedt = self%dedt + rt_bundle%q*self%w_lambda(iw)*self%w_bundle(ib)
      deallocate (dtau)                                                 ! deallocate temporary
      !$omp end task
    end do
    !$omp taskwait
  end do
  call trace_end
END SUBROUTINE solve

!> =====================================================================
!> Compute the optical depth along rays
!> =====================================================================
SUBROUTINE optical_depth (self, ds, rk, dtau)
  class (rt_t) :: self
  real, dimension(:),     pointer :: ds
  real, dimension(:,:,:), pointer :: rk
  real, dimension(:,:,:), pointer :: dtau
  !.....................................................................
  integer                         :: iy, iz
  !---------------------------------------------------------------------
  do iz=1,size(rk,3)
  do iy=1,size(rk,2)
    dtau(:,iy,iz) = ds(iy)*rk(:,iy,iz)
  end do
  end do
END SUBROUTINE optical_depth

!> =====================================================================
!> Compute the optical depth along rays
!> =====================================================================
SUBROUTINE interpolate (self, time)
  class (rt_t) :: self
  real(8):: time
  !.....................................................................
  !---------------------------------------------------------------------
END SUBROUTINE interpolate

END MODULE rt_mod
