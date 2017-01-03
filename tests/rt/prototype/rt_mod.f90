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
  USE patch_mod
  USE link_mod
  USE rt_bundle_mod
  USE eos_mod
  implicit none
  private
  ! --------------------------------------------------------------------
  type, public, extends(patch_t):: rt_t
    integer                           :: dim(3)                         ! patch dimensions
    integer                           :: n_lambda                       ! number of wavelengths
    integer                           :: n_bundle                       ! number of bundles
    real, dimension(:), pointer       :: w_lambda                       ! wavelength weights
    real, dimension(:), pointer       :: w_bundle                       ! space angle weights
    type(rt_bundle_t)                 :: bundle
    real, dimension(:,:,:),   pointer :: dtau
    real, dimension(:,:,:,:), pointer :: rk
    real, dimension(:,:,:,:), pointer :: src
    real, dimension(:,:,:),   pointer :: dedt
    real, dimension(:,:),     pointer :: qfw0
    real, dimension(:,:),     pointer :: qrv1
    type (link_t),            pointer :: upstream
    type (link_t),            pointer :: dnstream
  contains
    procedure :: init
    procedure :: update
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
  self%nv = 1
  self%nw = 1
  call self%patch_t%init
  self%i_e = 1
  self%dedt => self%mem(:,:,:,self%i_e,self%it,1)
  rt_bundle%dim = self%dim
  allocate (self%w_lambda(self%n_lambda))
  allocate (self%w_bundle(self%n_bundle))
  allocate (self%rk (self%dim(1), self%dim(2), self%dim(3), self%n_lambda))
  allocate (self%src(self%dim(1), self%dim(2), self%dim(3), self%n_lambda))
  self%grace = 1.0                                                      ! use values for one dtime
  self%dtime = 0.001                            ! FIXME
  call rt_bundle%init
  call trace_end
END SUBROUTINE init

!> =====================================================================
!> Solve the RT problem, by first obtaining the thermodynamics at the
!> designated RT time, and using them to compute the opacity and source
!> function values.  Then solve for the net heating rate, storing it in
!> a slot that rotates with time.
!> =====================================================================
SUBROUTINE update (self)
  class (rt_t) :: self
  real, dimension(:,:,:), pointer:: lnd, lne
  class(task_t),  pointer:: task
  class(patch_t), pointer:: mhd
  integer, save:: itimer=0
  !---------------------------------------------------------------------
  call trace%begin ('rmhd_t%update', itimer=itimer)
  if (self%time < 0d0) then
    self%dtime = 0.01                                                   ! 1 second
    self%time = self%dtime
    self%grace = 0.0                                                    ! no grace
  end if
  task => self%link%nbor%task                                           ! assumes the only nbor is RMHD
  select type (task)
  class is (patch_t)
  mhd => task
  allocate (lnd (mhd%gn(1), mhd%gn(2), mhd%gn(3)))                      ! log density scratch
  allocate (lne (mhd%gn(1), mhd%gn(2), mhd%gn(3)))                      ! log internal energy scratch
  call mhd%log_interpolate (self%time, mhd%i_d, lnd)                    ! interpolate in time and space
  call mhd%log_interpolate (self%time, mhd%i_e, lne)                    ! interpolate in time and space
  end select
  lne = lne - lnd                                                       ! internal energy per unit mass
  call eos%lookup (shape(lne), lne=lne, lnd=lnd, &                      ! call the table lookup, ...
                   rk=self%rk, src=self%src)                            ! ... results go into rt_t
  deallocate (lnd, lne)                                                 ! deallocate scratch
  call trace%end (itimer)
  call self%solve                                                       ! call the RT solver to get Q
  !-----------------------------------------------------------------------------
  ! Update the RT time.  When first starting, the RT time is negative, to force
  ! a first RT solution, which is then assigned time=dtime, which is the time
  ! interval over which the RT solution is deemed valid (with grace=0.0). The
  ! next RT solution is triggered when the MHD time exceeds the RT time, and
  ! after the RT is done the RT time is advanced with the RT time step dtime.
  !-----------------------------------------------------------------------------
END SUBROUTINE update

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
  real, dimension(:,:),   pointer :: qfw0, qfw1, qrv0, qrv1
  real, dimension(:,:,:), pointer :: s, rk, dtau, exp0, q
  class (patch_t),        pointer :: task
  integer, save                   :: itimer=0
  ! --------------------------------------------------------------------
  call trace_begin ('rt_t%solve', itimer=itimer)
  self%dedt = 0.0                                                       ! reset net cooling
  do ib=1,self%n_bundle                                                 ! loop over bundles
    !$omp single
    do iw=1,self%n_lambda                                               ! loop over wavelength
      !$omp task private(iw,s,rk,dtau,exp0)                             ! spawn forward task
      s  => self%src(:,:,:,iw)                                          ! pick source function
      rk => self%rk (:,:,:,iw)                                          ! pick opacity
      allocate (dtau(self%dim(1), self%dim(2), self%dim(3)))            ! temporary array
      allocate (exp0(self%dim(1), self%dim(2), self%dim(3)))            ! temporary array
      allocate (q   (self%dim(1), self%dim(2), self%dim(3)))            ! temporary array
      allocate (qfw0(self%dim(1), self%dim(3)))                         ! temporary array
      allocate (qfw1(self%dim(1), self%dim(3)))                         ! temporary array
      allocate (qrv0(self%dim(1), self%dim(3)))                         ! temporary array
      allocate (qrv1(self%dim(1), self%dim(3)))                         ! temporary array
      ds => self%mesh(2)%dsup                                           ! along this direction
      call self%optical_depth (ds, rk, dtau)                            ! compute dtau
      !
      exp0 = exp(-dtau)                                                 ! expopential factors
      qfw0 = -s(:,5,:)                                                  ! input net intensity
      call rt_bundle%forward (dtau, exp0, s, qfw0, q, qfw1)             ! solve forward
      !task => self%dnstream%task
      !select type (task)
      !class is (rt_t)
      !  task%qfw0 = qfw1                                                ! output net intensity
      !end select
      self%dedt = self%dedt + q*self%w_lambda(iw)*self%w_bundle(ib)     ! add net heat contrib
      !
      qrv1 = 0.0                                                        ! input net intensity
      call rt_bundle%reverse (dtau, exp0, s, qrv1, q, qrv0)             ! solve reverse
      !task => self%upstream%task
      !select type (task)
      !class is (rt_t)
      !  task%qrv1 = qrv0                                                ! output net intensity
      !end select
      self%dedt = self%dedt + q*self%w_lambda(iw)*self%w_bundle(ib)     ! add net heat contrib
      !
      deallocate (dtau, exp0, q)                                        ! deallocate temporaries
      !$omp end task
    end do
    !$omp end single
    !$omp taskwait
  end do
  call trace_end (itimer)
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
