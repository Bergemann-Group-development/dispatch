!=======================================================================
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
!> the exponential factor only neeas to be computed once, and can be
!> used in both solutions.
!=======================================================================
MODULE rt_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  USE rt_bundle_mod
  USE eos_mod
  USE bits_mod
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
    real, dimension(:,:,:),   pointer :: dtau => null()
    real, dimension(:,:,:,:), pointer :: rk   => null()
    real, dimension(:,:,:,:), pointer :: src  => null()
    real, dimension(:,:,:),   pointer :: dedt => null()
    real, dimension(:,:),     pointer :: qfw0 => null()
    real, dimension(:,:),     pointer :: qrv1 => null()
    real, dimension(:,:,:), pointer   :: lnds => null()                 ! scratch pointer to lnd
    real, dimension(:,:,:), pointer   :: lnes => null()                 ! scratch pointer to lne

    character(len=16)                 :: solver ='integral'             ! choice of solver
    type (link_t),            pointer :: upstream
    type (link_t),            pointer :: dnstream
    procedure(solve_integral), pointer:: solve => null()
  contains
    procedure :: init
    procedure :: dnload
    procedure :: update
    procedure :: solve_integral
    procedure :: solve_feautrier
    procedure :: optical_depth_centered
    procedure :: optical_depth_staggered
    procedure :: density_weight_factor
    procedure :: interpolate_rt
    procedure :: boundary_conditions
  end type
  real, parameter:: pi4 = 16.*atan(1.)
CONTAINS

!=======================================================================
!> We may or may not need locally allocated arrays.  If the calling
!> procedure holds the necessary data (which it does, since it computes
!> them and is picking up the output), it just needs to set the pointers
!> in this module to point to the existing arrays.  However, the optical
!> depth increments and the corresponding exp(-dtau) factors are
!> computed here, so they need array allocation.
!=======================================================================
SUBROUTINE init (self)
  class (rt_t) :: self
  integer:: ignore
  logical, save:: first_time=.true.
  real(8), save:: dtime=0.01_8, grace=0.0_8
  character(len=16) :: solver = 'integral'
  namelist /rt_params/ dtime, grace, solver
  ! --------------------------------------------------------------------
  call trace_begin ('rt_t%init')
  call io%lock%set
  if (first_time) then
    rewind (io%input); read (io%input, rt_params, iostat=ignore)
    if (io%master) write (*, rt_params)
  end if
  call io%lock%unset
  self%nv = 1                                                           ! time slots for dedt
  self%nw = 1                                                           ! no time deriv
  call self%patch_t%init                                                ! initialize patch
  call self%set (bits%static+bits%no_download)                          ! special bits
  self%idx%e = 1                                                          ! dedt slot
  self%n_lambda = eos%n_lambda                                          ! number of wavelengths
  self%n_bundle = 1                          ! FIXME                    ! one direction
  self%dedt => self%mem(:,:,:,self%idx%e,self%it,1)                       ! map dedt to mem
  rt_bundle%dim = self%dim                                              ! pass on dims
  allocate (self%w_lambda(self%n_lambda))
  allocate (self%w_bundle(self%n_bundle))
  allocate (self%rk  (self%dim(1), self%dim(2), self%dim(3), self%n_lambda))
  allocate (self%src (self%dim(1), self%dim(2), self%dim(3), self%n_lambda))
  self%w_lambda = eos%w_lambda
  self%w_bundle = pi4/2./self%n_bundle          ! FIXME                 ! space angle weights
  self%dtime = dtime                                                    ! 0.1 second default
  self%grace = grace                                                    ! grace off by default
  call rt_bundle%init
  self%solver = solver
  select case(self%solver)
  case('integral')
    self%solve => solve_integral
  case('feautrier')
    self%solve => solve_feautrier
  case default
    print *, 'WARNING: rt solver not recognised. Using integral solver.'
    self%solve => solve_integral
  end select
  call trace_end
END SUBROUTINE init

!=======================================================================
!=======================================================================
SUBROUTINE dnload (self)
  class (rt_t) :: self
END SUBROUTINE dnload

!=======================================================================
!> Solve the RT problem, by first obtaining the thermodynamics at the
!> designated RT time, and using them to compute the opacity and source
!> function values.  Then solve for the net heating rate, storing it in
!> a slot that rotates with time.
!=======================================================================
SUBROUTINE update (self)
  class (rt_t) :: self
  real, dimension(:,:,:), pointer:: lnd, lne, ne, tt, pg
  class(task_t),  pointer:: task
  class(patch_t), pointer:: mhd
  integer, save:: itimer=0
  integer :: ii, iy
  !---------------------------------------------------------------------
  call trace%begin ('rmhd_t%update', itimer=itimer)
  if (io%verbose > 1) print *, 'RADIATIVE TRANSFER TIME', self%time
  if (self%time < 0d0) then
    self%time  = 0.0                                                    ! set current time
  end if
  task => self%link%nbor%task                                           ! assumes the only nbor is RMHD
  select type (task)
  class is (patch_t)
    mhd => task
  end select
  allocate (lnd (mhd%gn(1), mhd%gn(2), mhd%gn(3)))                      ! log density scratch
 ! allocate(ne(mhd%gn(1), mhd%gn(2), mhd%gn(3)))
 ! allocate(tt(mhd%gn(1), mhd%gn(2), mhd%gn(3)))
 ! allocate(pg(mhd%gn(1), mhd%gn(2), mhd%gn(3)))
  allocate (lne (mhd%gn(1), mhd%gn(2), mhd%gn(3)))                      ! log internal energy scratch
  self%lnds => lnd
  self%lnes => lne
  call mhd%log_interpolate (self%time, mhd%idx%d, lnd)                    ! interpolate in time and space
  call mhd%log_interpolate (self%time, mhd%idx%e, lne)                    ! interpolate in time and space
  if (io%verbose > 1) then
    print *, self%id, 'lnd', minval(lnd), maxval(lnd)
    print *, self%id, 'lne', minval(lne), maxval(lne)
  end if
  lne = lne - lnd                                                       ! internal energy per unit mass
  
  call eos%lookup (shape(lne), lne=lne, lnd=lnd, &                      ! call the table lookup, ...
                   rk=self%rk, src=self%src)                            ! ... results go into rt_t      
  call trace%end (itimer)
  call self%solve                                                       ! call the RT solver to get Q
  call self%boundary_conditions
  call self%output
  deallocate (lnd, lne)                                                 ! deallocate scratch
END SUBROUTINE update

!=======================================================================
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
!=======================================================================
SUBROUTINE solve_integral (self)
  class (rt_t) :: self
  !.....................................................................
  integer                         :: ib, iw, iy
  real, dimension(:),     pointer :: ds
  real, dimension(:,:),   pointer :: qfw0, qfw1, qrv0, qrv1
  real, dimension(:,:,:), pointer :: s, rk, dtau, exp0, q
  class (patch_t),        pointer :: task
  integer, save                   :: itimer=0
  ! --------------------------------------------------------------------
  call trace_begin ('rt_t%solve', itimer=itimer)
  self%dedt = 0.0                                                       ! reset net cooling
  associate (mx=>self%dim(1), my=>self%dim(2), mz=>self%dim(3))
  do ib=1,self%n_bundle                                                 ! loop over bundles
    !$omp single
    do iw=1,self%n_lambda                                               ! loop over wavelength
      !$omp task private(s,rk,dtau,exp0)                                ! spawn forward task
      s  => self%src(:,6:my-4,:,iw)                                     ! pick source function
      rk => self%rk (:,6:my-4,:,iw)                                     ! pick opacity
      allocate (dtau(mx,1:my-9,mz))                                     ! temporary array
      allocate (exp0(mx,1:my-9,mz))                                     ! temporary array
      allocate (q   (mx,1:my-9,mz))                                     ! temporary array
      allocate (qfw0(mx,mz))                                            ! temporary array
      allocate (qfw1(mx,mz))                                            ! temporary array
      allocate (qrv0(mx,mz))                                            ! temporary array
      allocate (qrv1(mx,mz))                                            ! temporary array
      ds => self%mesh(2)%dsup(6:my-4)                                   ! along this direction
      call self%optical_depth_centered (ds, rk, dtau)                   ! compute dtau
      if (io%verbose>1) &
        print *,'dtau',iw,self%id, minval(dtau), maxval(dtau)
      !
      
#if defined __INTEL_COMPILER
      exp0 = 2**(-dtau*1.44269504) 
#else
      exp0 = exp(-dtau)                                                 ! expopential factors
#endif
    
      qfw0 = -s(:,6,:)                                                  ! input net intensity
      call rt_bundle%forward (dtau, exp0, s, qfw0, q, qfw1)             ! solve forward
      !task => self%dnstream%task
      !select type (task)
      !class is (rt_t)
      !  task%qfw0 = qfw1                                                ! output net intensity
      !end select
      self%dedt(:,6:my-4,:) = &
      self%dedt(:,6:my-4,:) + q*self%w_lambda(iw)*self%w_bundle(ib)     ! add net heat contrib
      !
      qrv1 = 0.0                                                        ! input net intensity
      call rt_bundle%reverse (dtau, exp0, s, qrv1, q, qrv0)             ! solve reverse
      !task => self%upstream%task
      !select type (task)
      !class is (rt_t)
      !  task%qrv1 = qrv0                                                ! output net intensity
      !end select
      self%dedt(:,6:my-4,:) = &
      self%dedt(:,6:my-4,:) + q*self%w_lambda(iw)*self%w_bundle(ib)     ! add net heat contrib
      if ((io%verbose>1 .and. iw==self%n_lambda) .or. io%verbose>2) then
        do iy=1,my-9
          print *,iy,sum(dtau(:,iy,:))/(mx*mz),sum(self%dedt(:,5+iy,:))/(mx*mz)
        end do
      end if
      !
      deallocate (dtau, exp0, q, qfw0, qfw1, qrv0, qrv1)                ! deallocate temporaries
      !$omp end task
    end do
    !$omp end single
    !$omp taskwait
  end do
  end associate
  call trace_end (itimer)
END SUBROUTINE solve_integral

SUBROUTINE solve_feautrier (self)
  class (rt_t) :: self
  !.....................................................................
  integer                         :: ib, iw, iy
  real, dimension(:),     pointer :: ds
  real, dimension(:,:),   pointer :: Ifw0, Ifw1, Irv0, Irv1
  real, dimension(:,:,:), pointer :: s, rk, dtau, exp0, q, dwf
  class (patch_t),        pointer :: task
  integer, save                   :: itimer=0
  ! --------------------------------------------------------------------
  call trace_begin ('rt_t%solve', itimer=itimer)
  self%dedt = 0.0                                                       ! reset net cooling
  associate (mx=>self%dim(1), my=>self%dim(2), mz=>self%dim(3))
  do ib=1,self%n_bundle                                                 ! loop over bundles
    !$omp single
    do iw=1,self%n_lambda                                               ! loop over wavelength
      s  => self%src(:,6:my-4,:,iw)                                     ! pick source function
      rk => self%rk (:,6:my-4,:,iw)                                     ! pick opacity
      allocate (dtau(mx,1:my-9,mz))                                     ! temporary array
      allocate (q   (mx,1:my-9,mz))                                     ! temporary array
      allocate (Ifw0(mx,mz))                                            ! temporary array
      allocate (Ifw1(mx,mz))                                            ! temporary array
      allocate (Irv0(mx,mz))                                            ! temporary array
      allocate (Irv1(mx,mz))                                            ! temporary array
      allocate (dwf (mx,1:my-9,mz))                                     ! temporary array
      ds => self%mesh(2)%dsup(6:my-4)                                   ! along this direction
      call self%optical_depth_staggered (ds, rk, dtau)                  ! compute dtau
      call self%density_weight_factor(ds,dwf)
      if (io%verbose>1) &
        print *,'dtau',iw,self%id, minval(dtau), maxval(dtau)
          
      Ifw0 = -s(:,6,:)                                                  ! input net intensity
      Irv1 = 0.0
      
      call rt_bundle%transfer (dtau, s, Ifw0, Ifw1, q, Irv1, Irv0)      ! solve bi-directional

      self%dedt(:,6:my-4,:) = &
      self%dedt(:,6:my-4,:) + q*self%w_lambda(iw)*self%w_bundle(ib) * dwf ! add net heat contrib

      deallocate (dtau, q, Ifw0, Ifw1, Irv0, Irv1, dwf)                ! deallocate temporaries
    end do
    !$omp end single
  end do
  end associate
  call trace_end (itimer)
END SUBROUTINE solve_feautrier

!=======================================================================
!> Compute the optical depth along rays
!=======================================================================
SUBROUTINE optical_depth_centered (self, ds, rk, dtau)
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
END SUBROUTINE optical_depth_centered

!=======================================================================
!> Compute the optical depth along rays. Dtau is staggered, therefore it
!> has only a half value on the lower edge of the patch.
!=======================================================================
SUBROUTINE density_weight_factor (self, ds, dwf)
  class (rt_t) :: self
  real, dimension(:),     pointer :: ds
  real, dimension(:,:,:), pointer :: dwf
  !.....................................................................
  integer                         :: iy, iz
  !---------------------------------------------------------------------
  do iz=1,size(dwf,3)
    do iy=1,size(dwf,2)
      dwf(:,iy,iz) = self%lnds(:,iy,iz) * ds(iy)
    end do
  end do
END SUBROUTINE density_weight_factor

!=======================================================================
!> Compute the optical depth along rays. Dtau is staggered, therefore it
!> has only a half value on the lower edge of the patch.
!=======================================================================
SUBROUTINE optical_depth_staggered (self, ds, rk, dtau)
  class (rt_t) :: self
  real, dimension(:),     pointer :: ds
  real, dimension(:,:,:), pointer :: rk
  real, dimension(:,:,:), pointer :: dtau,dt
  !.....................................................................
  integer                         :: iy, iz
  !---------------------------------------------------------------------
  allocate (dt(size(dtau,1),size(dtau,2),size(dtau,3)))
  do iz=1,size(rk,3)
    do iy=1,size(rk,2)
      dtau(:,iy,iz) = rk(:,iy,iz)* ds(iy)
    end do
    dt=dtau
    do iy=2,size(rk,2)-1
      dtau(:,iy,iz) = (dtau(:,iy,iz) + dt(:,iy-1,iz)) * .5
    end do
  end do
  deallocate(dt)
END SUBROUTINE optical_depth_staggered

!=======================================================================
!> Compute the optical depth along rays
!=======================================================================
SUBROUTINE interpolate_rt (self, time)
  class (rt_t) :: self
  real(8):: time
  !.....................................................................
  !---------------------------------------------------------------------
END SUBROUTINE interpolate_rt

!===============================================================================
!> Apply boundary conditions
!===============================================================================
SUBROUTINE boundary_conditions (self)
  class(rt_t):: self
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%boundary_conditions')
  associate (my => self%mesh(2)%gn)
  self%mem(:,  01:06,:,:,:,2) = 0.0
  self%mem(:,my-4:my,:,:,:,2) = 0.0
  end associate
  call trace%end()
END SUBROUTINE boundary_conditions

END MODULE rt_mod
