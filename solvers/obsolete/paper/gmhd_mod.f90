!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the Poisson problem in 3-D.  In a refined setting, where a smaller
!> patch exists inside a large one, we rely on the larger patch to provide
!> updates of the largest scales -- including at the boundaries, and we propagate
!> these changes to the smaller scale patch by adding it to the existing values.
!> The time and memory slots are used it the following way:
!>
!>                 level  initial guess of phi       iterative correction
!> time = t(now):   l     mem(:,:,:,i_phi,now,1)     mem(:,:,:,i_phi,now,2)
!>
!> The initial estimate for t(new) at the root grid is obtained from extrapolation
!> in time, using the previous time slot values.  Test this separately, first
!>
!> The initial value for the smaller patch cannot be obtained directly from just
!> time extrapolation, since it needs to match the boundary conditions set by
!> the larger patch, after it finishes iterating.  The best initial estimate
!> for smaller scale patches is thus obtained by adding the downloaded correction
!> from the larger scale to the previous solution.  Symbolically:
!>
!> mem(:,:,:,i_phi,prv,1) contains the previous solution, at time=t(prv)
!> mem(:,:,:,i_phi,now,2) is where we download the correction from lower levels to
!> mem(:,:,:,i_phi,now,1) = mem(:,:,:,i_phi,prv,1) + mem(:,:,:,i_phi,now,2) = IC
!> mem(:,:,:,i_phi,now,2) is finally filled with the difference, to be downloaded
!>
!> To support this from the larger patch requires that it too uses the same slots,
!> leaving when finished the correction on its scale in mem(:,:,:,i_phi,now,2).
!>
!> Steps        largest scale                   smaller scale
!> ----------   ---------------------------     ----------------------------------
!> initial:     mem1 = previous                 mem1 = previous
!>              mem2 = previous change          mem2 = previous local change
!> prediction:  mem1 = extrapolation            mem1 = extrapolation
!>              mem2 = change extrapolation     mem2 = change extrapolation
!> download:    mem1 = total extrapolation      mem1 = total extrapolation
!>              mem2 = none                     mem2 = downloaded large scale change
!> iteration:   mem1 = relaxation               mem1 = relaxation
!>              mem2 = new - old                mem2 = new-old
!>
!> The mem2 local correction is amenable to time extrapolation, which could
!> speed up convergence of the local iterations even more.  The constraints for
!> prediction on the smaller scale are: 1) it must not change the boundary
!> values (to be) imposed by the larger scale through addition, and 2) it must
!> NOT already contain (part of) the larger scale change to be added from the
!> larger scale.  It should thus entirely be referring to the local scale, which
!> is true if it is based entirely on the correction in step "iteration".  To be
!> clear, let's assume two scales that are entirely separate: A large (disk?)
!> scale, and a local (planetary) scale.  Changes on the small scale coming from
!> the large scale are slow, and are systematic throughout the small scale.  If
!> they are imposed only on the boundary the residual consist of a spectrum of
!> changes that includes the largest of the small scales (the whole small box).
!> If they are imposed as an additive term they are instead nearly exact, and
!> cause only minor extra work.
!>
!> The two accelerations should be tested separately.  First the addition of
!> the iterative change from the level above, and then the extrapolation in time
!> of the local change.
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE gmhd_mod
  USE io_mod
  USE omp_timer_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE scalar_mod
  USE stagger_mod
  USE laplace_mod
  USE mhd_mod
  USE boundaries_mod
  USE index_mod
  implicit none
  private
  type, public, extends(mhd_t):: gmhd_t
    real, dimension(:,:,:), pointer:: phi => null()
    real(kind=8):: G                                                                    ! Strength of gravity
    real(kind=8):: tolerance
    integer:: min_iter, max_iter, iter
    type(laplace_t)  :: laplace
  contains
    procedure:: init
    procedure:: update
    procedure:: solve
  end type
  !type(operator_t):: o
CONTAINS

!===============================================================================
!> Setup parameters for the Poisson solver
!===============================================================================
SUBROUTINE init (self)
  class(gmhd_t):: self
  real         :: tolerance=1e-4, G=1., omega_sor=0.0
  integer      :: min_iter=0, max_iter=100
  character(len=12):: solver='sor'
  logical:: first_time = .true.
  namelist /ghd_params/ G, tolerance, min_iter, max_iter, solver, omega_sor
  !-----------------------------------------------------------------------------
  ! Read standard parameters and copy to structure
  !-----------------------------------------------------------------------------
  call trace_begin ('gmhd_mod::init', 2)
  call self%mhd_t%init
  !call o%init (self%gn, self%ds)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, ghd_params)
    if (io%master) write (*,ghd_params)
  end if
  self%G         = G
  self%tolerance = tolerance
  self%min_iter  = min_iter
  self%max_iter  = max_iter
  self%nv        = 6
  !-----------------------------------------------------------------------------
  ! Set up the right dimensions from current patch, etc for laplace iteration method
  !-----------------------------------------------------------------------------
  call self%laplace%init (self, solver, omega_sor)
  !-----------------------------------------------------------------------------
  ! Memory mapping of density and potentiall field
  !-----------------------------------------------------------------------------
  self%d   => self%mem(:,:,:,self%idx%d  ,self%it,1)
  self%phi => self%mem(:,:,:,self%idx%phi,self%it,1)
  call trace_end
END SUBROUTINE init

!===============================================================================
!===============================================================================
!> Intercept the call from patch_list_mod::update to mhd_t%download, which
!> normally falls through to patch_t%update.  Instrument an extra flag, which
!> forces a download from lower levels to higher levels (only) for that mem
!> slot
!===============================================================================
!SUBROUTINE download_to (self, task, special, all_cells)
!  class(gmhd_t):: self
!  class(task_t), pointer:: task
!  integer:: special
!  logical, optional:: all_cells
!  !.............................................................................
!  special = 6
!  call self%task_t%download_to (task, special, all_cells)
!END SUBROUTINE download_to

!===============================================================================
!> Compute iteratively the solution to the Poisson equation, stopping when
!> reaching the tolerance tol
!===============================================================================
SUBROUTINE solve (self)
  class(gmhd_t):: self
  !
  real, pointer, dimension(:,:,:) :: FourPiGrho
  real, pointer, dimension(:,:,:) :: phi_1, phi_2
  real, pointer, dimension(:,:,:) :: phi_old, phi_new, phi_swap
  real(8), parameter :: FourPi = 4.0 * 3.14159265358979323846
  real               :: error
  integer            :: it
  !.............................................................................
  call trace_begin ('gmhd_mod::update', 2)
  associate(gn=>self%gn)
  allocate(FourPiGrho(gn(1),gn(2),gn(3)), &
           phi_1(gn(1),gn(2),gn(3)),      &
           phi_2(gn(1),gn(2),gn(3)))            ! allocate scratch memory
  end associate
  !-----------------------------------------------------------------------------
  ! Set up source of gravity
  !-----------------------------------------------------------------------------
  FourPiGrho = FourPi * self%G * self%d
  self%phi => self%mem(:,:,:,self%idx%phi,self%it,1)
  phi_1 => self%phi
  !phi_1 = 0.0
  phi_old => phi_1                              ! initial pointer arrangement
  phi_new => phi_2
  error = 2.*self%tolerance
  it    = 0
  !-----------------------------------------------------------------------------
  ! Do at least min_iter iterations, stop if error<=tolerance or max iterations
  ! reached.  Swap phi_old and phi_new pointers before the next iteration.
  !-----------------------------------------------------------------------------
  do while((it<self%min_iter) .or. &
           (error >= self%tolerance .and. it<self%max_iter))
    call self%laplace%solver (self, phi_old, phi_new, FourPiGrho, error)
    call boundaries%sphere (self, phi_new, 0.0, 0.0)
    if (io%verbose>3) write(*,'(a,i4,2(a,e10.3))'), 'Poisson: Iter=', it, &
      ' Tolerance=', self%tolerance, ' Error=', error
    phi_swap => phi_new; phi_new=>phi_old; phi_old=>phi_swap
    it = it + 1
  enddo
  self%iter = it
  !-----------------------------------------------------------------------------
  ! Check if we reached solution in time
  !-----------------------------------------------------------------------------
  if (it == self%max_iter .and. io%verbose>0) &
    print '(a,i5,2(a,e10.3))', &
     'Poisson: Did not reach convergence in patch. Max iter=', it, &
     ' Tolerance=', self%tolerance, ' Error=', error
  if (it < self%max_iter .and. io%verbose>2) &
    print '(2(a,i4,2x),2(a,e10.3,2x))', &
     'Poisson: Reached convergence in patch. Iter=', it, &
     'Max iter=', self%max_iter, 'Tolerance=', self%tolerance, 'Error=', error
  !-----------------------------------------------------------------------------
  ! Because we swap pointers each iteration, newest solution is in phi_old
  !-----------------------------------------------------------------------------
  self%phi = phi_old
  call trace_end
END SUBROUTINE solve

!===============================================================================
!> Add the gravitional force to the equation of motion
!===============================================================================
SUBROUTINE update (self)
  class(gmhd_t):: self
  real, pointer:: d(:,:,:), phi(:,:,:), lnd(:,:,:), ddn(:,:,:,:)
  real:: ds(3)
  !.............................................................................
  call trace_begin ('gmhd_mod::update')
  if (self%G == 0.0) then
    call self%mhd_t%update
  else
    call self%mhd_t%update
    call self%solve
    d   => self%mem(:,:,:,self%idx%d  ,self%it,1)
    phi => self%mem(:,:,:,self%idx%phi,self%it,1)
    associate (mx=>size(phi,1), my=>size(phi,2), mz=>size(phi,3))
    allocate (lnd(mx,my,mz), ddn(mx,my,mz,3))
    lnd = log(d)
    ddn = down(lnd)
    ddn = exp(ddn)
    ds = self%mesh%d
    self%dpdt(:,:,:,1) = self%dpdt(:,:,:,1) - ddxdn(ds,phi)*ddn(:,:,:,1)
    self%dpdt(:,:,:,2) = self%dpdt(:,:,:,2) - ddydn(ds,phi)*ddn(:,:,:,2)
    self%dpdt(:,:,:,3) = self%dpdt(:,:,:,3) - ddzdn(ds,phi)*ddn(:,:,:,3)
    deallocate (lnd,ddn)
    end associate
  end if
  call trace_end
END SUBROUTINE update

!===============================================================================
SUBROUTINE grav_potential (self, out)
  class(gmhd_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: out
  !...............................................................................
  call trace_begin ('hd_mod::log_density')
  out(1:self%n(1), &
      1:self%n(2), &
      1:self%n(3)) = log(self%mem(self%li(1):self%ui(1), &
       self%li(2):self%ui(2), self%li(3):self%ui(3), self%idx%phi, self%it, 1))
  call trace_end
END SUBROUTINE grav_potential
!===============================================================================
END MODULE gmhd_mod
