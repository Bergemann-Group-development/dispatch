! *******************************************************************
!=============================================================================
!> Ng acceleration of the statistical-equilibrium iteration.
!>
!> The extrapolation is the standard three-difference Ng in log-population
!> space.  With n1..n4 the last four iterates of ln b and
!>
!>   d0 = n4 - n3,   d1 = d0 - (n3 - n2),   d2 = d0 - (n2 - n1)
!>
!> the coefficients (a,b) minimise |d0 - a*d1 - b*d2| in the least-squares
!> sense, and the accelerated iterate is
!>
!>   n_new = (1 - a - b)*n4 + a*n3 + b*n2.
!>
!> Working in ln b is what makes this safe for a multiplicative update: the
!> extrapolated population is strictly positive for any (a,b).
!>
!> The normal equations are built from inner products of the difference
!> vectors, so each patch contributes a partial sum and the coefficients are
!> formed once from the total.  A patch cannot know when the total is complete,
!> so accumulation and application are one iteration apart: every iteration
!> fills a *current* Gram matrix, and at the start of the next one the filled
!> matrix is promoted to *done* and becomes eligible to be used.
!>
!> Whether it is actually used is the job of ng_should_fire.  Ng extrapolates a
!> linear fixed-point map, and the iteration is only that map in its asymptotic
!> regime; applying it elsewhere is at best wasted work and at worst a large
!> step in the wrong direction.  See the gate conditions there.
!=============================================================================
MODULE ng_mod
  USE io_mod
  USE trace_mod
  USE console_mod
  implicit none
  private

  !-----------------------------------------------------------------------------
  ! Gate tunables.  Exposed through &nlte_params so a run can be made more or
  ! less eager without recompiling; the defaults are the ones argued for in
  ! fe607_convergence.md.
  !-----------------------------------------------------------------------------
  real(8), public, save :: ng_rate_min=0.5d0  !< don't bother below this contraction rate
  real(8), public, save :: ng_rate_tol=0.15d0 !< |r - r_prev|/r must be under this
  real(8), public, save :: ng_clamp_max=1d-3  !< max fraction of corrections at the limiter
  !> Minimum iterations that must remain, at the measured rate, for the
  !> extrapolation to be worth its risk.  Expressing the headroom in iterations
  !> rather than as a bare multiple of convlim makes it adapt to the rate: a slow
  !> mode still has far to go even close to convlim, a fast one does not.
  integer, public, save :: ng_min_iters_left=3
  real(8), public, save :: ng_cond_min=1d-6  !< scaled Gram determinant floor
  real(8), public, save :: ng_max_change=1d1 !< reject a step moving a level by more
  integer, public, save :: ng_min_gap=4      !< iterations between applications

  !-----------------------------------------------------------------------------
  !> State shared by every patch copy of one atom.  nlte_patch allocates a single
  !> target on the master atom; the pointer component survives the atom copy, so
  !> all patches accumulate into the same object (as the three separate pointers
  !> this replaces did).
  !-----------------------------------------------------------------------------
  type, public :: ng_state_t
    !> Gram matrix being filled this iteration, and the completed one from last
    real(8) :: d_cur(0:2,1:2) = 0d0, d_done(0:2,1:2) = 0d0
    !> |d0|^2 over all patches: a true global norm of the last population change,
    !> and the only convergence signal here that is neither lagged (like the
    !> mpi_counter-reduced `global`) nor a single cell (like `local`)
    real(8) :: s_cur=0d0, s_done=0d0, s_prev=-1d0
    real(8) :: r=-1d0, r_prev=-1d0 !< contraction rate and its predecessor
    real(8) :: a=0d0, b=0d0        !< extrapolation coefficients
    !> corrections that hit the update limiter, counted over all patches; the
    !> iteration is not a linear map where this is a large fraction
    integer(8) :: nclamp_cur=0, ntot_cur=0, nclamp_done=0, ntot_done=0
    integer :: iter_promoted=-1      !< last iteration the roll-over ran
    integer :: iter_decided=-1       !< last iteration the gate was evaluated
    integer :: iter_fired=-huge(1)/2 !< last iteration an extrapolation was applied
    logical :: coeff_valid=.false.   !< are a,b usable this iteration
  end type ng_state_t

  public :: ng_promote, ng_should_fire, ng_coefficients, ng_accumulate
  public :: ng_extrapolate, ng_count_clamped, ng_condition

CONTAINS

!===============================================================================
!> Conditioning of the least-squares fit: the Gram determinant of (d1,d2)
!> normalised by the product of their squared norms, i.e. 1 - cos^2 of the angle
!> between them.  It is >= 0 by Cauchy-Schwarz and approaches 0 as the two
!> difference vectors become parallel -- a single dominant mode, where the
!> two-parameter fit is degenerate and the 2x2 solve amplifies roundoff.
!===============================================================================
REAL(8) FUNCTION ng_condition (self) RESULT(cond)
  type(ng_state_t), intent(in) :: self
  real(8) :: den, scale
  !-----------------------------------------------------------------------------
  associate (d => self%d_done)
    den   = d(1,1)*d(2,2) - d(1,2)*d(2,1)
    scale = d(1,1)*d(2,2)
  end associate
  if (scale > 0d0) then
    cond = den / scale
  else
    cond = 0d0
  endif
END FUNCTION ng_condition

!===============================================================================
!> Roll the current accumulation into the completed one and derive the
!> contraction rate from it.  Idempotent within an iteration: the first patch to
!> arrive does the work, the rest see iter_promoted and return.  Must be called
!> under a lock.
!===============================================================================
SUBROUTINE ng_promote (self, iter)
  type(ng_state_t), intent(inout) :: self
  integer, intent(in) :: iter
  !-----------------------------------------------------------------------------
  ! Patches run their own iteration counts and reach this at different times, so
  ! the roll has to be monotone: a patch arriving late, still on an earlier
  ! iteration, must not roll a second time and discard a live accumulation.
  if (iter <= self%iter_promoted) return
  self%iter_promoted = iter

  self%d_done      = self%d_cur
  self%s_done      = self%s_cur
  self%nclamp_done = self%nclamp_cur
  self%ntot_done   = self%ntot_cur
  self%d_cur       = 0d0
  self%s_cur       = 0d0
  self%nclamp_cur  = 0
  self%ntot_cur    = 0

  ! |d0| falls at the contraction rate, so the ratio of successive norms is the
  ! rate itself.  An iteration that extrapolated contributes nothing (the
  ! differences would span two different sequences), and its empty accumulation
  ! must not be written into the history as a zero -- leave the rate stale and
  ! let the next real accumulation resume from the last genuine norm.
  if (self%s_done > 0d0) then
    self%r_prev = self%r
    if (self%s_prev > 0d0) then
      self%r = sqrt(self%s_done / self%s_prev)
    else
      self%r = -1d0
    endif
    self%s_prev = self%s_done
  endif
  self%coeff_valid = .false.
END SUBROUTINE ng_promote

!===============================================================================
!> The automatic trigger.  Every condition below is a statement about whether
!> the Ng model applies, not about how well the run is going.
!===============================================================================
LOGICAL FUNCTION ng_should_fire (self, iter, err, convlim, reason) RESULT(fire)
  type(ng_state_t), intent(in) :: self
  integer, intent(in) :: iter
  real(8), intent(in) :: err, convlim
  ! which condition declined, for the DIAGNG trace; 0 = fired
  character(len=*), intent(out), optional :: reason
  real(8) :: clamp_frac
  !-----------------------------------------------------------------------------
  fire = .false.
  if (present(reason)) reason = 'none'

  ! Nothing accumulated yet, or the history is not four iterates deep.
  if (present(reason)) reason = 'no-history'
  if (self%s_done <= 0d0 .or. self%s_prev <= 0d0) return

  ! An extrapolation invalidates the difference history: oldb still holds
  ! un-extrapolated iterates while b holds an extrapolated one, so d0 spans two
  ! different sequences.  Three clean iterations are needed to refill it.
  if (present(reason)) reason = 'min-gap'
  if (iter - self%iter_fired < ng_min_gap) return

  ! Is the update a linear map?  While a significant fraction of the level
  ! corrections is being cut by the limiter, it is not, and the difference
  ! vectors describe the clamp rather than the iteration.  This is the condition
  ! that rejects the Fe III limiter plateau, which a rate test alone accepts
  ! because a pinned error looks perfectly stationary.
  if (present(reason)) reason = 'clamped'
  if (self%ntot_done <= 0) return
  clamp_frac = real(self%nclamp_done, 8) / real(self%ntot_done, 8)
  if (clamp_frac > ng_clamp_max) return

  ! Is it in the asymptotic regime?  A stationary rate is what "asymptotic"
  ! means: the error has collapsed onto the dominant modes.
  if (present(reason)) reason = 'not-stationary'
  if (self%r <= 0d0 .or. self%r_prev <= 0d0) return
  if (abs(self%r - self%r_prev) > ng_rate_tol * self%r) return

  ! Is it slow enough to be worth accelerating?  Below ng_rate_min the plain
  ! iteration finishes before the extrapolation repays its risk.
  if (present(reason)) reason = 'fast-enough'
  if (self%r < ng_rate_min) return

  ! Diverging, or on the edge of it: extrapolating along a growing mode
  ! amplifies it.  Leave that to the guard rails.
  if (present(reason)) reason = 'diverging'
  if (self%r >= 1d0) return

  ! Is there enough left to gain?  |err| decays as r per iteration, so reaching
  ! convlim takes ln(convlim/|err|)/ln(r) more of them; require at least
  ! ng_min_iters_left.  Near convlim an over-correction can push the run back out
  ! of the convergence band for longer than the acceleration saves.
  if (present(reason)) reason = 'no-headroom'
  if (abs(err) < convlim * self%r**(-ng_min_iters_left)) return

  ! Is the least-squares problem well conditioned?  A bare /= 0 test on the
  ! determinant almost never catches a degenerate fit; compare it to its own
  ! scale instead.
  if (present(reason)) reason = 'ill-conditioned'
  if (ng_condition(self) <= ng_cond_min) return

  if (present(reason)) reason = 'FIRE'
  fire = .true.
END FUNCTION ng_should_fire

!===============================================================================
!> Solve the 2x2 normal equations for the extrapolation coefficients.  Only
!> called once ng_should_fire has established that den is safely non-degenerate.
!===============================================================================
SUBROUTINE ng_coefficients (self)
  type(ng_state_t), intent(inout) :: self
  real(8) :: den
  !-----------------------------------------------------------------------------
  associate (d => self%d_done)
    den = d(1,1)*d(2,2) - d(1,2)*d(2,1)
    self%a = (d(0,1)*d(2,2) - d(0,2)*d(2,1)) / den
    self%b = (d(0,2)*d(1,1) - d(0,1)*d(1,2)) / den
  end associate
  self%coeff_valid = .true.
END SUBROUTINE ng_coefficients

!===============================================================================
!> Add this patch's contribution to the Gram matrix and to |d0|^2.
!>
!> d(1,2) and d(2,1) are the same inner product; both are kept so the solve
!> above reads as the textbook 2x2 system.  Must be called under a lock: several
!> patches reach this concurrently and every line is a read-modify-write on
!> shared state.
!===============================================================================
SUBROUTINE ng_accumulate (self, bpop, oldb, nlevel)
  type(ng_state_t), intent(inout) :: self
  real(8), dimension(:,:,:,:), intent(in) :: bpop
  real(8), dimension(:,:,:,:,:), intent(in) :: oldb
  integer, intent(in) :: nlevel
  real(8), dimension(:), allocatable :: d0, d1, d2
  real(8) :: n1, n2, n3, n4
  integer :: i, x, y, z, k, nvec, m(3)
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('ng_mod%accumulate', itimer=itimer)
  m = [size(bpop,1), size(bpop,2), size(bpop,3)]
  nvec = nlevel * product(m)
  allocate(d0(nvec), d1(nvec), d2(nvec))

  i = 0
  do z=1, m(3)
  do y=1, m(2)
  do x=1, m(1)
  do k=1, nlevel
    i = i + 1
    n1 = log(oldb(x,y,z,k,1))
    n2 = log(oldb(x,y,z,k,2))
    n3 = log(oldb(x,y,z,k,3))
    n4 = log(bpop(x,y,z,k))
    d0(i) = n4 - n3
    d1(i) = d0(i) - (n3 - n2)
    d2(i) = d0(i) - (n2 - n1)
  end do
  end do
  end do
  end do

  !$omp critical (ng_shared)
  associate (d => self%d_cur)
    d(0,1) = d(0,1) + dot_product(d0, d1)
    d(0,2) = d(0,2) + dot_product(d0, d2)
    d(1,1) = d(1,1) + dot_product(d1, d1)
    d(1,2) = d(1,2) + dot_product(d1, d2)
    d(2,1) = d(2,1) + dot_product(d2, d1)
    d(2,2) = d(2,2) + dot_product(d2, d2)
  end associate
  self%s_cur = self%s_cur + dot_product(d0, d0)
  !$omp end critical (ng_shared)

  call trace%end(itimer)
END SUBROUTINE ng_accumulate

!===============================================================================
!> Record how many level corrections this patch had cut by the update limiter.
!===============================================================================
SUBROUTINE ng_count_clamped (self, nclamp, ntot)
  type(ng_state_t), intent(inout) :: self
  integer(8), intent(in) :: nclamp, ntot
  !-----------------------------------------------------------------------------
  !$omp critical (ng_shared)
  self%nclamp_cur = self%nclamp_cur + nclamp
  self%ntot_cur   = self%ntot_cur   + ntot
  !$omp end critical (ng_shared)
END SUBROUTINE ng_count_clamped

!===============================================================================
!> Apply the extrapolation to one patch, renormalising to conserve particles.
!>
!> The step is measured before anything is written, and abandoned if any level
!> would move by more than ng_max_change.  A converging sequence does not need
!> an order-of-magnitude correction; if the fit asks for one, the linear model
!> behind it does not hold and the honest response is to keep the unaccelerated
!> populations.  This is a local veto -- it can leave one patch unaccelerated
!> while its neighbours are -- which is deliberate: the alternative is
!> propagating a step already known to be unphysical.
!===============================================================================
SUBROUTINE ng_extrapolate (self, bpop, oldb, nstar, totn, nlevel, emax, epos, applied)
  type(ng_state_t), intent(in) :: self
  real(8), dimension(:,:,:,:), intent(inout) :: bpop
  real(8), dimension(:,:,:,:,:), intent(in) :: oldb
  real(8), dimension(:,:,:,:), intent(in) :: nstar
  real(8), dimension(:,:,:), intent(in) :: totn
  integer, intent(in) :: nlevel
  real(8), intent(out) :: emax
  integer, intent(out) :: epos(4)
  logical, intent(out) :: applied
  real(8) :: n2, n3, n4, e, enew, worst
  integer :: x, y, z, k, m(3)
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('ng_mod%extrapolate', itimer=itimer)
  emax = 0d0
  epos = 0
  applied = .false.
  m = [size(bpop,1), size(bpop,2), size(bpop,3)]

  associate (a => self%a, b => self%b)
    ! pass 1: how far does the fit want to move things
    worst = 0d0
    do z=1, m(3)
    do y=1, m(2)
    do x=1, m(1)
    do k=1, nlevel
      n2 = log(oldb(x,y,z,k,2))
      n3 = log(oldb(x,y,z,k,3))
      n4 = log(bpop(x,y,z,k))
      e  = (1d0 - a - b)*n4 + a*n3 + b*n2
      worst = max(worst, abs(e - n4))
    end do
    end do
    end do
    end do
    if (exp(worst) - 1d0 > ng_max_change) then
      if (console_on(ch_m3d, 1)) then
        !$omp critical (ng_print)
        write(console_unit,'(A,ES10.3,A,ES10.3,A,2F9.4)') &
          'NG vetoed: step ', exp(worst) - 1d0, ' exceeds ng_max_change ', &
          ng_max_change, ', a,b =', a, b
        flush(console_unit)
        !$omp end critical (ng_print)
      end if
      call trace%end(itimer)
      return
    endif

    ! pass 2: commit
    do z=1, m(3)
    do y=1, m(2)
    do x=1, m(1)
      do k=1, nlevel
        n2 = log(oldb(x,y,z,k,2))
        n3 = log(oldb(x,y,z,k,3))
        n4 = log(bpop(x,y,z,k))
        e  = (1d0 - a - b)*n4 + a*n3 + b*n2
        bpop(x,y,z,k) = exp(e)
        enew = exp(e - n4) - 1d0
        if (abs(enew) > abs(emax)) then
          emax = enew
          epos = [k, x, y, z]
        endif
      end do
      ! strict macroscopic particle conservation, as the SE update does
      bpop(x,y,z,:) = bpop(x,y,z,:) * (totn(x,y,z) / sum(bpop(x,y,z,:) * nstar(x,y,z,:)))
    end do
    end do
    end do
  end associate

  applied = .true.
  call trace%end(itimer)
END SUBROUTINE ng_extrapolate

END MODULE ng_mod
