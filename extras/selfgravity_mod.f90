!> =============================================================================
!> Add selfgravity to any solver
!> =============================================================================
MODULE selfgravity_mod
  USE io_mod
  USE os_mod
  USE trace_mod
  USE patch_mod
  USE download_mod
  USE initial_mod
  USE kinds_mod
  USE data_hub_mod
  USE mesh_mod
  USE units_mod
  USE poisson_mod
  USE timer_mod
  USE omp_timer_mod
  USE bits_mod
  implicit none
  private
  type, public, extends(patch_t):: selfgravity_t
    type(initial_t)  :: initial
    type(data_hub_t) :: data_hub
    type(poisson_t)  :: poisson
    logical:: selfgravity=.false.
    real(8):: d0
  contains
    procedure:: init
    procedure:: output
    procedure:: pre_update
    procedure:: post_update
    procedure:: ahead_of
  end type
  real(8), save :: total_mass=0d0, total_volume=0d0
  real(8), save :: restart_time = 0d0
  logical, save :: selfgravity=.false., prediction=.true., detailed_timer=.false.
  integer, save :: verbose=0
  integer, save :: order=1
CONTAINS

!> =============================================================================
!> Add restart functionality to patch initialization.  This %init routine is
!> called for both internal and virtual patches, which both can be initialized
!> fully by the self%ininitial%condition() procedure. When instead calling the
!> self%input() procedure, we need to make a decision:  Either we must make sure
!> that the data_io%input() routine can handle off-rank reading, or we need to
!> nevertheless return ok=.true., and instead rely on a negative time to force
!> any task that relies on the patch for guard zone values to wait until data
!> has been recieved with MPI from the owner rank.
!> =============================================================================
SUBROUTINE init (self)
  class(selfgravity_t):: self
  real(kind=KindScalarVar), pointer :: ff(:,:,:,:)
  character(len=64)                 :: filename
  logical                           :: ok
  integer                           :: iostat
  logical, save                     :: first_time=.true., first_poisson=.true.
  real, save                        :: d0=-1.0
  real(kind=KindScalarVar), pointer :: d(:,:,:)
  namelist /selfgravity_params/ selfgravity, prediction, order, d0, verbose
  !-----------------------------------------------------------------------------
  call trace%begin ('selfgravity_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, selfgravity_params, iostat=iostat)
    if (io%master) write (io%output, selfgravity_params)
  end if
  !$omp end critical (input_cr)
  if (selfgravity) then
    self%selfgravity = selfgravity
    if (trim(self%kind) /= 'zeus_mhd_patch') self%nv = self%nv + 1
    self%idx%phi = self%nv
    if (self%nw == 1) then
      self%nv = self%nv + 1
      self%idx%dphi = self%nv
    end if
    call self%patch_t%init
    call self%poisson%init
  else
    self%idx%phi = -1
    call self%patch_t%init
  end if
  allocate (self%force_per_unit_volume(self%gn(1),self%gn(2),self%gn(3),3))
  self%d0 = d0
  !$omp critical (poisson_cr)
  if (first_poisson) then
    first_poisson = .false.
    if (self%use_data_hub) then
      !call data_hub%test
    else
      call download%test
    end if
  end if
  !$omp end critical (poisson_cr)
  if (selfgravity) then
    !-------------------------------------------------------------------------
    ! Sum up the total mass and total volume, to get the average density d0.
    ! FIXME: This needs to be modified for MPI
    !-------------------------------------------------------------------------
    if (self%id==1) then
      d => self%mem(:,:,:,self%idx%d,self%it,1)
      !$omp atomic
      total_mass = total_mass + self%fsum(d)*product(self%ds)
      !$omp atomic
      total_volume = total_volume + product(self%size)
    end if
    !-------------------------------------------------------------------------
    ! Sum up the total mass and total volume, to get the average density d0.
    ! FIXME: This needs to be modified for MPI
    !-------------------------------------------------------------------------
    if (self%time == 0.0 .and. self%id==1) then
      d => self%mem(:,:,:,self%idx%d,self%it,1)
      !$omp atomic
      total_mass = total_mass + self%fsum(d)*product(self%ds)
      !$omp atomic
      total_volume = total_volume + product(self%size)
    end if
  end if
  call trace%end()
END SUBROUTINE init

!> =============================================================================
!> Interface routine to the data_io output procedure
!> =============================================================================
SUBROUTINE update (self)
  class(selfgravity_t):: self
  !-----------------------------------------------------------------------------
END SUBROUTINE update

!> =============================================================================
!> Interface routine to the data_io output procedure
!> =============================================================================
SUBROUTINE output (self)
  class(selfgravity_t):: self
  !-----------------------------------------------------------------------------
END SUBROUTINE output

!> =============================================================================
!> Prepare the external forces for the pde call, in particular selfgravity. At
!> this point the nbors are declared up-to-date with respect to everything,
!> including their potentials, which we need for guard zones.  They may in fact
!> be providing only their preliminary potentials in their (nt-1) time slots,
!> while their older time slots should by now contain the converged potential.
!> =============================================================================
SUBROUTINE pre_update (self)
  class(selfgravity_t):: self
  !.............................................................................
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, phi
  real:: dmin, dmax
  integer:: iter, ix ,iy, iz
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. selfgravity) return
  call trace%begin ('selfgravity_t%pre_update', itimer=itimer)
  d   => self%mem(:,:,:,self%idx%d  ,self%it,1)
  phi => self%mem(:,:,:,self%idx%phi,self%it,1)
  if (self%d0==-1.0) then
    self%d0 = total_mass/total_volume
    print *, self%id, 'setting self%d0 =', self%d0, self%faver(d)
  end if
  dmin = self%fminval (d)
  dmax = self%fmaxval (d)
  if (dmin==0.0) &
    print *, self%id, 'WARNING: selfgravity_t%gravity, dmin, dmax, it =', dmin, dmax, self%it
  call self%poisson%update (self, phi, d, self%d0, detailed_timer)
  !-----------------------------------------------------------------------------
  ! compute estimate of dphi/dt
  !-----------------------------------------------------------------------------
  call time_derivs (self)
  call trace%end (itimer)
END SUBROUTINE pre_update

!> =============================================================================
!> Having advanced the MHD state to the new time, we make a preliminary
!> call to the Poisson solver, based on extrapolation of the nbor potentials.
!> We then provide this preliminary potential to nbors, until we have a
!> better offer.  The poisson_time controls if an update affects the self%time
!> or the poisson_time; the two are updated alternatingly
!> =============================================================================
SUBROUTINE post_update (self)
  class(selfgravity_t):: self
  !-----------------------------------------------------------------------------
  if (.not. selfgravity) return
  call trace%begin ('selfgravity_t%post_update')
  if (selfgravity) then
    if (self%poisson_time == self%time) then
      !call timestep%update(self%id, self%iit, self%dt, self%mem, self%mesh)
      call self%rotate
      self%rotated = .true.
      if (verbose>0) &
        print *, self%id, 'selfgravity_t%update updated self%time to', self%time
    else
      call self%dnload (only=self%idx%d)
      call self%dnload (only=self%idx%phi)
      call self%pre_update
      self%poisson_time = self%time
      self%rotated = .true.
      if (verbose>0) &
        print *, self%id, 'selfgravity_t%update updated self%poisson_time to', self%poisson_time
    end if
  end if
  call trace%end()
END SUBROUTINE post_update

!> =============================================================================
!> Compute an approximate, conservative estimate of dphi/dt, which will cause
!> the timestep routine to make a good guess at the next phi, at the same time
!> as it is updating the other variables.  The 1st order prediction is in
!> practice better at reducing the number of iterations than the 2nd order.
!> =============================================================================
SUBROUTINE time_derivs (self)
  class(selfgravity_t):: self
  integer:: it1, it2, it3, nt, l(3), u(3)
  real(8):: t1, t2, t3
  integer, save:: itimer=0
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: dphidt
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('selfgravity_t%time_derivs', itimer=itimer)
  if (self%nw == 2) then
    dphidt => self%mem(:,:,:,self%idx%phi,self%it,2)
  else
    dphidt => self%mem(:,:,:,self%idx%dphi,self%it,1)
  end if

  if (prediction) then
    nt = self%nt
    it1 = self%iit(nt-1)
    it2 = self%iit(nt-2)
    it3 = self%iit(nt-3)
    t1 = self%t(it1)
    t2 = self%t(it2)
    t3 = self%t(it3)
    dphidt = 0.0
    if (it2/=it1 .and. t1/=t2) then
      l = self%mesh%li
      u = self%mesh%ui
      if (order==2 .and. t2/=t3) then
        dphidt(l(1):u(1),l(2):u(2),l(3):u(3)) = &
                 (self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%phi,it1,1) - &
                  self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%phi,it2,1))/(t1-t2)*1.5 - &
                 (self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%phi,it2,1) - &
                  self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%phi,it3,1))/(t2-t3)*0.5
      else
        dphidt(l(1):u(1),l(2):u(2),l(3):u(3)) = &
                 (self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%phi,it1,1) - &
                  self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%phi,it2,1))/(t1-t2)
      end if
      if (verbose>2) print *,self%id,'d(phi)/dt min,max,nv =', &
          self%fminval(dphidt), self%fmaxval(dphidt), size(self%mem,4)
    end if
  else
    dphidt = 0.0
  end if
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE time_derivs

!> =============================================================================
!> This functions checks if self is too much ahead of the nbor 'patch'.  If it
!> returns .true. it means the 'self' task is not ready to be updated.
!> =============================================================================
LOGICAL FUNCTION ahead_of (self, patch)
  class(selfgravity_t):: self
  class(patch_t):: patch
  !-----------------------------------------------------------------------------
  ! Level differences larger than 1 should never happen!
  !-----------------------------------------------------------------------------
  select type (patch)
  class is (selfgravity_t)
  if (abs(self%level-patch%level)>1) then
    ahead_of = .false.
    print *,'WARNING: level difference larger than one:', &
      self%id, self%level, patch%id, patch%level
  !-----------------------------------------------------------------------------
  ! If wee need to solve the Poisson equation, we ignore nbors with higher ot
  ! the same resolution, and require nbords with lower resolution to be in front
  ! with their Poisson solutions.
  !-----------------------------------------------------------------------------
  else if (self%poisson_time < self%time) then
    if (self%level <= patch%level) then
      ahead_of = .false.
    else
      !ahead_of = self%time > patch%poisson_time + patch%grace*patch%dtime
      ahead_of = self%time > patch%time + patch%grace*patch%dtime
    end if
    if (ahead_of .and. verbose>1) &
      print 1, self%id, self%level, 'at', self%time, &
        'cannot update because', patch%id, patch%level, 'has too old Poisson time', &
        patch%poisson_time + patch%grace*patch%dtime
      1 format(i6,i3,2x,a,g15.6,2x,a,i6,i3,2x,a,g15.6)
  !-----------------------------------------------------------------------------
  ! If we are going to update the dynamic time, we need the nbor to be in front.
  !-----------------------------------------------------------------------------
  else
    ahead_of = self%time > patch%time + patch%grace*patch%dtime
    if (ahead_of .and. verbose>1) &
      print 1, self%id, self%level, 'at', self%time, &
        'cannot update because', patch%id, patch%level, 'has too old dynamic time', &
        patch%time + patch%grace*patch%dtime
  end if
  end select
END FUNCTION ahead_of

END MODULE selfgravity_mod
