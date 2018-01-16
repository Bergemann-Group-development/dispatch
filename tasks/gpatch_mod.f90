!> =============================================================================
!> The gpatch_t data type extends the patch_t with a generic guard zone dnload
!> procedure, which covers the need for all simple, Cartesian mesh based
!> solvers, and eliminates the need for a dnload procedure in each one of them.
!> It also adds restart functionality when initializing patches.
!> =============================================================================
MODULE gpatch_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE download_mod
  USE initial_mod
  USE kinds_mod
  USE data_hub_mod
  USE mesh_mod
  USE units_mod
  USE poisson_mod
  USE omp_timer_mod
  USE bits_mod
  USE timestep_mod
  implicit none
  private
  type, public, extends(patch_t):: gpatch_t
    type(initial_t)  :: initial
    type(data_hub_t) :: data_hub
    type(poisson_t)  :: poisson
    logical:: selfgravity=.false.
    real(8):: d0
  contains
    procedure:: init
    procedure:: dnload
    procedure:: gravity
    procedure:: update
    procedure:: ahead_of
  end type
  real(8), save :: total_mass=0d0, total_volume=0d0
  real(8), save :: restart_time = 0d0
  logical, save :: selfgravity=.false., prediction=.true., detailed_timer=.false.
  integer, save :: verbose=0
  integer, save :: order=1
CONTAINS

!> =============================================================================
!> Add restart functionality to patch initialization
!> =============================================================================
SUBROUTINE init (self)
  class(gpatch_t):: self
  real(kind=KindScalarVar), pointer :: ff(:,:,:,:)
  character(len=64)                 :: filename
  logical                           :: ok
  integer                           :: iostat
  logical, save                     :: first_time=.true.
  real, save                        :: d0=-1.0
  real(kind=KindScalarVar), pointer :: d(:,:,:)
  namelist /selfgravity_params/ selfgravity, prediction, order, d0, verbose
  !-----------------------------------------------------------------------------
  call trace%begin ('gpatch_t%init')
  !$omp critical (poisson_cr)
  if (first_time) then
    rewind (io%input)
    read (io%input, selfgravity_params, iostat=iostat)
    if (io%master) write (io%output, selfgravity_params)
  end if
  !$omp end critical (poisson_cr)
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
    call self%patch_t%init
  end if
  if (io%verbose>0) & 
    print *, self%id, 'gpatch_t%init: ', trim(self%kind), self%nv, self%idx%d, self%idx%phi, real(self%size)
  self%d0 = d0
  call self%initial%init (self%kind)
  !$omp critical (poisson_cr)
  if (first_time) then
    first_time = .false.
    if (self%use_data_hub) then
      !call data_hub%test
    else
      call download%test
    end if
  end if
  !$omp end critical (poisson_cr)
  !-----------------------------------------------------------------------------
  ! Try reading a restart snapshot.
  !-----------------------------------------------------------------------------
  call self%input (ok)
  if (ok) then
    restart_time = self%time
  !-----------------------------------------------------------------------------
  ! If the restart parameter is not set, or the file does not exist, initialize
  ! the patch from scratch, resetting other necessary counters only if positive
  ! restart_time indicates that restart files were read for previous patches.
  !-----------------------------------------------------------------------------
  else
    self%mem = 0.0
    ff => self%mem(:,:,:,:,1,1)
    call self%initial%condition (self%mesh, ff, self%idx)
    !---------------------------------------------------------------------------
    ! The first time, all patches that exist also have good GC values, and we
    ! can freely dnload both densities and potential and then get a first
    ! solution of the Poisson equation immediately.  The root grid does this
    ! first, which means subsequent patches from initial refinement get their
    ! BCs set correctly, for an immediately good startup.
    !---------------------------------------------------------------------------
    if (selfgravity) then
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
      call self%dnload (only=self%idx%d)
      call self%dnload (only=self%idx%phi)
      call self%gravity
    end if
  end if
  !-----------------------------------------------------------------------------
  ! If restart succeeded for some patch, use it to set time and output cadence.
  ! FIXME: This may fail in some corner cases when running on multiple MPI ranks,
  ! if one rank happens to handle only patches that do not exist yet.
  !-----------------------------------------------------------------------------
  if (restart_time > 0d0) then
    self%time = restart_time
    self%t(self%it) = self%time
    self%out_next = (int(self%time/io%out_time)+1)*io%out_time
    self%iout = self%restart+1
    !$omp master
    write (filename,'(a,i5.5,"/")') trim(io%outputname), self%iout
    if (io%iodir/=self%iout) call io%mkdir (trim(filename))
    io%iodir = self%iout
    !$omp end master
  end if
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Add default guard zone download functionality.  This may be overloaded in
!> experiment_mod, if more specialized guard zone loading is required.
!> =============================================================================
SUBROUTINE dnload (self, only)
  class(gpatch_t):: self
  integer, optional:: only
  !-----------------------------------------------------------------------------
  if (self%use_data_hub) then
    call self%data_hub%update (self%link, only=only)
  else
    call download%download_link (self%link, only=only)
  end if
END SUBROUTINE dnload

!> =============================================================================
!> Prepare the external forces for the pde call, in particular selfgravity. At
!> this point the nbors are declared up-to-date with respect to everything,
!> including their potentials, which we need for guard zones.  They may in fact
!> be providing only their preliminary potentials in their (nt-1) time slots,
!> while their older time slots should by now contain the converged potential.
!> =============================================================================
SUBROUTINE gravity (self)
  class(gpatch_t):: self
  !.............................................................................
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, phi
  real:: dmin, dmax
  integer:: iter, ix ,iy, iz
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. selfgravity) return
  call trace%begin ('gpatch_t%gravity', itimer=itimer)
  d   => self%mem(:,:,:,self%idx%d  ,self%it,1)
  phi => self%mem(:,:,:,self%idx%phi,self%it,1)
  if (self%d0==-1.0) then
    self%d0 = total_mass/total_volume
    print *, self%id, 'setting self%d0 =', self%d0, self%faver(d)
  end if
  dmin = self%fminval (d)
  dmax = self%fmaxval (d)
  if (dmin==0.0) &
    print *, self%id, 'WARNING: gpatch_t%gravity, dmin, dmax, it =', dmin, dmax, self%it
  call self%poisson%update (self, phi, d, self%d0, detailed_timer)
  !-----------------------------------------------------------------------------
  ! compute estimate of dphi/dt
  !-----------------------------------------------------------------------------
  call time_derivs (self)
  call trace%end (itimer)
END SUBROUTINE gravity

!> =============================================================================
!> Having advanced the MHD state to the new time, we make a preliminary
!> call to the Poisson solver, based on extrapolation of the nbor potentials.
!> We then provide this preliminary potential to nbors, until we have a
!> better offer.  The poisson_time controls if an update affects the self%time
!> or the poisson_time; the two are updated alternatingly
!> =============================================================================
SUBROUTINE update (self)
  class(gpatch_t):: self
  !-----------------------------------------------------------------------------
  if (.not. selfgravity) return
  call trace%begin ('gpatch_t%update')
  if (self%poisson_time == self%time) then
    call timestep%update(self%id, self%iit, self%dt, self%mem, self%mesh)
    call self%rotate
    self%rotated = .true.
    if (verbose>0) &
      print *, self%id, 'gpatch_t%update updated self%time to', self%time
  else
    call self%dnload (only=self%idx%d)
    call self%dnload (only=self%idx%phi)
    call self%gravity
    self%poisson_time = self%time
    self%rotated = .true.
    if (verbose>0) &
      print *, self%id, 'gpatch_t%update updated self%poisson_time to', self%poisson_time
  end if
  call trace%end
END SUBROUTINE update

!> =============================================================================
!> Compute an approximate, conservative estimate of dphi/dt, which will cause
!> the timestep routine to make a good guess at the next phi, at the same time
!> as it is updating the other variables.  The 1st order prediction is in
!> practice better at reducing the number of iterations than the 2nd order.
!> =============================================================================
SUBROUTINE time_derivs (self)
  class(gpatch_t):: self
  integer:: it1, it2, it3, nt, l(3), u(3)
  real(8):: t1, t2, t3
  integer, save:: itimer=0
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: dphidt
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('gpatch_t%time_derivs', itimer=itimer)
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
  class(gpatch_t):: self
  class(patch_t):: patch
  !-----------------------------------------------------------------------------
  ! Level differences larger than 1 should never happen!
  !-----------------------------------------------------------------------------
  select type (patch)
  class is (gpatch_t)
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

END MODULE gpatch_mod
