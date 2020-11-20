!===============================================================================
!> Add selfgravity to any solver
!===============================================================================
MODULE selfgravity_mod
  USE io_mod
  USE io_unit_mod
  USE os_mod
  USE trace_mod
  USE patch_mod
  USE gpatch_mod
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
  USE vector_ops
  USE timestep_mod
  USE scratch_mod
  implicit none
  private
  type, public:: selfgravity_t
    type(initial_t)  :: initial
    type(poisson_t)  :: poisson
    logical:: selfgravity=.false.
    real(8):: d0
    real(8):: poisson_time=-1.0
    class(gpatch_t), pointer:: patch=>null()
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: pre_update
    procedure:: post_update
    procedure:: ahead_of
  end type
  real(8), save :: total_mass=0d0, total_volume=0d0
  real(8), save :: restart_time = 0d0
  logical, save :: selfgravity=.false., prediction=.true.
  integer, save :: verbose=0
  integer, save :: order=1
CONTAINS

!===============================================================================
!> Further initialisation of the self-gravity module
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(selfgravity_t):: self
  class(gpatch_t), pointer:: patch
  real(kind=KindScalarVar), pointer :: ff(:,:,:,:)
  character(len=64)                 :: filename
  logical                           :: ok
  integer                           :: iostat
  logical, save                     :: first_time=.true., first_poisson=.true.
  real, save                        :: d0=-1.0
  real(kind=KindScalarVar), pointer :: d(:,:,:)
  namelist /selfgravity_params/ selfgravity, prediction, order, d0, verbose
  !-----------------------------------------------------------------------------
  call trace%begin ('selfgravity_t%pre_init')
  self%patch => patch
  !
  call io%lock%set ('selfgravity_t%pre_init')
  if (first_time) then
    rewind (io%input)
    read (io%input, selfgravity_params, iostat=iostat)
    call io%namelist_warning('selfgravity_params', iostat)
    write (io%namelist, selfgravity_params)
      first_time = .false.
  end if
  call io%lock%unset ('selfgravity_t%pre_init')
  !.............................................................................
  self%d0 = d0
  if (selfgravity) then
    self%selfgravity = selfgravity
    patch%idx%phi = 0
    call patch%idx%request_index (patch%idx%phi, patch%nv)
    if (patch%nw == 2) then
      patch%idx%dphi = patch%idx%phi
    else
      patch%idx%dphi = 0
      call patch%idx%request_index (patch%idx%dphi, patch%nv)
    end if
    !...........................................................................
    if (verbose > 0) then
      write (io_unit%log,'(a,i6,4i4)') &
       'selfgravity_t%pre_init: id, nv, nw, phi, dphi =', &
        patch%id, patch%nv, patch%nw, patch%idx%phi, patch%idx%dphi
    end if
    !...........................................................................
    call self%poisson%init
  else
    patch%idx%phi = -1
  end if
  !
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> Further initialisation of the self-gravity module
!===============================================================================
SUBROUTINE init (self)
  class(selfgravity_t):: self
  !.............................................................................
  class(patch_t), pointer:: patch
  real(kind=KindScalarVar), pointer :: ff(:,:,:,:)
  character(len=64)                 :: filename
  logical                           :: ok
  integer                           :: iostat
  logical, save                     :: first_poisson=.true.
  real(kind=KindScalarVar), pointer :: d(:,:,:)
  real                              :: dm
  !-----------------------------------------------------------------------------
  call trace%begin ('selfgravity_t%init')
  patch => self%patch
  !$omp critical (poisson_cr)
  if (first_poisson) then
    first_poisson = .false.
    call download%test
  end if
  !$omp end critical (poisson_cr)
  if (selfgravity) then
    !---------------------------------------------------------------------------
    ! Sum up the total mass and total volume, to get the average density d0.
    ! FIXME: This needs to be modified for MPI
    !---------------------------------------------------------------------------
    if (patch%time == 0.0) then
      d => patch%mem(:,:,:,patch%idx%d,patch%it,1)
      dm = patch%fsum(d)*product(patch%ds)
      !$omp atomic
      total_mass = total_mass + dm
      !$omp atomic
      total_volume = total_volume + product(patch%size)
      !.........................................................................
      if (verbose > 0) then
        write (io_unit%log,'(a,i6,1p,3e12.2)') &
         'selfgravity_t%init: id, dm =', patch%id, dm, total_mass, total_volume
      end if
      !.........................................................................
    end if
  end if
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Prepare the external forces for the pde call, in particular selfgravity. At
!> this point the nbors are declared up-to-date with respect to everything,
!> including their potentials, which we need for guard zones.  They may in fact
!> be providing only their preliminary potentials in their (nt-1) time slots,
!> while their older time slots should by now contain the converged potential.
!===============================================================================
SUBROUTINE pre_update (self)
  class(selfgravity_t):: self
 !------------------------------------------------------------------------------
  if (.not. selfgravity) return
  call trace%begin ('selfgravity_t%pre_update')
  call solve (self, self%patch%it)
  call trace%end ()
END SUBROUTINE pre_update

!===============================================================================
!> Solve for the potential, for a given time slot
!===============================================================================
SUBROUTINE solve (self, it)
  class(selfgravity_t):: self
  integer:: it
  !.............................................................................
  class(gpatch_t), pointer:: patch
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, phi
  real:: dmin, dmax
  integer:: iter, ix ,iy, iz
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. selfgravity) return
  patch => self%patch
  call trace%begin ('selfgravity_t%solve', itimer=itimer)
  d   => patch%mem(:,:,:,patch%idx%d  ,it,1)
  phi => patch%mem(:,:,:,patch%idx%phi,it,1)
  if (self%d0==-1.0) then
    self%d0 = total_mass/total_volume
    if (verbose > 0) then
      write(io_unit%log,'(a,i6,1p,2e11.2)') &
       'selfgravity_t%pre_update: id, self%d0, faver(d) = ', &
        patch%id, self%d0, patch%faver(d)
    end if
  end if
  !.............................................................................
  dmin = patch%fminval (d)
  dmax = patch%fmaxval (d)
  if (dmin==0.0) then
    print *, patch%id, 'WARNING: selfgravity_t%gravity, dmin, dmax, it =', &
      dmin, dmax, patch%new
  end if
  !.............................................................................
  call self%poisson%update (patch, phi, d, self%d0)
  if (verbose > 1) then
    call patch%stats (d  , '  d', io_unit%log)
    call patch%stats (phi, 'phi', io_unit%log)
  end if
  !-----------------------------------------------------------------------------
  ! Compute the force per unit volume (except for ZEUS solvers, which take the
  ! gradient of the potential internally).
  !-----------------------------------------------------------------------------
  if (patch%kind(1:4) /= 'zeus') then
    call scratch%allocate (scratch%force_per_unit_mass, [patch%gn, 3])
    scratch%force_per_unit_mass = -grad (patch%mesh, phi)
  endif
  !-----------------------------------------------------------------------------
  ! compute estimate of dphi/dt
  !-----------------------------------------------------------------------------
  call time_derivs (patch)
  call trace%end (itimer)
END SUBROUTINE solve

!===============================================================================
!> Having advanced the MHD state to the new time, we make a preliminary
!> call to the Poisson solver, based on extrapolation of the nbor potentials.
!> We then provide this preliminary potential to nbors, until we have a
!> better offer.  The poisson_time controls if an update affects the self%time
!> or the poisson_time; the two are updated alternatingly
!===============================================================================
SUBROUTINE post_update (self)
  class(selfgravity_t):: self
  !.............................................................................
  class(gpatch_t), pointer:: patch
  real(kind=KindScalarVar), pointer, dimension(:,:,:,:):: phi, dphi
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. selfgravity) return
  call trace%begin ('selfgravity_t%post_update', itimer=itimer)
  if (selfgravity) then
    patch => self%patch
    if (self%poisson_time == patch%time) then
      if (verbose > 0) then
        write (io_unit%log,'(a,i6,1p,e12.4)') &
         'selfgravity_t%post_update: id, patch%time =', patch%id, patch%time
      end if
    else
      !-------------------------------------------------------------------------
      ! Extrapolate phi forward in time, if not done in the task update
      !-------------------------------------------------------------------------
      if (prediction .and. patch%nw==1) then
        phi  => patch%mem(:,:,:,patch%idx%phi ,:,1)
        dphi => patch%mem(:,:,:,patch%idx%dphi,:,1)
        call timestep%update_single (patch%id, patch%iit, patch%dt, phi, dphi, &
          patch%mesh, patch%lock, order)
      end if
      !-------------------------------------------------------------------------
      ! Advance the time to time+dtime, to prepare for post_update of phi
      !-------------------------------------------------------------------------
      call patch%rotate
      !-------------------------------------------------------------------------
      ! Fill guard zones with time extrapolated values and solve
      !-------------------------------------------------------------------------
      patch%allow_extrapolate = .true.
      call patch%dnload (only=patch%idx%d)
      call patch%dnload (only=patch%idx%phi)
      patch%allow_extrapolate = .false.
      call solve (self, patch%it)
      self%poisson_time = patch%time
      if (verbose > 0) then
        write (io_unit%log,'(a,i6,1p,e12.4)') &
         'selfgravity_t%post_update: id, poisson_time =', &
         patch%id, self%poisson_time
      end if
    end if
  end if
  call trace%end(itimer)
END SUBROUTINE post_update

!===============================================================================
!> Compute an approximate, conservative estimate of dphi/dt, which will cause
!> the timestep routine to make a good guess at the next phi, at the same time
!> as it is updating the other variables.  The 1st order prediction is in
!> practice better at reducing the number of iterations than the 2nd order.
!===============================================================================
SUBROUTINE time_derivs (self)
  class(patch_t):: self
  integer:: it1, it2, it3, nt, l(3), u(3)
  real(8):: t1, t2, t3
  integer, save:: itimer=0
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: dphidt
  !-----------------------------------------------------------------------------
  call trace%begin ('selfgravity_t%time_derivs', itimer=itimer)
  if (self%nw == 2) then
    dphidt => self%mem(:,:,:,self%idx%phi,self%new,2)
  else
    dphidt => self%mem(:,:,:,self%idx%dphi,self%new,1)
  end if
  dphidt = 0.0

  if (prediction) then
    nt = self%nt
    it1 = self%iit(nt-1)
    it2 = self%iit(nt-2)
    it3 = self%iit(nt-3)
    t1 = self%t(it1)
    t2 = self%t(it2)
    t3 = self%t(it3)
    if (it2/=it1 .and. t1/=t2) then
      l = self%mesh%li
      u = self%mesh%ui
      if (order==2 .and. t2/=t3) then
        dphidt(l(1):u(1),l(2):u(2),l(3):u(3)) = &
        (self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%phi,it1,1) - &
        self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%phi,it2,1))/(t1-t2)*1.5 &
       -(self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%phi,it2,1) - &
        self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%phi,it3,1))/(t2-t3)*0.5
      else
        dphidt(l(1):u(1),l(2):u(2),l(3):u(3)) = &
          (self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%phi,it1,1) - &
          self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%phi,it2,1))/(t1-t2)
      end if
      !.........................................................................
      if (verbose > 2) then
        write(io_unit%log,'(a,i6,1p,2e12.3,i4)') &
         'selfgravity_t%time_derivs: id, dphidt min,max, nv =', &
          self%id, self%fminval(dphidt), self%fmaxval(dphidt), size(self%mem,4)
      end if
      !.........................................................................
    end if
  else
    dphidt = 0.0
  end if
  call trace%end (itimer)
END SUBROUTINE time_derivs

!===============================================================================
!> This functions checks if 'target' is too much ahead of the nbor 'source'.
!> If it returns .true. it means the target task is not ready to be updated.
!===============================================================================
LOGICAL FUNCTION ahead_of (self, target, source)
  class(selfgravity_t):: self
  class(patch_t):: target, source
  !-----------------------------------------------------------------------------
  ! Level differences larger than 1 should never happen!
  !-----------------------------------------------------------------------------
  if (abs(target%level-source%level)>1) then
    ahead_of = .false.
    print *,'WARNING: level difference larger than one:', &
      target%id, target%level, source%id, source%level
  !-----------------------------------------------------------------------------
  ! If wee need to solve the Poisson equation, we ignore nbors with higher ot
  ! the same resolution, and require nbords with lower resolution to be in front
  ! with their Poisson solutions.
  !-----------------------------------------------------------------------------
  else if (self%poisson_time < target%time) then
    if (target%level <= source%level) then
      ahead_of = .false.
    else
      !ahead_of = target%time > source%poisson_time + source%grace*source%dtime
      ahead_of = target%time > source%time + source%grace*source%dtime
    end if
    !...........................................................................
    if (ahead_of .and. verbose > 1) then
      print 1, target%id, target%level, 'at', target%time, &
       'cannot update because', source%id, source%level, &
       'has too old Poisson time', self%poisson_time + source%grace*source%dtime
      1 format("selfgravity_t%ahead_of:",i6,i3,2x,a,g15.6,2x,a,i6,i3,2x,a,g15.6)
    end if
  !-----------------------------------------------------------------------------
  ! If we are going to update the dynamic time, we need the nbor to be in front.
  !-----------------------------------------------------------------------------
  else
    ahead_of = target%time > source%time + source%grace*source%dtime
    !...........................................................................
    if (ahead_of .and. verbose > 1) then
      print 1, target%id, target%level, 'at', target%time, &
       'cannot update because', source%id, source%level, &
       'has too old dynamic time', source%time + source%grace*source%dtime
    end if
    !...........................................................................
  end if
END FUNCTION ahead_of

END MODULE selfgravity_mod
