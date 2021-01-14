!===============================================================================
!> Add selfgravity to any solver.
!>
!> Use residuals r = source - laplace(phi) to generate delta-Phi values.
!> Corrections computed by a patch should then by pushed back to the higher
!> level child patches, recursively.  The point is to use the lower level
!> corrections to bring the longer wavelength parts of the residuals at 
!> higher levels down, before the local solutions are iterated.
!>
!> Residuals are computed at the current level for internal cells, and are then
!> overloaded by restriction from child tasks, if any.
!>
!> Increments dphi to the potential are then downloaded from the parent into
!> the interior and guard zones, and are overloaded by better values from
!> tasks on the same level.
!>
!> Each patch computes, from the given residuals rphi, correctioins dphi in
!> its interior, while also making the corrections available for download.
!>
!> pre_update step:
!> ----------------
!> The pre_update step has accurate BCs for the potential, and obtains the
!> best possible values for the residual from the %downloads procedure.
!> The Poisson solver produces %dphi values, which are applied to the %phi
!> values, and are made available for time- and space-interpolation to nbors.
!>
!> post_update step:
!> -----------------
!> The purpose of this step is to correct the potential in response to the
!> local changes in the density, while using potential BCs obtained from
!> time extrapolation in the download procedure.
!>
!> Both solutions start by computing the residual as the difference between
!> the local source function and the laplacian of the local potential.  The
!> potential BCs and ICs have been modified by nbors at the same and lower
!> levels, which help reduce longer range residuals, while the residuals
!> include contributions from higher level task, which need help to correct
!> their potentials.
!>
!> The scheme calls for three variable slots: %phi, %dphi, and %rphi, but
!> only %dphi and %rphi need to exist as time slices
!===============================================================================
MODULE selfgravity_mod
  USE io_mod
  USE io_unit_mod
  USE os_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE gpatch_mod
  USE link_mod
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
  USE shared_mod
  USE scaling_mod
  !
  USE download_mod
  USE guard_zones_mod
  implicit none
  private
  type, public:: selfgravity_t
    type(initial_t)  :: initial
    type(poisson_t)  :: poisson
    logical:: selfgravity=.false.
    real(8):: d0
    real(8):: poisson_time=-1.0
    class(gpatch_t), pointer :: patch=>null()
    real(var_k), allocatable ::   src(:,:,:)
    real(var_k), allocatable ::  dphi(:,:,:)
    real(var_k), pointer     ::     d(:,:,:) => null()
    real(var_k), pointer     :: d_pre(:,:,:) => null()
    real(var_k), pointer     ::   res(:,:,:) => null()
    real(var_k), pointer     ::   phi(:,:,:) => null()
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: pre_update
    procedure:: downloads
    procedure:: uploads
    procedure:: restrict
    procedure:: prolong
    procedure:: post_update
    !procedure:: ahead_of
    procedure:: stats
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
    if (verbose > 0) &
      write (stderr,*) 'id,nv,mv=', patch%id,patch%nv,patch%mv
    call patch%idx%request_index (patch%idx%phi, patch%mv)
    if (patch%nw == 2) then
      patch%idx%rphi = patch%idx%phi
    else
      patch%idx%rphi = 0
      call patch%idx%request_index (patch%idx%rphi, patch%mv)
    end if
    !...........................................................................
    if (verbose > 0) then
      write (io_unit%log,'(a,i6,4i4)') &
       'selfgravity_t%pre_init: id, nv, nw, rphi, dphi =', &
        patch%id, patch%nv, patch%nw, patch%idx%rphi, patch%idx%phi
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
!> Actual initialisation of the self-gravity module
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
    self%d => patch%mem(:,:,:,patch%idx%d,patch%it,1)
    if (patch%time == 0.0) then
      dm = patch%fsum(self%d)*product(patch%ds)
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
  self%poisson%id = patch%id
  !-----------------------------------------------------------------------------
  if (.not.allocated(self%dphi)) &
    allocate (self%dphi(patch%gn(1),patch%gn(2),patch%gn(3)))
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
  class(gpatch_t), pointer:: patch
  integer, dimension(3):: l, u
  integer:: iy, iz
 !------------------------------------------------------------------------------
  if (.not. selfgravity) return
  call trace%begin ('selfgravity_t%pre_update')
  patch => self%patch
  self%d   => patch%mem(:,:,:,patch%idx%d  ,patch%it,1)
  self%phi => patch%mem(:,:,:,patch%idx%phi,patch%it,1)
  self%d_pre => self%d
  allocate (self%src(patch%gn(1),patch%gn(2),patch%gn(3)))
  !-----------------------------------------------------------------------------
  ! If the patch is the root patch, compute the mean density
  !-----------------------------------------------------------------------------
  call sentinel (self, 'pre beg')
  if (patch%is_periodic()) then
    call patch%make_periodic (self%phi)
    l = patch%li
    u = patch%ui
    shared%mean_density = &
      sum(self%d(l(1):u(1),l(2):u(2),l(3):u(3)))/product(u-l+1)
  end if
  !-----------------------------------------------------------------------------
  ! Initialize potential first time around, then upload potentials from lower
  ! levels as boundary conditions, compute residuals locally, amd download
  ! residuals from overlying higher levels
  !*****************************************************************************
  if (patch%istep == 0) then
    self%phi = 0.0
  end if
  call self%uploads
  call self%downloads
  !*****************************************************************************
  if (verbose > 1) then
    call yz_indices (patch, iy, iz)
    write (stderr,*) 'pre1 src' , patch%id
    write (stderr,'((18f8.4))') self%src (:,iy,iz)
    write (stderr,*) 'pre1 res' , patch%id
    write (stderr,'((18f8.4))') self%res (:,iy,iz)
    write (stderr,*) 'pre1 phi', patch%id
    write (stderr,'((18f8.4))') self%phi (:,iy,iz)
  end if
  !-----------------------------------------------------------------------------
  ! The solver adds the correction to phi that corresponds to the source 
  !*****************************************************************************
  call self%poisson%update (self%patch, self%res, self%dphi, self%d)
  self%phi = self%phi + self%dphi
  !*****************************************************************************
  call sentinel (self, 'pre upd')
  if (verbose > 0) then
    call yz_indices (patch, iy, iz)
    write (stderr,*) 'pre2 res' , patch%id
    write (stderr,'((18f8.4))') self%res(:,iy,iz)
    write (stderr,*) 'pre2 src', patch%id
    write (stderr,'((18f8.4))') self%src(:,iy,iz)
    write (stderr,*) 'pre2 lap', patch%id
    write (stderr,'((18f8.4))') self%src(:,iy,iz) - self%res(:,iy,iz)
    write (stderr,*) 'pre2 phi', patch%id
    write (stderr,'((18f8.4))') self%phi(:,iy,iz)
  end if
  !-----------------------------------------------------------------------------
  ! Compute the force per unit volume (except for ZEUS solvers, which take the
  ! gradient of the potential internally).
  !-----------------------------------------------------------------------------
  if (patch%kind(1:4) /= 'zeus') then
    call scratch%allocate (scratch%force_per_unit_mass, [patch%gn, 3])
    scratch%force_per_unit_mass = -grad (patch%mesh, self%phi)
  endif
  !-----------------------------------------------------------------------------
  ! Compute estimate of dphi/dt
  !-----------------------------------------------------------------------------
!  call time_derivs (patch)
  call sentinel (self, 'pre end')
  call trace%end ()
END SUBROUTINE pre_update

!===============================================================================
!> Update the potential to the new time, using time extrapolated potential BCs
!> and the increment of internal densities; this picks up most of the changes
!> of interest to nbors looking for BCs.  The potential is sub-sequently
!> corrected in the next pre_update, when the potenial BCs are better known.
!===============================================================================
SUBROUTINE post_update (self)
  class(selfgravity_t):: self
  !.............................................................................
  class(gpatch_t), pointer:: patch
  real(var_k), dimension(:,:,:), pointer:: phi_now, phi_prv, d_now, d_prv 
  real(var_k), dimension(:,:,:,:), pointer:: phi_t, dphi_t
  integer:: it_now, it_prv  
  integer:: iy, iz
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. selfgravity) return
  call sentinel (self, 'post beg')
  call trace%begin ('selfgravity_t%post_update', itimer=itimer)
  patch => self%patch
  !-----------------------------------------------------------------------------
  ! To improve the potential IC, extrapolate phi forward in time, if not done
  ! in the task update
  !-----------------------------------------------------------------------------
  if (prediction .and. patch%nw==1) then
    phi_t  => patch%mem(:,:,:,patch%idx%phi ,:,1)
    dphi_t => patch%mem(:,:,:,patch%idx%dphi,:,1)
    call timestep%update_single (patch%id, patch%iit, patch%dt, phi_t, dphi_t, &
      patch%mesh, patch%lock, order)
  end if
  !-----------------------------------------------------------------------------
  ! Advance the patch%time to time+dtime, to prepare for post_update of phi
  !-----------------------------------------------------------------------------
  patch%mem(:,:,:,patch%idx%phi ,patch%new,1) = &
  patch%mem(:,:,:,patch%idx%phi ,patch%it ,1)
  patch%mem(:,:,:,patch%idx%rphi,patch%new,1) = &
  patch%mem(:,:,:,patch%idx%rphi,patch%it ,1)
  call patch%rotate
  self%d   => patch%mem(:,:,:,patch%idx%d  ,patch%it,1)
  self%phi => patch%mem(:,:,:,patch%idx%phi,patch%it,1)
  !-----------------------------------------------------------------------------
  ! Compute residuals locally, and download overlying higher levels and
  ! pick up the potential from lower level
  !-----------------------------------------------------------------------------
  if (verbose > 1) then
    call yz_indices (patch, iy, iz)
    write (stderr,*) 'post0 phi', patch%id
    write (stderr,'((18f8.4))') self%phi (:,iy,iz)
  end if
  !*****************************************************************************
  patch%allow_extrapolate = .true.
  call self%uploads
  patch%allow_extrapolate = .false.
  !*****************************************************************************
  call sentinel (self, 'post dnl')
  if (verbose > 1) then
    call yz_indices (patch, iy, iz)
    write (stderr,*) 'post1 src' , patch%id
    write (stderr,'((18f8.4))') self%src (:,iy,iz)
    write (stderr,*) 'post1 res' , patch%id
    write (stderr,'((18f8.4))') self%res (:,iy,iz)
    write (stderr,*) 'post1 phi', patch%id
    write (stderr,'((18f8.4))') self%phi (:,iy,iz)
  end if
  !-----------------------------------------------------------------------------
  ! The solver incremnts the potential, based on the new BC values and the
  ! residual, now corrected for the local changes of the density
  !*****************************************************************************
  call poisson%laplace_2nd (patch%ds, self%phi, self%src, self%res)
  self%res = self%res + scaling%grav*(self%d-self%d_pre)
  call self%poisson%update (self%patch, self%res, self%dphi, self%d)
  self%phi = self%phi + self%dphi
  !*****************************************************************************
  call sentinel (self, 'post upd')
  if (verbose > 0) then
    call yz_indices (patch, iy, iz)
    write (stderr,*) 'post2 res' , patch%id
    write (stderr,'((18f8.4))') self%res (:,iy,iz)
    write (stderr,*) 'post2 dphi', patch%id
    write (stderr,'((18f8.4))') self%dphi(:,iy,iz)
    write (stderr,*) 'post2 phi', patch%id
    write (stderr,'((18f8.4))') self%phi (:,iy,iz)
  end if
  !*****************************************************************************
  deallocate (self%src)
  call sentinel (self, 'post end')
  call trace%end(itimer)
END SUBROUTINE post_update

!===============================================================================
!> Downlaod residuals (only) from higher level patches that overlap with the
!> internal cells of the patch
!===============================================================================
SUBROUTINE downloads (self)
  class(selfgravity_t):: self
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: patch, nbpatch
  real(var_k), pointer:: res(:,:,:)
  integer:: iy, iz
  !-----------------------------------------------------------------------------
  ! Download values from the L+1 tasks contained in this task
  !-----------------------------------------------------------------------------
  patch => self%patch
  select type (patch)
  class is (patch_t)
  self%d   => patch%mem(:,:,:,patch%idx%d   ,patch%it,1)
  self%phi => patch%mem(:,:,:,patch%idx%phi ,patch%it,1)
  self%res => patch%mem(:,:,:,patch%idx%rphi,patch%it,1)
  self%src = scaling%grav*(self%d-shared%mean_density)
  call sentinel (self, 'downloads')
  call poisson%laplace_2nd (patch%ds, self%phi, self%src, self%res)
  nbor => patch%link%nbor
  do while (associated(nbor))
    nbpatch => nbor%task
    select type (nbpatch)
    class is (patch_t)
    if (nbpatch%level == patch%level+1 ) then
      if (all(abs(patch%distance (nbpatch)) < 0.5*patch%size)) then
        call restrict (self, nbpatch, patch, patch%idx%rphi)
!        write (stdout,'(2(a,i6,":",i2.2))') 'downloads', &
!          nbpatch%id, nbpatch%level, '  to', patch%id, patch%level
!        call yz_indices (patch, iy, iz)
!        write (stderr,*) 'res', patch%id, iy, iz
!        write (stderr,'((18f8.4))') self%res(:,iy,iz)
      end if
    end if
    end select
    nbor => nbor%next
  end do
  end select
END SUBROUTINE downloads

!===============================================================================
!> Uplaod potentials corrections (only) from higher level patches that overlap
!> with the internal cells of the patch
!===============================================================================
SUBROUTINE uploads (self)
  class(selfgravity_t):: self
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: patch, nbpatch
  integer:: iy, iz, iit(9)
  !-----------------------------------------------------------------------------
  ! Upload values from the L-1 tasks on the nbor list
  !-----------------------------------------------------------------------------
  patch => self%patch
  select type (patch)
  class is (patch_t)
  self%d   => patch%mem(:,:,:,patch%idx%d  ,patch%it,1)
  self%phi => patch%mem(:,:,:,patch%idx%phi,patch%it,1)
  nbor => patch%link%nbor
  do while (associated(nbor))
    nbpatch => nbor%task
    select type (nbpatch)
    class is (patch_t)
      if (nbpatch%level == patch%level-1) then
        if (all(abs(patch%distance (nbpatch)) < 0.5*nbpatch%size)) then
          call self%prolong (nbpatch, patch, patch%idx%phi)
        end if
      end if
    end select
    nbor => nbor%next
  end do
  end select
END SUBROUTINE uploads

!===============================================================================
!> Indices through center
!===============================================================================
SUBROUTINE yz_indices (patch, iy, iz)
  class(patch_t):: patch
  integer:: iy, iz
  if (patch%position(2) < -0.05) then
    iy = 34
  else if (patch%position(2) > 0.05) then
    iy = 3
  else
    iy = 19
  end if
  if (patch%position(3) < -0.05) then
    iz = 34
  else if (patch%position(3) > 0.05) then
    iz = 3
  else
    iz = 19
  end if
END SUBROUTINE yz_indices

!===============================================================================
!> Restrict
!===============================================================================
SUBROUTINE restrict0 (self, source, target, iv)
  class(selfgravity_t):: self
  class(patch_t), pointer:: source, target
  integer:: iv, jt(2)
  real:: pt(2)
  !-----------------------------------------------------------------------------
  call source%time_interval (target%time, jt, pt, &
                             extrap=target%allow_extrapolate)
  if (verbose > 1) then
    write (stderr,'(a,i6,2f7.3)') 'extrap: src,pt=', source%id, pt
  end if
  call download%different (target, source, jt, pt, only=iv)
END SUBROUTINE restrict0

!===============================================================================
!> Restrict
!===============================================================================
SUBROUTINE restrict (self, source, target, iv)
  class(selfgravity_t)    :: self
  class(patch_t), pointer :: source, target
  integer                 :: iv
  real, pointer           :: src(:,:,:,:), tgt(:,:,:)
  integer                 :: i1, i2, i3, j1, j2, j3, k(2), l(3), u(3)
  real                    :: pt(2), p1, p2, p3
  real(8)                 :: d(3), fl(3), fu(3)
  logical                 :: outside(3)
  class(mesh_t), pointer  :: sm, tm
  !-----------------------------------------------------------------------------
  call source%time_interval (target%time, k, pt, &
                             extrap=target%allow_extrapolate)
  src => source%mem(:,:,:,iv,:,1)
  tgt => target%mem(:,:,:,iv,target%it,1)
  if (verbose > 1) then
    write (stderr,'("source: id,t=",i6,f8.4)') source%id, source%time
    write (stderr,'(10f8.4)') source%t(source%iit)
    write (stderr,'(11i7)'  ) source%iit, k
    write (stderr,'(9f8.4)' ) src(18,18,18,source%iit)
  end if
  d = target%distance (source)
  fl = d - source%size/2d0
  fu = d + source%size/2d0
  l = ceiling(fl/target%mesh%d + target%mesh%o) + 1
  u =   floor(fu/target%mesh%d + target%mesh%o) - 1
  do i3=l(3),u(3)
    tm => target%mesh(3)
    sm => source%mesh(3)
    outside(3) = i3 < tm%li .or. i3 > tm%ui
    p3 = (tm%r(i3) + d(3))/sm%d + sm%o
    j3 = p3
    p3 = p3 - j3
    do i2=l(2),u(2)
      tm => target%mesh(2)
      sm => source%mesh(2)
      outside(2) = i2 < tm%li .or. i2 > tm%ui
      p2 = (tm%r(i2) + d(2))/sm%d + sm%o
      j2 = p2
      p2 = p2 - j2
      do i1=l(1),u(1)
        outside(1) = i1 < tm%li .or. i1 > tm%ui
        if (any(outside)) then
          tm => target%mesh(1)
          sm => source%mesh(1)
          p1 = (tm%r(i1) + d(1))/sm%d + sm%o
          j1 = p1
          p1 = p1 - j1
          tgt(i1,i2,i3) = &
              pt(1)*((1.-p3)*((1.-p2)*((1.0-p1)*src(j1  ,j2  ,j3  ,k(1))    + &
                                       (    p1)*src(j1+1,j2  ,j3  ,k(1)))   + &
                              (   p2)*((1.0-p1)*src(j1  ,j2+1,j3  ,k(1))    + &
                                       (    p1)*src(j1+1,j2+1,j3  ,k(1))))  + &
                     (   p3)*((1.-p2)*((1.0-p1)*src(j1  ,j2  ,j3+1,k(1))    + &
                                       (    p1)*src(j1+1,j2  ,j3+1,k(1)))   + &
                              (   p2)*((1.0-p1)*src(j1  ,j2+1,j3+1,k(1))    + &
                                       (    p1)*src(j1+1,j2+1,j3+1,k(1))))) + &
              pt(2)*((1.-p3)*((1.-p2)*((1.0-p1)*src(j1  ,j2  ,j3  ,k(2))    + &
                                       (    p1)*src(j1+1,j2  ,j3  ,k(2)))   + &
                              (   p2)*((1.0-p1)*src(j1  ,j2+1,j3  ,k(2))    + &
                                       (    p1)*src(j1+1,j2+1,j3  ,k(2))))  + &
                     (   p3)*((1.-p2)*((1.0-p1)*src(j1  ,j2  ,j3+1,k(2))    + &
                                       (    p1)*src(j1+1,j2  ,j3+1,k(2)))   + &
                              (   p2)*((1.0-p1)*src(j1  ,j2+1,j3+1,k(2))    + &
                                      (    p1)*src(j1+1,j2+1,j3+1,k(2)))))
        end if
      end do
    end do
  end do
END SUBROUTINE restrict

!===============================================================================
!> Prolong
!===============================================================================
SUBROUTINE prolong (self, source, target, iv)
  class(selfgravity_t)    :: self
  class(patch_t), pointer :: source, target
  integer                 :: iv
  real, pointer           :: src(:,:,:,:), tgt(:,:,:)
  integer                 :: i1, i2, i3, j1, j2, j3, k(2)
  real                    :: pt(2), p1, p2, p3
  real(8)                 :: d(3)
  logical                 :: outside(3)
  class(mesh_t), pointer  :: sm, tm
  integer                 :: tid, sid
  !-----------------------------------------------------------------------------
  call source%time_interval (target%time, k, pt, &
                             extrap=target%allow_extrapolate)
  src => source%mem(:,:,:,iv,:,1)
  tgt => target%mem(:,:,:,iv,target%it,1)
  tid = target%id
  sid = source%id
  if (verbose > 1) then
    write (stderr,'("source: id,t=",i6,f8.4)') source%id, source%time
    write (stderr,'(10f8.4)') source%t(source%iit)
    write (stderr,'(11i7)'  ) source%iit, k
    write (stderr,'(9f8.4)' ) src(18,18,18,source%iit)
  end if
  d = target%distance (source)
  do i3=target%mesh(3)%lb+1,target%mesh(3)%ub-1
    tm => target%mesh(3)
    sm => source%mesh(3)
    outside(3) = i3 < tm%li .or. i3 > tm%ui
    p3 = (tm%r(i3) + d(3))/sm%d + sm%o
    j3 = p3
    if (j3 > sm%gn-1) then
      j3 = sm%gn - 1
    end if
    p3 = p3 - j3
    do i2=target%mesh(2)%lb+1,target%mesh(2)%ub-1
      tm => target%mesh(2)
      sm => source%mesh(2)
      outside(2) = i2 < tm%li .or. i2 > tm%ui
      p2 = (tm%r(i2) + d(2))/sm%d + sm%o
      j2 = p2
      if (j2 > sm%gn-1) then
        j2 = sm%gn - 1
      end if
      p2 = p2 - j2
      do i1=target%mesh(1)%lb+1,target%mesh(1)%ub-1
        outside(1) = i1 < tm%li .or. i1 > tm%ui
        if (any(outside)) then
          tm => target%mesh(1)
          sm => source%mesh(1)
          p1 = (tm%r(i1) + d(1))/sm%d + sm%o
          j1 = p1
          if (j1 > sm%gn-1) then
            j1 = sm%gn - 1
          end if
          if (j1 < 1) then
            j1 = 1
          end if
          p1 = p1 - j1
          tgt(i1,i2,i3) = &
              pt(1)*((1.-p3)*((1.-p2)*((1.0-p1)*src(j1  ,j2  ,j3  ,k(1))    + &
                                       (    p1)*src(j1+1,j2  ,j3  ,k(1)))   + &
                              (   p2)*((1.0-p1)*src(j1  ,j2+1,j3  ,k(1))    + &
                                       (    p1)*src(j1+1,j2+1,j3  ,k(1))))  + &
                     (   p3)*((1.-p2)*((1.0-p1)*src(j1  ,j2  ,j3+1,k(1))    + &
                                       (    p1)*src(j1+1,j2  ,j3+1,k(1)))   + &
                              (   p2)*((1.0-p1)*src(j1  ,j2+1,j3+1,k(1))    + &
                                       (    p1)*src(j1+1,j2+1,j3+1,k(1))))) + &
              pt(2)*((1.-p3)*((1.-p2)*((1.0-p1)*src(j1  ,j2  ,j3  ,k(2))    + &
                                       (    p1)*src(j1+1,j2  ,j3  ,k(2)))   + &
                              (   p2)*((1.0-p1)*src(j1  ,j2+1,j3  ,k(2))    + &
                                       (    p1)*src(j1+1,j2+1,j3  ,k(2))))  + &
                     (   p3)*((1.-p2)*((1.0-p1)*src(j1  ,j2  ,j3+1,k(2))    + &
                                       (    p1)*src(j1+1,j2  ,j3+1,k(2)))   + &
                              (   p2)*((1.0-p1)*src(j1  ,j2+1,j3+1,k(2))    + &
                                      (    p1)*src(j1+1,j2+1,j3+1,k(2)))))
        end if
      end do
    end do
  end do
  if (verbose > 99) print *, sid, tid
END SUBROUTINE prolong

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

!===============================================================================
SUBROUTINE stats (self, f, label, unit)
  class(selfgravity_t)   :: self
  real, dimension(:,:,:) :: f
  character(*)           :: label
  integer                :: unit
  !-----------------------------------------------------------------------------
  if (verbose > 1) then
    call self%patch%stats (f, label, unit)
  end if
END SUBROUTINE stats

SUBROUTINE sentinel (self, label)
  class(selfgravity_t):: self
  character(*):: label
  class(patch_t), pointer:: patch => null()
  real:: phi, phi_prv=0.0
  save:: phi_prv, patch
  !-----------------------------------------------------------------------------
  if (verbose > 1) then
    write (stderr,'(i6,2x,"average density:",f10.6,2x,a)') self%patch%id, &
      self%patch%faver (self%patch%idx%d), label
  end if
  if (self%patch%id == 1) then
    patch => self%patch
    if (verbose > 0) then
      write (stderr,'(i6,2x,"average density:",f10.6,2x,a)') self%patch%id, &
        self%patch%faver (self%patch%idx%d), label
    end if
  end if
  if (associated(patch)) then
    phi = patch%mem(18,18,18,patch%idx%phi,patch%it,1)
    if (phi /= phi_prv .and. verbose > 0) then
      !write (stderr,'(a,i6,i4,f8.4,2x,a)') &
      write (stderr,*) &
       "phi changed(1): id,istep,phi=", self%patch%id, patch%istep, phi, label
      phi_prv = phi
    end if
    if (self%patch%id == 1 .and. verbose > 0) then
      phi = self%phi(18,18,18)
      if (phi /= phi_prv) then
        !write (stderr,'(a,i6,i4,f8.4,2x,a)') &
        write (stderr,*) &
        "phi changed(2): id,istep,phi=", self%patch%id, patch%istep, phi, label
        phi_prv = phi
      end if
    end if
  end if
END SUBROUTINE sentinel

END MODULE selfgravity_mod
