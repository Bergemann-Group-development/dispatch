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
  USE halo_mod
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
  USE math_mod
  USE guards_mod
  implicit none
  private
  type, public:: selfgravity_t
    type(poisson_t)  :: poisson
    logical:: on=.false.
    !logical:: selfgravity=.false.
    real(8):: fourPiG
    real(8):: d0
    real(8):: poisson_time=-1.0
    class(gpatch_t), pointer         :: patch      => null()
    real(var_k), pointer, contiguous ::   d(:,:,:) => null()
    real(var_k), pointer, contiguous :: phi(:,:,:) => null()
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: pre_update
    procedure:: post_update
    procedure, private:: download
    procedure:: solve
    procedure:: stats
  end type
  real(8), save :: total_mass=0d0, total_volume=0d0
  real(8), save :: restart_time = 0d0
  logical, save :: on=.false., prediction=.true.
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
  namelist /selfgravity_params/ on, prediction, order, d0, verbose
  !-----------------------------------------------------------------------------
  call trace%begin ('selfgravity_t%pre_init')
  self%patch => patch
  !
  call io%lock%set ('selfgravity_t%pre_init')
  if (first_time) then
    rewind (io%input)
    read (io%input, selfgravity_params, iostat=iostat)
    call io%namelist_warning('selfgravity_params', iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, selfgravity_params)
    end if
    write (io%namelist, selfgravity_params)
    first_time = .false.
  end if
  call io%lock%unset ('selfgravity_t%pre_init')
  !.............................................................................
  self%d0 = d0
  if (on) then
    self%on = on
    patch%idx%phi = 0
    if (verbose > 0) &
      write (stderr,*) 'id,nv,mv=', patch%id,patch%nv,patch%mv
    call patch%idx%request_index (patch%idx%phi, patch%mv)
    !---------------------------------------------------------------------------
    ! We need to keep previous time slots of the potential, in order to be able
    ! to do time extrapolation in the post_update step.  However, the potential
    ! should not be included among the %nv variables, since it needs special
    ! treatment of the guard zones; cf. self%download() below.
    !---------------------------------------------------------------------------
    call self%poisson%init
  else
    patch%idx%phi = -1
  end if
  self%fourPiG = math%pi4*scaling%grav
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
  end if
  !$omp end critical (poisson_cr)
  if (on) then
    !---------------------------------------------------------------------------
    ! Sum up the total mass and total volume, to get the average density d0.
    ! FIXME: This needs to be modified for MPI
    !---------------------------------------------------------------------------
    if (patch%time == 0.0) then
      self%d => patch%mem(:,:,:,patch%idx%d,patch%it,1)
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
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> The pre-update step starts with only density downloaded (done in halo_t),
!> and needs to download the potential in the guard zones (the interior already
!> has values consistent on the small scales).  This is enough to get the
!> desired solution, but the number of iterations can possibly be reduced by
!> adding the large scale change from the parent patch.  This would be cheaper
!> if the prolonged large scale solution was kept in a separate series of time
!> slots.
!===============================================================================
SUBROUTINE pre_update (self)
  class(selfgravity_t):: self
  class(gpatch_t), pointer:: patch
  integer, dimension(3):: l, u
  integer:: iy, iz
 !------------------------------------------------------------------------------
  if (.not. on) return
  call trace%begin ('selfgravity_t%pre_update')
  patch => self%patch
  self%d   => patch%mem(:,:,:,patch%idx%d  ,patch%it,1)
  self%phi => patch%mem(:,:,:,patch%idx%phi,patch%it,1)
  !.............................................................................
  if (verbose > 1) then
    call yz_indices (patch, iy, iz)
    write (stderr,*) 'pre0 phi', patch%id, patch%it
    write (stderr,'((18f8.4))') self%phi (:,iy,iz)
  end if
  call self%download ()
  if (verbose > 1) then
    call yz_indices (patch, iy, iz)
    write (stderr,*) 'pre1 phi', patch%id, patch%it
    write (stderr,'((18f8.4))') self%phi (:,iy,iz)
  end if
  !*****************************************************************************
  call self%solve ()
  !*****************************************************************************
  if (verbose > 0) then
    call yz_indices (patch, iy, iz)
    write (stderr,*) 'pre2 phi', patch%id, patch%it
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
  call trace%end ()
END SUBROUTINE pre_update

!===============================================================================
!> Get potential guard zone values from the same and lower levels, in reverse
!> order.  This ensures that changes from larger scales are efficiently
!> propagated to smaller scales.  To speed up convergence, prologong from
!> overlapping larger scales
!===============================================================================
SUBROUTINE download (self)
  class(selfgravity_t):: self
  class(patch_t), pointer:: tgt
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor, reverse, new, next
  class(task_t), pointer:: src
  !-----------------------------------------------------------------------------
  ! Constructe a reversed nbor list by appending to a new pointer, keeping only
  ! nbors at the same or lower levels
  !-----------------------------------------------------------------------------
  call trace%begin ('selfgravity_t%download')
  tgt => self%patch
  nullify (reverse)
  nbor => tgt%link%nbor
  do while (associated(nbor))
    if (nbor%task%level <= tgt%level) then
      !-------------------------------------------------------------------------
      ! In the special case where tgt%istep=0, use the parent to load phi
      !-------------------------------------------------------------------------
      if (tgt%istep == 0 .and. nbor%task%id == tgt%parent%id) then
        link => self%patch%link
        call link%download (link%nbor, only=tgt%idx%phi, &
                            all_tgt=.true., all_src=.true.)
      end if
      !-------------------------------------------------------------------------
      ! Prepend each task, reversing the order
      !-------------------------------------------------------------------------
      allocate (new, source=nbor)
      new%next => reverse
      reverse => new
    end if
    nbor => nbor%next
  end do
  !-----------------------------------------------------------------------------
  ! Download guard zone values of the potential from the reverse list
  !-----------------------------------------------------------------------------
  link => self%patch%link
  select type (link)
  class is (halo_t)
    call self%patch%link%download (reverse, only=tgt%idx%phi)
  end select
  !-----------------------------------------------------------------------------
  ! Deallocate reversed list
  !-----------------------------------------------------------------------------
  nbor => reverse
  do while (associated(nbor))
    next => nbor%next
    deallocate (nbor)
    nbor => next
  end do
  call trace%end ()
END SUBROUTINE download

!===============================================================================
!> The post-update step updates the potential to a preliminart value at time
!> t+dt, by time-extrapolating the potential to the 'new' timeslot.
!===============================================================================
SUBROUTINE post_update (self)
  class(selfgravity_t):: self
  !.............................................................................
  class(gpatch_t), pointer:: patch
  real(var_k), dimension(:,:,:), pointer:: phi_now, phi_prv, d_now, d_prv 
  real(var_k), dimension(:,:,:,:), pointer:: phi_t, dphi_t
  integer:: it_now, it_prv  
  integer:: iy, iz
  real(8):: time
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. on) return
  call trace%begin ('selfgravity_t%post_update', itimer=itimer)
  patch => self%patch
  !-----------------------------------------------------------------------------
  ! To improve the potential IC, extrapolate phi forward in time, if not done
  ! in the task update
  !-----------------------------------------------------------------------------
   if (verbose > 1) then
    call yz_indices (patch, iy, iz)
    write (stderr,*) 'post0 phi', patch%id, patch%it
    self%phi => patch%mem(:,:,:,patch%idx%phi,patch%it,1)
    write (stderr,'((18f8.4))') self%phi (:,iy,iz)
  end if
  time = patch%time + patch%dtime
  call patch%time_extrapolate (patch%idx%phi, order, time)
  self%d   => patch%mem(:,:,:,patch%idx%d  ,patch%new,1)
  self%phi => patch%mem(:,:,:,patch%idx%phi,patch%new,1)
  if (verbose > 1) then
    write (stderr,*) 'post1 phi', patch%id, patch%new
    write (stderr,'((18f8.4))') self%phi (:,iy,iz)
  end if
  call self%solve ()
  if (verbose > 1) then
    write (stderr,*) 'post2 phi', patch%id, patch%new
    write (stderr,'((18f8.4))') self%phi (:,iy,iz)
  end if
  !-----------------------------------------------------------------------------
  call trace%end(itimer)
END SUBROUTINE post_update

!===============================================================================
!> Solve the Poisson equation with given src function (the density self%d
!> is only used to gauge the residual error)
!===============================================================================
SUBROUTINE solve (self)
  class(selfgravity_t):: self
  real, pointer, contiguous:: src(:,:,:)
  integer:: m(3), l(3), u(3)
  !.............................................................................
  m = self%patch%gn(1:3)
  l = self%patch%mesh%li
  u = self%patch%mesh%ui
  allocate (src(m(1),m(2),m(3)))
  src = self%fourPiG*self%d
  call self%poisson%update (self%patch, src, self%phi, self%d, l, u)
  deallocate (src)
END SUBROUTINE solve

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
