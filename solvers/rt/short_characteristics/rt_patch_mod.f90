!> ============================================================================
!> Module to set up angles and bins for an RT patch with only absorption.
!> This means that the only extra memory with time slot that is needed is
!> one slot for the net radiation heating, \Sum_bins \Sum_angles Q_bin_angle
!> ============================================================================
MODULE rt_patch_mod
  USE io_mod
  USE trace_mod
  USE timer_mod
  USE index_mod
  USE radau_mod
  USE rt_sc_mod
  USE solver_mod
  USE task_mod
  USE link_mod
  USE eos_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! RT patch data type, extended from MHD, with extra RT memory and parameters
  !-----------------------------------------------------------------------------
  type, public, extends(solver_t):: rt_patch_t
    logical                    :: needs_nbors=.true.
    integer                    :: n_omega=0
    integer                    :: n_bin=0
    real(8)                    :: rt_time=0.0
    real(8)                    :: rt_grace_time=0.0
    real(8), pointer           :: rt_t(:)
    type(rt_sc_t), pointer     :: rt_sc(:)
    real, pointer              :: vars(:,:,:,:,:)
  contains
    procedure:: init
    procedure:: update 
    procedure:: upstream_is_ok
    procedure:: init_nbors
    procedure:: eos_volume
  end type
  integer:: verbose=1
CONTAINS

!> ============================================================================
!> Read parameters
!> ============================================================================
SUBROUTINE init (self)
  class (rt_patch_t):: self
  integer, save:: n_mu=3, n_phi=4, n_bin=4
  integer:: i_mu, i_phi, i_omega, n_omega, iostat
  real, dimension(:), allocatable:: mu, w_mu, phi
  logical, save:: first_time=.true.
  namelist /rt_patch_params/ n_mu, n_phi, n_bin, verbose
  !-----------------------------------------------------------------------------
  ! Read parameters
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_patch_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, rt_patch_params, iostat=iostat)
    if (io%master) write (io%output, rt_patch_params)
  end if
  !$omp end critical (input_cr)
  if (.not. io%master) verbose=0
  n_omega      = 1 + n_phi*(n_mu-1)
  self%n_omega = n_omega
  self%n_bin   = n_bin
  !-----------------------------------------------------------------------------
  ! Assign angles and integration weight and allocate the required RT solvers
  ! for all angles
  !-----------------------------------------------------------------------------
  allocate (mu(n_mu),w_mu(n_mu), phi(n_phi))
  allocate (self%rt_sc(self%n_omega))
  call radau%init (n_mu, mu, w_mu)
  i_omega = 0
  do i_phi=1,n_phi
    do i_mu=1,n_mu-1
      i_omega = i_omega+1
      self%rt_sc(i_omega)%id = self%id
      self%rt_sc(i_omega)%omega%id = i_omega
      self%rt_sc(i_omega)%omega%mu  = mu(i_mu)
      self%rt_sc(i_omega)%omega%w   = w_mu(i_mu)
      self%rt_sc(i_omega)%omega%phi = phi(i_phi)
    end do
  end do
  i_omega = i_omega+1
  self%rt_sc(i_omega)%id = self%id
  self%rt_sc(i_omega)%omega%id = i_omega
  self%rt_sc(i_omega)%omega%mu  = mu(n_mu)
  self%rt_sc(i_omega)%omega%w   = w_mu(n_mu)
  self%rt_sc(i_omega)%omega%phi = 0.0
  !-----------------------------------------------------------------------------
  ! Punt to the next level, initializing patch memory and meshes
  !-----------------------------------------------------------------------------
  self%mhd = .true.
  self%nv  = 10
  self%nt  = 5
  allocate (self%rt_t(self%nt))
  call self%solver_t%init
  !-----------------------------------------------------------------------------
  ! Initialize variable indices, using slot 9 for temperature and slot 10 for Q
  !-----------------------------------------------------------------------------
  call self%idx%init (5, self%mhd)
  self%idx%tt = 9
  self%idx%qr = 10
  call self%idx%copy
  !-----------------------------------------------------------------------------
  ! Initialize the directional solvers
  !-----------------------------------------------------------------------------
  allocate (self%vars(self%gn(1),self%gn(2),self%gn(3),self%n_bin,3))
  call io%bits_mem(storage_size(self%vars),product(self%gn)*3)
  do i_omega=1,self%n_omega
    self%rt_sc(i_omega)%verbose = verbose
    self%rt_sc(i_omega)%vars => self%vars
    self%rt_sc(i_omega)%n_bin = self%n_bin
    self%rt_sc(i_omega)%mesh => self%mesh
    call self%rt_sc(i_omega)%init
  end do
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Run through all angular solvers, and attempt to update them.  Solvers that
!> lack upstream info (based on internal criteria) will do nothing, while 
!> solvers that have enough information will be updated.  The angular solvers
!> are so lightweight (taking of the order 100-200 ns) that it would be too
!> much overhead to consider them separate tasks.
!>
!> It would be easy to spawn independent OMP tasks for the ones that are
!> ready, and then wait on them, but for now we avoid complications.
!>
!> The update step computes the current RT patch time, based on the current
!> time in the attached RT solvers for each angle
!> =============================================================================
SUBROUTINE update (self)
  class (rt_patch_t):: self
  integer:: i_omega, dir
  real, dimension(:,:,:), pointer:: q
  !-----------------------------------------------------------------------------
  ! Make sure nbor lists exit
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_patch_t%update')
  if (self%needs_nbors) then
    call self%init_nbors
    self%needs_nbors = .false.
  end if
  !-----------------------------------------------------------------------------
  ! Compute the opacity and source function for all bins
  !-----------------------------------------------------------------------------
  call self%eos_volume
  !-----------------------------------------------------------------------------
  ! Compute the oldest time
  !-----------------------------------------------------------------------------
  q => self%mem(:,:,:,self%idx%qr,self%it,1)
  q = 0.0
  do i_omega=1,self%n_omega
    do dir=-1,1,2
      if (self%upstream_is_ok(i_omega,dir)) then
        self%rt_sc(i_omega)%vars => self%vars
        call self%rt_sc(i_omega)%update(dir)
      end if
    end do
    !---------------------------------------------------------------------------
    ! Sum the net heating over angle
    !---------------------------------------------------------------------------
    q = q + self%rt_sc(i_omega)%q*self%rt_sc(i_omega)%omega%w
  end do
  !----------------------------------------------------------------------------
  ! If the RT solution is far enough advanced, ask the MHD solver to update.
  ! If this update does not happen, the task time is left unchanged, and the
  ! task will again get a chance to update when the scheduler calls it next
  ! time.
  !----------------------------------------------------------------------------
  if (self%rt_time + self%rt_grace_time > self%time) then
    call self%solver_t%update
  end if
  call trace%end
END SUBROUTINE update

!> ============================================================================
!> Check upstream relations. Upstream tasks are added to the relevant nbor list
!> ============================================================================
SUBROUTINE init_nbors (self)
  class (rt_patch_t):: self
  integer:: i_omega, n_fw, n_rv
  type(link_t), pointer:: nbor, fw, fw_prev, rv, rv_prev
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_patch_t%init_nbors')
  do i_omega=1,self%n_omega
    nbor => self%link%nbor
    allocate (fw, source=nbor)
    allocate (rv, source=nbor)
    self%rt_sc(i_omega)%fw_nbor => fw
    self%rt_sc(i_omega)%rv_nbor => rv
    fw_prev => fw
    rv_prev => rv
    n_fw = 0
    n_rv = 0
    do while (associated(nbor))
      allocate (fw, source=nbor)
      allocate (rv, source=nbor)
      if (self%rt_sc(i_omega)%needs(nbor%task,+1)) then
        fw_prev%next => fw
        fw_prev => fw
        n_fw = n_fw+1
      end if
      if (self%rt_sc(i_omega)%needs(nbor%task,-1)) then
        rv_prev%next => rv
        rv_prev => rv
        n_rv = n_rv+1
      end if
      nbor => nbor%next
    end do
  end do
  if (verbose > 0) &
    print *, 'rt_patch_t%init_nbors: id, n_fw, n_rv =', self%id, n_fw, n_rv
  call trace%end
END SUBROUTINE init_nbors

!> ============================================================================
!> Check upstream times.  Upstream tasks are marked 'needed' in the nbor list
!> ============================================================================
FUNCTION upstream_is_ok (self, i_omega, dir)
  logical:: upstream_is_ok
  class (rt_patch_t):: self
  integer:: i_omega, dir
  type(link_t), pointer:: link, nbor, next
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_patch_t%upstream_is_ok')
  upstream_is_ok = .true.
  if (dir>0) then
    nbor => self%rt_sc(i_omega)%fw_nbor
  else
    nbor => self%rt_sc(i_omega)%rv_nbor
  end if
  do while (associated(nbor))
    task => nbor%task
    select type (task)
    type is (rt_patch_t)
    if (nbor%needed .and. task%rt_time < self%time - self%grace*self%dtime) then
      upstream_is_ok = .false.
      exit
    end if
    end select
    nbor => nbor%next
  end do
  call trace%end
END FUNCTION upstream_is_ok

!> ============================================================================
!> Compute the EOS quantities for the volume.  The 1,2 slots of the last index
!> correspond to opacity and source function, respectively.
!> ============================================================================
SUBROUTINE eos_volume (self)
  class (rt_patch_t):: self
  integer, save:: itimer=0
  real, dimension(:,:,:), pointer:: d, e
  real, dimension(:,:,:,:), pointer:: rk, src
  integer:: m(4), ib, i3
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_sc_t%eos_volume', itimer=itimer)
  d   => self%mem(:,:,:,self%idx%d,self%it,1)
  e   => self%mem(:,:,:,self%idx%e,self%it,1)
  rk  => self%vars(:,:,:,:,1)
  src => self%vars(:,:,:,:,2)
  m = shape(rk)
  call eos%lookup_table (m, d=d, e=e, rk=rk, src=src)
  call trace%end (itimer)
END SUBROUTINE eos_volume

END MODULE rt_patch_mod
