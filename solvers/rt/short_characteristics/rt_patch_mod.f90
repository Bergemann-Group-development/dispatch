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
  USE mhd_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! RT patch data type, extended from MHD, with extra RT memory and parameters
  !-----------------------------------------------------------------------------
  type, public, extends(mhd_t):: rt_patch_t
    real(8), allocatable:: rt_time(:)
    real(8):: rt_grace_time=0.0
    integer:: n_omega, n_bin
    type(rt_sc_t), dimension(:), allocatable:: rt_sc
  contains
    procedure:: init
    procedure:: update 
    procedure:: upstream_is_ok
    !procedure:: test
  end type
CONTAINS

!> ============================================================================
!> Read parameters
!> ============================================================================
SUBROUTINE init (self)
  class (rt_patch_t):: self
  integer, save:: n_mu=3, n_phi=4, n_bin=8
  integer:: i_mu, i_phi, i_omega, n_omega
  real, dimension(:), allocatable:: mu, w_mu, phi
  namelist /rt_patch_params/ n_mu, n_phi, n_bin
  !-----------------------------------------------------------------------------
  ! Read parameters
  !-----------------------------------------------------------------------------
  rewind (io%input)
  read (io%input, rt_patch_params)
  if (io%master) write (io%output, rt_patch_params)
  n_omega      = 1 + n_phi*(n_mu-1)
  self%n_omega = n_omega
  self%n_bin   = n_bin
  !-----------------------------------------------------------------------------
  ! Assign angles and integration weight
  !-----------------------------------------------------------------------------
  allocate (mu(n_mu),w_mu(n_mu), phi(n_phi))
  call radau%init (n_mu, mu, w_mu)
  print *, 'MU', mu
  print *, 'WMU', w_mu
  i_omega = 0
  do i_phi=1,n_phi
    i_omega = i_omega+1
    self%rt_sc(i_omega)%omega%mu  = mu(1)
    self%rt_sc(i_omega)%omega%w   = w_mu(1)
    self%rt_sc(i_omega)%omega%phi = 0.0
    do i_mu=2,n_mu
      i_omega = i_omega+1
      self%rt_sc(i_omega)%omega%mu  = mu(i_mu)
      self%rt_sc(i_omega)%omega%w   = w_mu(i_mu)
      self%rt_sc(i_omega)%omega%phi = phi(i_phi)
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Initialize variable indices, using slot 9 for temperature and slot 10 for Q
  !-----------------------------------------------------------------------------
  self%mhd = .true.
  self%nv  = 10
  self%nt  = 5
  call self%idx%init (5, self%mhd)
  self%idx%tt = 9
  self%idx%qr = 10
  call self%idx%copy
  !-----------------------------------------------------------------------------
  ! Allocate the required RT solvers for all angles
  !-----------------------------------------------------------------------------
  allocate (self%rt_sc(self%n_omega))
  allocate (self%rt_time(self%nt))
  !-----------------------------------------------------------------------------
  ! Punt to the next level, initializing patch memory and meshes
  !-----------------------------------------------------------------------------
  call self%mhd_t%init
END SUBROUTINE init

!> ============================================================================
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
!> ============================================================================
SUBROUTINE update (self)
  class (rt_patch_t):: self
  integer:: i_omega
  real, dimension(:,:,:), pointer:: q
  !----------------------------------------------------------------------------
  self%rt_time(self%it) = self%rt_sc(1)%time
  q => self%mem(:,:,:,self%idx%qr,self%it,1)
  do i_omega=2,self%n_omega
    if (self%upstream_is_ok(i_omega)) then
      call self%rt_sc(i_omega)%update
    end if
    q = q + self%rt_sc(i_omega)%q*self%rt_sc(i_omega)%omega%w
    self%rt_time(self%it) = min(self%rt_time(self%it), self%rt_sc(i_omega)%time)
  end do
  !----------------------------------------------------------------------------
  ! If the RT solution is far enough advanced, ask the MHD solver to update.
  ! If this update does not happen, the task time is left unchanged, and the
  ! task will again get a chance to update when the scheduler calls it next
  ! time.
  !----------------------------------------------------------------------------
  if (self%rt_time(self%it) + self%rt_grace_time > self%time) then
    call self%mhd_t%update
  end if
END SUBROUTINE update

!> ============================================================================
!> Check upstream times
!> ============================================================================
FUNCTION upstream_is_ok (self, i_omega)
  logical:: upstream_is_ok
  class (rt_patch_t):: self
  integer:: i_omega
  !-----------------------------------------------------------------------------
  upstream_is_ok = all(self%rt_sc(i_omega)%upstream_times &
                       > self%time - self%grace*self%dtime)
END FUNCTION

END MODULE rt_patch_mod
