!===============================================================================
!> This is a bare minimum interface to the new extras module, showing how to
!> interface with a pre-existing force_mod%selected() function.
!>
!> In this version, the force_per_unit_mass array is allocated in %init, and its
!> values are set by a call to the previous standard force_t interface, in order
!> for the new arrangement to be as non-intrusive as possible, relative to pre-
!> existing experiments.
!>
!> In the future, external forces should NOT use this interface.  Instead, each
!> type of force should use its own data type (primarily to hold parameters
!> specific to that force), while the resulting force should be addded to
!> (not overwrite) the generic scratch_t%force_per_unit_mass or _volume arrays.
!>
!> To comply with that policy, this interface also adds the force, and thus
!> assumes that the generic array is set to zero elsewhere, before each timestep.
!> The proper place to do that is in the extras_t data type, which is where the
!> user chooses which forces to include and combine.
!===============================================================================
MODULE forces_mod
  USE io_mod
  USE trace_mod
  USE gpatch_mod
  USE mesh_mod
  USE link_mod
  USE index_mod
  USE scratch_mod
  USE model_1d_mod
  USE interpolate_mod
  implicit none
  private
  type, public:: forces_t
    class(gpatch_t), pointer:: patch => null()
    real, dimension(:,:,:), pointer:: Ux=>null(), Uy=>null(), Uz=>null()
  contains
    procedure:: pre_init
    procedure:: pre_update
  end type
  real:: margin=0.0
CONTAINS

!===============================================================================
!> Initialize patch link
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(forces_t):: self
  class(gpatch_t), pointer:: patch
  !.............................................................................
  call trace%begin ('forces_t%pre_init')
  self%patch => patch
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> Use the pre-existing force%selected () function to set or add the force
!===============================================================================
SUBROUTINE pre_update (self)
  class(forces_t)  :: self
  !.............................................................................
  class(mesh_t), pointer :: m1, m2, m3
  real, pointer          :: gz(:,:,:)
  integer                :: i, i1, i2, i3
  real                   :: x1, x2, x3, rhat(3)
  real(8)                :: r, g
  integer, save          :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('forces_t%pre_update', itimer=itimer)
  call scratch%allocate (scratch%force_per_unit_mass, [self%patch%gn,3])
  i = 0
  m1 => self%patch%mesh(1)
  m2 => self%patch%mesh(2)
  m3 => self%patch%mesh(3)
  do i3=m3%lb,m3%ub
    do i2=m2%lb,m2%ub
      x2 = m2%r(i2)
      do i1=m1%lb,m1%ub
        x1 = m1%r(i1)
        r = sqrt(m1%p**2+m2%p**2+m3%p**2)
        x3 = m3%r(i3) + r
        r = sqrt(x1**2+x2**2+x3**2)
        rhat = [x1,x2,x3]/r
        g = interpolate%one_d (r, model_1d%r, model_1d%g, i)
        scratch%force_per_unit_mass(i1,i2,i3,:) = -g*rhat
      end do
    end do
  end do
  gz => scratch%force_per_unit_mass(:,:,:,3)
  call self%patch%aux%register('gz',gz)
  call trace%end (itimer)
END SUBROUTINE pre_update

END MODULE forces_mod
