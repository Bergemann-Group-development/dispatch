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
!> (not overwrite) the generic patch_t%force_per_unit_mass or _volume arrays.
!>
!> To comply with that policy, this interface also adds the force, and thus
!> assumes that the generic array is set to zero elsewhere, before each timestep.
!> The proper place to do that is in the extras_t data type, which is where the
!> user chooses which forces to include and combine.
!===============================================================================
MODULE forces_mod
  USE io_mod
  USE mpi_mod
  USE trace_mod
  USE patch_mod
  USE force_mod
  USE link_mod
  USE mpi_global_mod
  implicit none
  private
  type, public:: forces_t
    type(force_t):: force
    class(patch_t), pointer:: patch
    real, dimension(:,:,:), pointer:: Ux=>null(), Uy=>null(), Uz=>null()
    real(8):: p_glob(3)=0.0_8
    real(8):: l_aver(3)=0.0_8
  contains
    procedure:: init
    procedure:: dealloc
    procedure:: pre_update
    procedure:: post_update
  end type
  type(mpi_global_t):: p_global
  real(8):: p_aver(3)=0d0
  real(8):: p_prev(3)=0d0
CONTAINS

!===============================================================================
! Initialize forcing
!===============================================================================
SUBROUTINE init (self, link)
  class(forces_t):: self
  class(link_t):: link
  !.............................................................................
  call trace%begin ('forces_t%init')
  self%patch => task2patch (link%task)
  call self%force%init (self%patch%kind, self%patch%id, self%patch%mesh)
  !$omp critical (update_global_cr)
  call p_global%init (3)
  !$omp end critical (update_global_cr)
  call update_global (self)
  call trace%end
END SUBROUTINE init

!===============================================================================
!> Update the global momentum average.  p_aver (computed below) is the current
!> contribution to the global linear momentum.  self%p_prev is the previous
!> contribution, and self%p_glob is the resulting global linear momentum.
!===============================================================================
SUBROUTINE update_global (self)
  class(forces_t):: self
  class(patch_t), pointer:: patch
  integer:: i, iv, l(3), u(3)
  real(8):: w, l_aver(3)
  integer, save:: updates=0
  !-----------------------------------------------------------------------------
  if (self%force%t_damp <= 0.0) return
  call trace%begin ('forces_t%update_global')
  !-----------------------------------------------------------------------------
  patch => self%patch
  l = patch%mesh%li
  u = patch%mesh%ui
  w = product(patch%mesh%d/patch%mesh%b)
  !-----------------------------------------------------------------------------
  ! For every task update, update the patch-local contribution
  !-----------------------------------------------------------------------------
  do i=1,3
    iv = patch%idx%px+i-1
    l_aver(i) = w*sum(patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,patch%it,1))
    !$omp atomic update
    p_aver(i) = p_aver(i) + l_aver(i) - self%l_aver(i)
    self%l_aver(i) = l_aver(i)
  end do
  !-----------------------------------------------------------------------------
  ! Every io%nwrite updates, update the global average
  !-----------------------------------------------------------------------------
  !$omp atomic update
  updates = updates+1
  if (updates >= io%nwrite) then
    !$omp critical (update_global_cr)
    if (updates >= io%nwrite) then
      updates = 0    
      self%p_glob = p_global%update (p_aver-p_prev)
      if (self%force%verbose > 0) then
        write(*,'(a,i3, 2i6,f12.6,3(2x,1p,3e10.2))') 'update_global:', &
          mpi%rank, patch%id, patch%istep, patch%time, self%p_glob, p_aver-p_prev
      end if
      p_prev = p_aver
    end if
    !$omp end critical (update_global_cr)
  end if
  call trace%end
END SUBROUTINE update_global

!===============================================================================
!> Deallocate
!===============================================================================
SUBROUTINE dealloc (self)
  class(forces_t):: self
  !.............................................................................
  call trace%begin ('forces_t%dealloc')
  call self%force%dealloc
  call trace%end
END SUBROUTINE dealloc

!===============================================================================
!> Use the pre-existing force%selected () function to set or add the force
!===============================================================================
SUBROUTINE pre_update (self)
  class(forces_t):: self
  !.............................................................................
  real, dimension(:,:,:), pointer:: d
  real, dimension(:,:,:,:), pointer:: p
  integer:: i
  real(8):: p_glob
  real, save:: dtime=1e-3
  !-----------------------------------------------------------------------------
  ! Use existing function call
  !-----------------------------------------------------------------------------
  associate (patch => self%patch, m => self%patch%gn)
  d => patch%mem(:,:,:,patch%idx%d,patch%it,1)
  p => patch%mem(:,:,:,patch%idx%px:patch%idx%pz,patch%it,1)
  if (.not.allocated(patch%force_per_unit_mass)) then
    allocate (patch%force_per_unit_mass(m(1),m(2),m(3),3))
    call io%bits_mem (storage_size(patch%force_per_unit_mass), &
                     product(shape(patch%force_per_unit_mass)), 'fpm')
  end if
  call update_global (self)
  patch%force_per_unit_mass = &
  self%force%selected(patch%time, d, p, self%Ux, self%Uy, self%Uz, patch%mesh)
  if (patch%istep > 0 .and. self%force%t_damp > 0.0) then
    do i=1,3
      !$omp atomic read
      p_glob = self%p_glob(i)
      patch%force_per_unit_mass(:,:,:,i) = &
      patch%force_per_unit_mass(:,:,:,i) - p_glob/self%force%t_damp
    end do
  end if
  end associate
END SUBROUTINE pre_update

!===============================================================================
!> Save on node memory, by deallocating
!===============================================================================
SUBROUTINE post_update (self)
  class(forces_t):: self
  !-----------------------------------------------------------------------------
  if (allocated(self%patch%force_per_unit_mass)) then
    call io%bits_mem (-storage_size(self%patch%force_per_unit_mass), &
                      product(shape(self%patch%force_per_unit_mass)), '-fpm')
    deallocate (self%patch%force_per_unit_mass)
  end if
END SUBROUTINE post_update

END MODULE forces_mod
