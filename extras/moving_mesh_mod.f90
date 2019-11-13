!===============================================================================
!> Add the extra functionality needed to handle the motion of patches
!===============================================================================
MODULE moving_mesh_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE link_mod
  USE patch_mod
  implicit none
  private
  type, public:: moving_mesh_t
    class(patch_t), pointer:: patch
    real(8):: velocity(3)=0.0
  contains
    procedure:: pre_init
    procedure:: post_update
  end type
  logical, save:: on=.false.
  logical, save:: first_time=.true.
CONTAINS

!===============================================================================
!> Initialize
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class (moving_mesh_t), target:: self
  class (patch_t), pointer:: patch
  integer:: iostat
  namelist /moving_mesh_params/ on
  !-----------------------------------------------------------------------------
  call trace%begin ('moving_mesh_t%init')
  self%patch => patch
  !$omp critical (input_cr)
  if (first_time) then
    rewind (io%input)
    read (io%input, moving_mesh_params, iostat=iostat)
    write (io%output, moving_mesh_params)
    first_time = .false.
  end if
  !$omp end critical (input_cr)
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> After update
!===============================================================================
SUBROUTINE post_update (self)
  class(moving_mesh_t):: self
  real, dimension(:,:,:), pointer, contiguous:: d, p1, p2, p3
  real(8):: mass, momentum(3)
  integer:: i1, i2, i3, l(3), u(3)
  !----------------------------------------------------------------------------
  if (.not. on) return
  call trace%begin ('moving_mesh_t%post_update')
  !----------------------------------------------------------------------------
  ! Compute average momentum and velocity
  !----------------------------------------------------------------------------
  d  => self%patch%mem(:,:,:,self%patch%idx%d ,self%patch%it,1)
  p1 => self%patch%mem(:,:,:,self%patch%idx%px,self%patch%it,1)
  p2 => self%patch%mem(:,:,:,self%patch%idx%py,self%patch%it,1)
  p3 => self%patch%mem(:,:,:,self%patch%idx%pz,self%patch%it,1)
  mass = 0d0
  momentum = 0d0
  l = self%patch%li
  u = self%patch%ui
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    mass     = mass + d(i1,i2,i3)
    momentum = momentum + [p1(i1,i2,i3), &
                           p2(i1,i2,i3), &
                           p3(i1,i2,i3)]
  end do
  end do
  end do
  mass     = mass    /product(u-l+1)
  momentum = momentum/product(u-l+1)
  self%velocity = momentum/mass
  !----------------------------------------------------------------------------
  ! Make a Galilean transformation that cancels the average velocity
  !----------------------------------------------------------------------------
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    p1(i1,i2,i3) = p1(i1,i2,i3) - d(i1,i2,i3)*self%velocity(1)
    p2(i1,i2,i3) = p1(i1,i2,i3) - d(i1,i2,i3)*self%velocity(2)
    p3(i1,i2,i3) = p1(i1,i2,i3) - d(i1,i2,i3)*self%velocity(3)
  end do
  end do
  end do
  call trace%end ()
END SUBROUTINE post_update

END MODULE moving_mesh_mod
