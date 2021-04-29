!===============================================================================
!===============================================================================
MODULE viscosity_mod
  USE io_mod
  USE trace_mod
  USE gpatch_mod
  USE scratch_mod
  implicit none
  private
  type, public:: viscosity_t
    class(gpatch_t), pointer:: patch => null()
    real, dimension(:,:,:), pointer, contiguous:: d  => null()
    real, dimension(:,:,:), pointer, contiguous:: u1 => null()
    real, dimension(:,:,:), pointer, contiguous:: u2 => null()
    real, dimension(:,:,:), pointer, contiguous:: u3 => null()
  contains
    procedure:: pre_init
    procedure:: pre_update
    procedure:: set_pointers
  end type
CONTAINS

!===============================================================================
!> Initialize patch link
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(viscosity_t):: self
  class(gpatch_t), pointer:: patch
  !.............................................................................
  call trace%begin ('viscosity_t%pre_init')
  self%patch => patch
  call trace%end()
END SUBROUTINE pre_init

SUBROUTINE set_pointers (self, d, u1, u2, u3)
  class(viscosity_t):: self
  !.............................................................................
  real, dimension(:,:,:), pointer, contiguous:: d
  real, dimension(:,:,:), pointer, contiguous:: u1
  real, dimension(:,:,:), pointer, contiguous:: u2
  real, dimension(:,:,:), pointer, contiguous:: u3
  self%d  => d
  self%u1 => u1
  self%u2 => u2
  self%u3 => u3
END SUBROUTINE set_pointers

!===============================================================================
!> Use the pre-existing force%selected () function to set or add the force
!===============================================================================
SUBROUTINE pre_update (self)
  class(viscosity_t):: self
  !.............................................................................
  real, dimension(:,:,:), pointer, contiguous:: d
  real, dimension(:,:,:), pointer, contiguous:: u1
  real, dimension(:,:,:), pointer, contiguous:: u2
  real, dimension(:,:,:), pointer, contiguous:: u3
  real, dimension(:,:,:), allocatable:: du1   , du2   , du3
  real, dimension(:,:,:), allocatable:: d1u1  , d1u2  , d1u3
  real, dimension(:,:,:), allocatable:: d2u1  , d2u2  , d2u3
  real, dimension(:,:,:), allocatable:: d2u1u2, d2u1u3
  real, dimension(:,:,:), allocatable:: d2u2u3, d2u2u1
  real, dimension(:,:,:), allocatable:: d2u3u1, d2u3u2
  real, dimension(:,:,:,:), pointer, contiguous:: ff
  integer:: i1, i2, i3, m(3)
  real::c1, c2, c3
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('viscosity_t%pre_update', itimer=itimer)
  c1 = 1./self%patch%ds(1)
  c1 = 1./self%patch%ds(2)
  c1 = 1./self%patch%ds(3)
  m = self%patch%gn
  d  => self%d
  u1 => self%u1
  u2 => self%u2
  u3 => self%u3
  allocate (du1   (m(1),m(2),m(3)))
  allocate (du2   (m(1),m(2),m(3)))
  allocate (du3   (m(1),m(2),m(3)))
  allocate (d1u1  (m(1),m(2),m(3)))
  allocate (d1u2  (m(1),m(2),m(3)))
  allocate (d1u3  (m(1),m(2),m(3)))
  allocate (d2u1  (m(1),m(2),m(3)))
  allocate (d2u2  (m(1),m(2),m(3)))
  allocate (d2u3  (m(1),m(2),m(3)))
  allocate (d2u1u2(m(1),m(2),m(3)))
  allocate (d2u1u3(m(1),m(2),m(3)))
  allocate (d2u2u3(m(1),m(2),m(3)))
  allocate (d2u2u1(m(1),m(2),m(3)))
  allocate (d2u3u1(m(1),m(2),m(3)))
  allocate (d2u3u2(m(1),m(2),m(3)))
  do i3=2,m(3)-1
  do i2=2,m(2)-1
  do i1=2,m(1)-1
    du1(i1,i2,i3) = (u1(i1,i2,i3) - u1(i1-1,i2,i3))*c1
    du2(i1,i2,i3) = (u2(i1,i2,i3) - u1(i1,i2-1,i3))*c2
    du3(i1,i2,i3) = (u3(i1,i2,i3) - u1(i1,i2,i3-1))*c3
  end do
  end do
  end do
  do i3=2,m(3)-1
  do i2=2,m(2)-1
  do i1=2,m(1)-1
    d2u1(i1,i2,i3) =  (du1(i1+1,i2,i3) - du1(i1,i2,i3))*d(i1,i2,i3)*c1
    d2u2(i1,i2,i3) =  (du2(i1,i2+1,i3) - du1(i1,i2,i3))*d(i1,i2,i3)*c2
    d2u3(i1,i2,i3) =  (du3(i1,i2,i3+1) - du1(i1,i2,i3))*d(i1,i2,i3)*c3
    d1u1(i1,i2,i3) = ((du1(i1+1,i2,i3) + du1(i1,i2,i3))*d(i1,i2,i3))*0.5*c1
    d1u2(i1,i2,i3) = ((du1(i1,i2+1,i3) + du1(i1,i2,i3))*d(i1,i2,i3))*0.5*c2
    d1u3(i1,i2,i3) = ((du1(i1,i2,i3+1) + du1(i1,i2,i3))*d(i1,i2,i3))*0.5*c3
  end do
  end do
  end do
  do i3=3,m(3)-2
  do i2=3,m(2)-2
  do i1=3,m(1)-2
    d2u1u2(i1,i2,i3) = ((d1u1(i1,i2+1,i3) - du1(i1,i2-1,i3))*d(i1,i2,i3))*0.5*c2
    d2u1u3(i1,i2,i3) = ((d1u1(i1,i2,i3+1) - du1(i1,i2,i3-1))*d(i1,i2,i3))*0.5*c3
    d2u2u3(i1,i2,i3) = ((d1u2(i1,i2,i3+1) - du2(i1,i2,i3-1))*d(i1,i2,i3))*0.5*c3
    d2u2u1(i1,i2,i3) = ((d1u2(i1+1,i2,i3) - du2(i1-1,i2,i3))*d(i1,i2,i3))*0.5*c1
    d2u3u1(i1,i2,i3) = ((d1u3(i1+1,i2,i3) - du1(i1-1,i2,i3))*d(i1,i2,i3))*0.5*c1
    d2u3u2(i1,i2,i3) = ((d1u1(i1,i2+1,i3) - du1(i1,i2-1,i3))*d(i1,i2,i3))*0.5*c2
  end do
  end do
  end do
  ff => scratch%force_per_unit_mass
  do i3=3,m(3)-2
  do i2=3,m(2)-2
  do i1=3,m(1)-2
    ff(i1,i2,i3,1) =    0.5*d(i1  ,i2,i3)*( &
                      (d2u1  (i1+1,i2,i3) - d2u1  (i1-1,i2,i3)) &
                    + (d2u1u2(i1,i2+1,i3) - d2u1u2(i1,i2-1,i3)) &
                    + (d2u1u3(i1,i2,i3+1) - d2u1u3(i1,i2,i3-1)))*c1
    ff(i1,i2,i3,2) =    0.5*d(i1 ,i2,i3)*( &
                      (d2u2u1(i1+1,i2,i3) - d2u1u2(i1-1,i2,i3)) &
                    + (d2u2  (i1,i2+1,i3) - d2u1  (i1,i2-1,i3)) &
                    + (d2u2u3(i1,i2,i3+1) - d2u1u2(i1,i2,i3-1)))*c2
    ff(i1,i2,i3,3) =    0.5*d(i1  ,i2,i3)*( &
                      (d2u3u1(i1+1,i2,i3) - d2u1u2(i1-1,i2,i3)) &
                    + (d2u3u2(i1,i2+1,i3) - d2u1u2(i1,i2-1,i3)) &
                    + (d2u3  (i1,i2,i3+1) - d2u1  (i1,i2,i3-1)))*c3
  end do
  end do
  end do
  call trace%end (itimer)
END SUBROUTINE pre_update

END MODULE viscosity_mod
