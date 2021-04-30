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
    real:: mu=1.0
    real:: u_max
  contains
    procedure:: courant_condition
    procedure:: pre_init
    procedure:: pre_update
    procedure:: set_pointers
    procedure:: unset_pointers
  end type
  logical, save:: on=.false.
CONTAINS

!===============================================================================
!> Initialize patch link
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(viscosity_t):: self
  class(gpatch_t), target:: patch
  logical, save:: first_time=.true.
  real:: mu=1.0
  integer                           :: iostat
  namelist /viscosity_params/ on, mu
  !.............................................................................
  call trace%begin ('viscosity_t%pre_init')
  self%patch => patch
  call io%lock%set ('viscosity_t%pre_init')
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, viscosity_params, iostat=iostat)
    call io%namelist_warning ('viscosity_params', iostat)
    write (io%namelist, viscosity_params)
    self%mu = mu
  end if
  call io%lock%unset ('viscosity_t%pre_init')
  call trace%end()
END SUBROUTINE pre_init

SUBROUTINE set_pointers (self)
  class(viscosity_t):: self
  !.............................................................................
  integer:: m(3)
  ! ----------------------------------------------------------------------------
  m = self%patch%gn

  self%d  => self%patch%mem(:,:,:,self%patch%idx%d,self%patch%it,1)
  allocate(self%u1(m(1),m(2),m(3)))
  allocate(self%u2(m(1),m(2),m(3)))
  allocate(self%u3(m(1),m(2),m(3)))
  self%u1 = self%patch%mem(:,:,:,self%patch%idx%px,self%patch%it,1) / &
            self%patch%mem(:,:,:,self%patch%idx%d ,self%patch%it,1)
  self%u2 = self%patch%mem(:,:,:,self%patch%idx%py,self%patch%it,1) / &
            self%patch%mem(:,:,:,self%patch%idx%d ,self%patch%it,1)
  self%u3 = self%patch%mem(:,:,:,self%patch%idx%pz,self%patch%it,1) / &
            self%patch%mem(:,:,:,self%patch%idx%d ,self%patch%it,1)
END SUBROUTINE set_pointers

SUBROUTINE unset_pointers (self)
  class(viscosity_t):: self
  !.............................................................................
  real, dimension(:,:,:), pointer, contiguous:: d
  real, dimension(:,:,:), pointer, contiguous:: u1
  real, dimension(:,:,:), pointer, contiguous:: u2
  real, dimension(:,:,:), pointer, contiguous:: u3

  nullify (self%d)
  deallocate(self%u1)
  deallocate(self%u2)
  deallocate(self%u3)
END SUBROUTINE unset_pointers

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
  real(8):: ds(3), u_max, c1, c2, c3
  integer:: i1, i2, i3, m(3), nn(4), li(3), ui(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.on) return
  call trace%begin ('viscosity_t%pre_update', itimer=itimer)
  call self%set_pointers()
   c1 = 1./self%patch%ds(1)
   c2 = 1./self%patch%ds(2)
   c3 = 1./self%patch%ds(3)
  ds = self%patch%ds
  m = self%patch%gn
  li=self%patch%li
  ui=self%patch%ui
  nn(1:3)=m
  nn(4) = 3
  call scratch%allocate (scratch%force_per_unit_mass, nn)
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
  do i3=2,m(3)
  do i2=2,m(2)
  do i1=2,m(1)
    du1(i1,i2,i3) = u1(i1,i2,i3) - u1(i1-1,i2,i3)*c1 ! 1/t
    du2(i1,i2,i3) = u2(i1,i2,i3) - u2(i1,i2-1,i3)*c2 ! 1/t
    du3(i1,i2,i3) = u3(i1,i2,i3) - u3(i1,i2,i3-1)*c3 ! 1/t
  end do
  end do
  end do
  do i3=2,m(3)-1
  do i2=2,m(2)-1
  do i1=2,m(1)-1
    d2u1(i1,i2,i3) = (du1(i1+1,i2,i3) - du1(i1,i2,i3))*d(i1,i2,i3)*c1 ! m/l^4t
    d2u2(i1,i2,i3) = (du2(i1,i2+1,i3) - du2(i1,i2,i3))*d(i1,i2,i3)*c2 ! m/l^4t
    d2u3(i1,i2,i3) = (du3(i1,i2,i3+1) - du3(i1,i2,i3))*d(i1,i2,i3)*c3 ! m/l^4t
    d1u1(i1,i2,i3) = (du1(i1+1,i2,i3) + du1(i1,i2,i3))*0.5 ! 1/t *d(i1,i2,i3))*0.5*c1 ! m/l^4t
    d1u2(i1,i2,i3) = (du2(i1,i2+1,i3) + du2(i1,i2,i3))*0.5 ! 1/t*d(i1,i2,i3))*0.5*c2 ! m/l^4t
    d1u3(i1,i2,i3) = (du3(i1,i2,i3+1) + du3(i1,i2,i3))*0.5 ! 1/t*d(i1,i2,i3))*0.5*c3 ! m/l^4t
  end do
  end do
  end do
  do i3=3,m(3)-2
  do i2=3,m(2)-2
  do i1=3,m(1)-2
    d2u1u2(i1,i2,i3) = ((d1u1(i1,i2+1,i3) - du2(i1,i2-1,i3))*d(i1,i2,i3))*0.5*c2 ! m/l^4t
    d2u1u3(i1,i2,i3) = ((d1u1(i1,i2,i3+1) - du3(i1,i2,i3-1))*d(i1,i2,i3))*0.5*c3 ! m/l^4t
    d2u2u3(i1,i2,i3) = ((d1u2(i1,i2,i3+1) - du3(i1,i2,i3-1))*d(i1,i2,i3))*0.5*c3 ! m/l^4t
    d2u2u1(i1,i2,i3) = ((d1u2(i1+1,i2,i3) - du2(i1-1,i2,i3))*d(i1,i2,i3))*0.5*c1 ! m/l^4t
    d2u3u1(i1,i2,i3) = ((d1u3(i1+1,i2,i3) - du1(i1-1,i2,i3))*d(i1,i2,i3))*0.5*c1 ! m/l^4t
    d2u3u2(i1,i2,i3) = ((d1u1(i1,i2+1,i3) - du1(i1,i2-1,i3))*d(i1,i2,i3))*0.5*c2 ! m/l^4t
  end do
  end do
  end do
  ff => scratch%force_per_unit_mass
  do i3=li(3), ui(3)
  do i2=li(2), ui(2)
  do i1=li(1), ui(1)
    ff(i1,i2,i3,1) =     0.5*self%mu/d(i1,i2,i3)*self%patch%ds(1)*( &
                      (d2u1  (i1+1,i2,i3) - d2u1  (i1-1,i2,i3)) &
                    + (d2u1u2(i1,i2+1,i3) - d2u1u2(i1,i2-1,i3)) &
                    + (d2u1u3(i1,i2,i3+1) - d2u1u3(i1,i2,i3-1)))
    ff(i1,i2,i3,2) =     0.5*self%mu/d(i1,i2,i3)*self%patch%ds(2)*( &
                      (d2u2u1(i1+1,i2,i3) - d2u2u1(i1-1,i2,i3)) &
                    + (d2u2  (i1,i2+1,i3) - d2u2  (i1,i2-1,i3)) &
                    + (d2u2u3(i1,i2,i3+1) - d2u2u3(i1,i2,i3-1)))
    ff(i1,i2,i3,3) =     0.5*self%mu/d(i1,i2,i3)*self%patch%ds(3)*( &
                      (d2u3u1(i1+1,i2,i3) - d2u3u1(i1-1,i2,i3)) &
                    + (d2u3u2(i1,i2+1,i3) - d2u3u2(i1,i2-1,i3)) &
                    + (d2u3  (i1,i2,i3+1) - d2u3  (i1,i2,i3-1)))
  end do
  end do
  end do
  u_max = maxval(ff(li(1):ui(1),li(2):ui(2),li(3):ui(3),1))!*ds(1)
 ! u_max = u_max/(1.0+(self%mu*ds(1)**2))
  !print *, 'umax x', u_max
  self%u_max = max(self%u_max,u_max)
  u_max = maxval(ff(li(1):ui(1),li(2):ui(2),li(3):ui(3),2))!*ds(2)
  !u_max = u_max/(1.0+(self%mu*ds(2)**2))
 ! print *, 'umax y', u_max
  self%u_max = max(self%u_max,u_max)
  u_max = maxval(ff(li(1):ui(1),li(2):ui(2),li(3):ui(3),3))!*ds(3)
  !u_max = u_max/(1.0+(self%mu*ds(3)**2))
  !print *, 'umax z', u_max
  self%u_max = max(self%u_max,u_max)
  call self%unset_pointers()
  call trace%end (itimer)
END SUBROUTINE pre_update

SUBROUTINE courant_condition(self)
  class(viscosity_t):: self
  !print *, 'visc u_max', self%u_max, self%patch%u_max
  self%patch%u_max = max(self%patch%u_max,self%u_max)

END SUBROUTINE courant_condition

END MODULE viscosity_mod
