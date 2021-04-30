!===============================================================================
!===============================================================================
MODULE viscosity_mod
  USE io_mod
  USE trace_mod
  USE gpatch_mod
  USE scratch_mod
  USE scaling_mod
  USE kinds_mod
  implicit none
  private
  type, public:: viscosity_t
    class(gpatch_t), pointer:: patch => null()
    real, dimension(:,:,:), pointer, contiguous:: d  => null()
    real, dimension(:,:,:), pointer, contiguous:: u1 => null()
    real, dimension(:,:,:), pointer, contiguous:: u2 => null()
    real, dimension(:,:,:), pointer, contiguous:: u3 => null()
    real:: mu = 1.0
    real:: u_max = 0.0
    logical:: on = .true.
    logical:: dynamic = .true.
    logical:: allocated=.false.
  contains
    procedure:: courant_condition
    procedure:: pre_init
    procedure:: pre_update
    procedure:: set_pointers
    procedure, private:: deallocate
  end type
CONTAINS

!===============================================================================
!> Initialize patch link
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(viscosity_t):: self
  class(gpatch_t), target:: patch
  integer:: iostat
  real    :: mu = 1.0
  logical :: on = .false.
  logical :: dynamic = .true.
  logical :: first_time = .true.
  namelist /viscosity_params/ on, dynamic, mu
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
  end if
  self%on = on
  self%dynamic = dynamic
  self%mu = mu * scaling%t / scaling%l**2
  call io%lock%unset ('viscosity_t%pre_init')
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!===============================================================================
SUBROUTINE set_pointers (self, u1, u2, u3)
  class(viscosity_t):: self
  real(var_k), dimension(:,:,:), pointer, contiguous, optional:: u1, u2, u3
  !.............................................................................
  integer:: m(3)
  !-----------------------------------------------------------------------------
  m = self%patch%gn
  self%d  => self%patch%mem(:,:,:,self%patch%idx%d,self%patch%it,1)
  if (present(u1)) then
    self%allocated = .false.
    self%u1 => u1
    self%u2 => u2
    self%u3 => u3
  else
    self%allocated = .true.
    allocate(self%u1(m(1),m(2),m(3)))
    allocate(self%u2(m(1),m(2),m(3)))
    allocate(self%u3(m(1),m(2),m(3)))
    self%u1 = self%patch%mem(:,:,:,self%patch%idx%px,self%patch%it,1) / &
              self%patch%mem(:,:,:,self%patch%idx%d ,self%patch%it,1)
    self%u2 = self%patch%mem(:,:,:,self%patch%idx%py,self%patch%it,1) / &
              self%patch%mem(:,:,:,self%patch%idx%d ,self%patch%it,1)
    self%u3 = self%patch%mem(:,:,:,self%patch%idx%pz,self%patch%it,1) / &
              self%patch%mem(:,:,:,self%patch%idx%d ,self%patch%it,1)
  end if
END SUBROUTINE set_pointers

!===============================================================================
SUBROUTINE deallocate (self)
  class(viscosity_t):: self
  !.............................................................................
  real, dimension(:,:,:), pointer, contiguous:: d
  real, dimension(:,:,:), pointer, contiguous:: u1
  real, dimension(:,:,:), pointer, contiguous:: u2
  real, dimension(:,:,:), pointer, contiguous:: u3
  !-----------------------------------------------------------------------------
  nullify (self%d)
  deallocate(self%u1)
  deallocate(self%u2)
  deallocate(self%u3)
END SUBROUTINE deallocate

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
  real, dimension(:,:,:), allocatable:: du11, du12, du13
  real, dimension(:,:,:), allocatable:: du21, du22, du23
  real, dimension(:,:,:), allocatable:: du31, du32, du33
  real, dimension(:,:,:), allocatable:: divu
  real, dimension(:,:,:,:), pointer, contiguous:: ff
  real(8):: ds(3)
  real:: c1, c2, c3, c, u_max(3)
  integer:: i1, i2, i3, m(3), nn(4), li(3), ui(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  call trace%begin ('viscosity_t%pre_update', itimer=itimer)
  call self%set_pointers()
  ds = self%patch%ds
  m = self%patch%gn
  nn(1:3)=m
  nn(4) = 3
  call scratch%allocate (scratch%force_per_unit_mass, nn)
  d  => self%d
  u1 => self%u1
  u2 => self%u2
  u3 => self%u3
  allocate (du11(m(1),m(2),m(3)), du12(m(1),m(2),m(3)), du13(m(1),m(2),m(3)))
  allocate (du21(m(1),m(2),m(3)), du22(m(1),m(2),m(3)), du23(m(1),m(2),m(3)))
  allocate (du31(m(1),m(2),m(3)), du32(m(1),m(2),m(3)), du33(m(1),m(2),m(3)))
  allocate (divu(m(1),m(2),m(3)))
  !-----------------------------------------------------------------------------
  ! Centered partial derivatives
  !-----------------------------------------------------------------------------
  c1 = 0.5/ds(1)
  c2 = 0.5/ds(2)
  c3 = 0.5/ds(3)
  li=self%patch%li
  ui=self%patch%ui
  do i3=li(3)-1,ui(3)+1
  do i2=li(2)-1,ui(2)+1
  do i1=li(1)-1,ui(1)+1
    du11(i1,i2,i3) = (u1(i1+1,i2,i3) - u1(i1-1,i2,i3))*c1
    du21(i1,i2,i3) = (u2(i1+1,i2,i3) - u2(i1-1,i2,i3))*c1
    du31(i1,i2,i3) = (u2(i1+1,i2,i3) - u3(i1-1,i2,i3))*c1
    du12(i1,i2,i3) = (u1(i1,i2+1,i3) - u1(i1,i2-1,i3))*c2
    du22(i1,i2,i3) = (u2(i1,i2+1,i3) - u2(i1,i2-1,i3))*c2
    du32(i1,i2,i3) = (u3(i1,i2+1,i3) - u3(i1,i2-1,i3))*c2
    du13(i1,i2,i3) = (u1(i1,i2,i3+1) - u1(i1,i2,i3-1))*c3
    du23(i1,i2,i3) = (u2(i1,i2,i3+1) - u2(i1,i2,i3-1))*c3
    du33(i1,i2,i3) = (u3(i1,i2,i3+1) - u3(i1,i2,i3-1))*c3
  end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Compute divu and symmetrize
  !-----------------------------------------------------------------------------
  do i3=li(3)-1,ui(3)+1
  do i2=li(2)-1,ui(2)+1
  do i1=li(1)-1,ui(1)+1
    divu(i1,i2,i3) = du11(i1,i2,i3) &
                   + du22(i1,i3,i3) &
                   + du33(i1,i2,i3)
    du11(i1,i2,i3) = du11(i1,i2,i3) + du11(i1,i2,i3)
    du22(i1,i2,i3) = du22(i1,i2,i3) + du22(i1,i2,i3)
    du33(i1,i2,i3) = du33(i1,i2,i3) + du33(i1,i2,i3)
    du12(i1,i2,i3) = du21(i1,i2,i3) + du12(i1,i2,i3)
    du13(i1,i2,i3) = du31(i1,i2,i3) + du13(i1,i2,i3)
    du23(i1,i2,i3) = du32(i1,i2,i3) + du23(i1,i2,i3)
  end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Subtract compressive part and make per-volume
  !-----------------------------------------------------------------------------
  c = 2./3.
  do i3=li(3)-1,ui(3)+1
  do i2=li(2)-1,ui(2)+1
    if (self%dynamic) then
      do i1=li(1)-1,ui(1)+1
        du11(i1,i2,i3) = (du11(i1,i2,i3) - c*divu(i1,i2,i3))
        du22(i1,i2,i3) = (du22(i1,i2,i3) - c*divu(i1,i2,i3))
        du33(i1,i2,i3) = (du33(i1,i2,i3) - c*divu(i1,i2,i3))
        du12(i1,i2,i3) =  du12(i1,i2,i3)
        du13(i1,i2,i3) =  du13(i1,i2,i3)
        du23(i1,i2,i3) =  du23(i1,i2,i3)
      end do
    else
      do i1=li(1)-1,ui(1)+1
        du11(i1,i2,i3) = d(i1,i2,i3)*(du11(i1,i2,i3) - c*divu(i1,i2,i3))
        du22(i1,i2,i3) = d(i1,i2,i3)*(du22(i1,i2,i3) - c*divu(i1,i2,i3))
        du33(i1,i2,i3) = d(i1,i2,i3)*(du33(i1,i2,i3) - c*divu(i1,i2,i3))
        du12(i1,i2,i3) = d(i1,i2,i3)* du12(i1,i2,i3)
        du13(i1,i2,i3) = d(i1,i2,i3)* du13(i1,i2,i3)
        du23(i1,i2,i3) = d(i1,i2,i3)* du23(i1,i2,i3)
      end do
    end if
  end do
  end do
!-----------------------------------------------------------------------------
  ! Divergence of per-volume quantity, in per unit mass
  !-----------------------------------------------------------------------------
  ff => scratch%force_per_unit_mass
  li=self%patch%li
  ui=self%patch%ui
  do i3=li(3),ui(3)
  do i2=li(2),ui(2)
  do i1=li(1),ui(1)
    c = 0.5*self%mu/d(i1,i2,i3)
    ff(i1,i2,i3,1) = ((du11(i1+1,i2,i3) - du11(i1-1,i2,i3))*c1 &
                   +  (du12(i1,i2+1,i3) - du12(i1,i2-1,i3))*c2 &
                   +  (du13(i1,i2,i3+1) - du13(i1,i2,i3-1))*c3)*c
    ff(i1,i2,i3,2) = ((du12(i1+1,i2,i3) - du12(i1-1,i2,i3))*c1 &
                   +  (du22(i1,i2+1,i3) - du22(i1,i2-1,i3))*c2 &
                   +  (du23(i1,i2,i3+1) - du23(i1,i2,i3-1))*c3)*c
    ff(i1,i2,i3,3) = ((du13(i1+1,i2,i3) - du13(i1-1,i2,i3))*c1 &
                   +  (du23(i1,i2+1,i3) - du23(i1,i2-1,i3))*c2 &
                   +  (du33(i1,i2,i3+1) - du33(i1,i2,i3-1))*c3)*c
  end do
  end do
  end do
  self%u_max = self%mu/minval(ds)**2
  !-----------------------------------------------------------------------------
  deallocate (du11,du12,du13)
  deallocate (du21,du22,du23)
  deallocate (du31,du32,du33)
  if (self%allocated) then
    call self%deallocate ()
  end if
  call trace%end (itimer)
END SUBROUTINE pre_update

!===============================================================================
!===============================================================================
SUBROUTINE courant_condition(self)
  class(viscosity_t):: self
  !-----------------------------------------------------------------------------
  !print *, 'visc u_max', self%u_max, self%patch%u_max
  self%patch%u_max = max (self%patch%u_max, self%u_max)
END SUBROUTINE courant_condition

END MODULE viscosity_mod
