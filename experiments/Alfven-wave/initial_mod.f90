!> =============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!> =============================================================================
MODULE initial_mod
  !USE const
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE random_mod
  USE omp_mod
  USE index_mod
  USE io_unit_mod
  USE force_mod
  implicit none
  private
  type, public:: initial_t
    integer:: seed
    real:: gamma=5./3.
    real:: d0
    logical:: mhd
    real(8):: time=0d0
    character(len=64):: solver
  contains
    procedure:: init
    procedure:: condition
  end type
  real, save:: csound=1.
CONTAINS

!> =============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!> =============================================================================
SUBROUTINE init (self, solver, gamma)
  class(initial_t):: self
  character(len=64):: solver
  real, optional:: gamma
  !---------------------------------------------------------------------
  call trace_begin('initial_t%init')
  self%solver = trim(solver)
  
  call trace%end
END subroutine init


SUBROUTINE condition (self, m, ff, idx)
  class(initial_t) :: self
  type(index_t):: idx
  class(mesh_t)    :: m(3)
  real, dimension(:,:,:,:), pointer :: ff
  ! ....................................................................
  integer:: i, j, k
  real:: x, y, z
  logical:: first_time = .true.
  real, dimension(:,:,:,:), pointer:: b, p
  real, dimension(:,:,:), pointer:: d, vx, vy, vz, Bx, By, Bz, px, py, pz, s, e1, et
  real,save:: lambda,k_w
  
  !---------------------------------------------------------------------
  ! Allocate memory and set default values
  !---------------------------------------------------------------------
  if (trim(self%solver(1:6)) == 'ramses') then 
    et => ff(:,:,:,idx%e)
    allocate (e1(m(1)%gn,m(2)%gn,m(3)%gn))
  end if

  d => ff(:,:,:,idx%d)
  if (trim(self%solver(1:7)) == 'stagger') s => ff(:,:,:,idx%s)
  p => ff(:,:,:,idx%px:idx%pz)
  B => ff(:,:,:,idx%bx:idx%bz)
  px=>p(:,:,:,1); py=>p(:,:,:,2); pz=>p(:,:,:,3);
  Bx=>B(:,:,:,1); By=>B(:,:,:,2); Bz=>B(:,:,:,3);
  !
  lambda=1.
  k_w=2.*acos(-1.0d0)/lambda
  !
  do i=m(1)%lb,m(1)%ub
    do j=m(2)%lb,m(2)%ub
      do k=m(3)%lb,m(3)%ub
        x = m(1)%p + m(1)%r(i)
        y = m(2)%p + m(2)%r(j)
        z = m(3)%p + m(3)%r(k)
        d (i,j,k)= 1.
        px(i,j,k)= d(i,j,k) * 1.*cos(k_w*z)
        py(i,j,k)=-d(i,j,k) * 1.*sin(k_w*z)
        pz(i,j,k)= 0.
        Bx(i,j,k)= 1.*cos(k_w*z) 
        By(i,j,k)=-1.*sin(k_w*z)
        Bz(i,j,k)= 1.
        if (trim(self%solver(1:6)) == 'ramses') then
          et(i,j,k)= 1./(self%gamma-1d0) + 0.5*d(i,j,k)*(cos(k_w*z)**2 + sin(k_w*z)**2) + 0.5 * (Bx(i,j,k)**2 + By(i,j,k)**2 + Bz(i,j,k)**2)  
        end if
        if (trim(self%solver(1:7)) == 'stagger') then
          s(i,j,k)= (log(1.)-log(d(i,j,k))*self%gamma)/(self%gamma-1.0) * d(i,j,k)
        end if
      end do
    end do
  end do
  call trace%end
END SUBROUTINE condition

END MODULE initial_mod
