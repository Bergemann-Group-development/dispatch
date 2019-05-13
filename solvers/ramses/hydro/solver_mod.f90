!===============================================================================
!> The solver_t data type extends the patch_t with a generic guard zone dnload
!> procedure, which covers the need for all simple, Cartesian mesh based 
!> solvers, and eliminates the need for a dnload procedure in each one of them.
!> It also adds restart functionality when initializing patches.
!===============================================================================
MODULE solver_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE task_mod
  USE kinds_mod
  USE hd_mod
  USE stagger_mod
  implicit none
  private
  type, public, extends(hd_t):: solver_t
  contains
    procedure:: init
    procedure:: dealloc
    procedure, nopass:: cast2solver
    procedure:: update
    procedure:: dnload
    procedure:: p2U
    procedure:: U2p
    procedure:: e2E_th
    procedure:: E_th2e
    procedure:: log_density => void
    procedure:: log_pressure => void
    procedure:: gas_pressure_ngz
    procedure:: gas_temperature
    procedure:: gas_temperature_ngz
    procedure:: velocity_magnitude => void
    procedure:: grav_potential => void
    procedure:: apply_heating
    procedure:: gas_velocity_vector
    procedure:: gas_velocity_scalar
    procedure:: compression_magnitude
    procedure:: vorticity_magnitude
    generic:: gas_velocity => gas_velocity_vector, gas_velocity_scalar
  end type
  type(solver_t), public:: solver
CONTAINS

!===============================================================================
!> Organize calls to the extras and the hydro solver
!===============================================================================
SUBROUTINE init (self)
  class(solver_t) :: self
  !.............................................................................
  call trace%begin ('solver_t%init')
  call self%hd_t%init
  call self%extras_t%init
  call trace%end
END SUBROUTINE init

!===============================================================================
!> Call a hiearchy of %dealloc procedures, which should deallocate everything
!> allocated in the task
!===============================================================================
SUBROUTINE dealloc (self)
  class(solver_t) :: self
  !.............................................................................
  call trace%begin ('solver_t%dealloc')
  call self%hd_t%dealloc
  call trace%end
END SUBROUTINE dealloc

!===============================================================================
!> Cast a generic task_t to patch_t
!===============================================================================
FUNCTION cast2solver (task) RESULT(solver)
  class(task_t), pointer:: task
  class(solver_t), pointer:: solver
  !.............................................................................
  select type (task)
  class is (solver_t)
  solver => task
  class default
  nullify(solver)
  call io%abort ('patch_t%cast: failed to cast a task to patch_t')
  end select
END FUNCTION cast2solver

!===============================================================================
!> Organize calls to the extras and the hydro solver
!===============================================================================
SUBROUTINE update (self)
  class(solver_t) :: self
  associate (d=>self%mem(:,:,:,self%idx%d,self%it,1))
  !.............................................................................
  call self%extras_t%pre_update
  call self%hd_t%update
  call self%extras_t%post_update
  end associate
END SUBROUTINE update

!===============================================================================
!> Organize calls to the extras and the hydro solver
!===============================================================================
SUBROUTINE dnload (self, only)
  class(solver_t) :: self
  integer, optional:: only
  !.............................................................................
  call self%extras_t%dnload (only)
END SUBROUTINE dnload

!===============================================================================
!> Get velocities from momenta
!===============================================================================
SUBROUTINE p2U(self, U, it)
  class(solver_t) :: self
  real, dimension(:,:,:,:), pointer:: U
  integer:: it
  !.............................................................................
  real, dimension(:,:,:,:), pointer:: dd
  integer:: i
  !-----------------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,           it,1), &
            p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1))
  do i=1,3
    U(:,:,:,i) = p(:,:,:,i)/d
  end do
  end associate
END SUBROUTINE p2U

!===============================================================================
!> Put velocities to momenta
!===============================================================================
SUBROUTINE U2p(self, U, it)
  class(solver_t) :: self
  real, dimension(:,:,:,:), pointer:: U
  integer:: it
  !.............................................................................
  real, dimension(:,:,:,:), pointer:: dd
  integer:: i
  !-----------------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d           ,it,1), &
            p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1))
  do i=1,3
    p(:,:,:,i) = U(:,:,:,i)*d
  end do
  end associate
END SUBROUTINE U2p

!===============================================================================
!> Get thermal energy per unit mass from entropy
!===============================================================================
SUBROUTINE e2E_th(self, E_th, it)
  class(solver_t) :: self
  integer:: it
  !.............................................................................
  real, dimension(:,:,:), pointer:: E_th
  real:: g1
  !-----------------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,it,1), &
            e => self%mem(:,:,:,self%idx%e,it,1))
  E_th = e/d
  end associate
END SUBROUTINE e2E_th

!===============================================================================
!> Put thermal energy to entropy per unit volume
!===============================================================================
SUBROUTINE E_th2e(self, E_th, it)
  class(solver_t) :: self
  integer:: it
  !.............................................................................
  real, dimension(:,:,:), pointer:: E_th
  real:: g1
  !-----------------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,it,1), &
            e => self%mem(:,:,:,self%idx%e,it,1))
  e = d*E_th
  end associate
END SUBROUTINE E_th2e

!===============================================================================
FUNCTION up(f,i) RESULT (g)
  real, dimension(:,:,:), intent(in):: f
  real, dimension(size(f,1),size(f,2),size(f,3)):: g
  integer:: i
  !.............................................................................
  g = 0.5*(cshift(f,1,i)+f)
END FUNCTION

!===============================================================================
FUNCTION gas_pressure_ngz (self) RESULT (pg)
  class(solver_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: pg
  integer :: l(3), u(3)
  !-----------------------------------------------------------------------------
  l = self%li; u = self%ui
  associate (d  => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%d ,self%it,1), &
             px => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%px,self%it,1), &
             py => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%py,self%it,1), &
             pz => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%pz,self%it,1), &
             e  => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e ,self%it,1))
  if (self%gamma==1.0) then
    pg = d(l(1):u(1),l(2):u(2),l(3):u(3))
  else
    pg = (self%gamma-1.0) &
       *(e(l(1):u(1),l(2):u(2),l(3):u(3)) -0.5*( &
           (px(l(1):u(1),l(2):u(2),l(3):u(3))**2 + &
            py(l(1):u(1),l(2):u(2),l(3):u(3))**2 + &
            pz(l(1):u(1),l(2):u(2),l(3):u(3))**2) &
             /d(l(1):u(1),l(2):u(2),l(3):u(3))))
  end if
  end associate
END FUNCTION gas_pressure_ngz

!===============================================================================
!> temperature from ideal gas law, taking into account that the temperature is
!> normalised by k_B and mu (so that specific gas constant = 1)
!===============================================================================
FUNCTION gas_temperature (self, lnd, ss) RESULT (tmp)
  class(solver_t):: self
  real, dimension(:,:,:), pointer:: lnd, ss
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: tmp
  !-----------------------------------------------------------------------------
  associate (d => self%mem(:,:,:,  1,self%it,1))
  tmp = self%gas_pressure()/d
  if (io%verbose > 2) then
    if (any(tmp<0)) then
      error stop "solver_t%gas_temperature: T<0"
    end if
  end if
  end associate
END FUNCTION gas_temperature

!===============================================================================
FUNCTION gas_temperature_ngz (self) RESULT (tmp)
  class(solver_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: tmp
  integer :: l(3), u(3)
  !---------------------------------------------------------------------
  l = self%li; u = self%ui
  associate (d => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),  1,self%it,1))
  tmp = self%gas_pressure_ngz()/d

  if (io%verbose > 2) then
    if (any(tmp<0)) then
      error stop "solver_t%gas_temperature: T<0"
    end if
  end if
  end associate
END FUNCTION gas_temperature_ngz

!===============================================================================
SUBROUTINE log_density (self, v)
  class(solver_t):: self
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = log(self%mem(:,:,:,self%id,self%it,1))
END SUBROUTINE

!===============================================================================
SUBROUTINE log_pressure (self, lnd, ss, v)
  class(solver_t):: self
  real, dimension(:,:,:), pointer:: lnd, ss, v
  !-----------------------------------------------------------------------------
  v = log(self%gas_pressure())
END SUBROUTINE

!===============================================================================
SUBROUTINE velocity_magnitude (self, v)
  class (solver_t):: self
  real, dimension(:,:,:), pointer:: v, d
  real, dimension(:,:,:,:), pointer:: p, U
  integer:: i
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d,self%it,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%px,self%it,1)
  v = 0.0
  do i=1,3
    v(:,:,:) = v(:,:,:) + (p(:,:,:,i)/d)**2
  end do
  v = sqrt(v)
END SUBROUTINE

!===============================================================================
SUBROUTINE void (self, v)
  class (solver_t):: self
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = 0.0
END SUBROUTINE

!===============================================================================
SUBROUTINE apply_heating (self, q)
  class (solver_t):: self
  real, dimension(:,:,:):: q
END SUBROUTINE

!===============================================================================
FUNCTION gas_velocity_vector (self) RESULT (v)
  class (solver_t):: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: d
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: p
  real(kind=KindScalarVar), dimension(self%gn(1),self%gn(2),self%gn(3),3):: v
  integer:: i
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d,self%it,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  do i=1,3
    v(:,:,:,i) = p(:,:,:,i)/d
  end do
END FUNCTION gas_velocity_vector

!===============================================================================
FUNCTION gas_velocity_scalar (self, idir) RESULT (v)
  class (solver_t):: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: d
  real(kind=KindScalarVar), dimension(self%gn(1),self%gn(2),self%gn(3)):: v
  integer:: idir
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d,self%it,1)
  select case (idir)
  case (1)
    v = self%mem(:,:,:,self%idx%px,self%it,1) / d
  case (2)
    v = self%mem(:,:,:,self%idx%py,self%it,1) / d
  case (3)
    v = self%mem(:,:,:,self%idx%pz,self%it,1) / d
  case default
    error stop "solver_mod::gas_velocity_scalar:: invalid value of idir"
  end select
END FUNCTION gas_velocity_scalar

!===============================================================================
!> Compute the compression
!===============================================================================
SUBROUTINE compression_magnitude (self, w)
  class (solver_t):: self
  real(kind=KindScalarVar), dimension(:,:,:):: w
  real(kind=KindScalarVar), dimension(:,:,:), allocatable:: vx, vy, vz
  integer:: i
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('solver_t%compression_magnitude', itimer=itimer)
  allocate (vx(self%gn(1),self%gn(2),self%gn(3)))
  allocate (vy(self%gn(1),self%gn(2),self%gn(3)))
  allocate (vz(self%gn(1),self%gn(2),self%gn(3)))
  vx = self%mem(:,:,:,self%idx%px,self%it,1)/self%mem(:,:,:,self%idx%d,self%it,1)
  vy = self%mem(:,:,:,self%idx%py,self%it,1)/self%mem(:,:,:,self%idx%d,self%it,1)
  vz = self%mem(:,:,:,self%idx%pz,self%it,1)/self%mem(:,:,:,self%idx%d,self%it,1)
  w = max(- stagger%ddx(self%ds(1),vx) &
          - stagger%ddy(self%ds(2),vy) &
          - stagger%ddz(self%ds(3),vz), 0.0)
  deallocate (vx, vy, vz)
  call trace%end (itimer)
END SUBROUTINE compression_magnitude

!===============================================================================
!> Compute the vorticity
!===============================================================================
SUBROUTINE vorticity_magnitude (self, w)
  class (solver_t):: self
  real(kind=KindScalarVar), dimension(:,:,:):: w
  real(kind=KindScalarVar), dimension(:,:,:), allocatable:: vx, vy, vz
  integer:: i
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('solver_t%vorticity_magnitude', itimer=itimer)
  allocate (vx(self%gn(1),self%gn(2),self%gn(3)))
  allocate (vy(self%gn(1),self%gn(2),self%gn(3)))
  allocate (vz(self%gn(1),self%gn(2),self%gn(3)))
  vx = self%mem(:,:,:,self%idx%px,self%it,1)/self%mem(:,:,:,self%idx%d,self%it,1)
  vy = self%mem(:,:,:,self%idx%py,self%it,1)/self%mem(:,:,:,self%idx%d,self%it,1)
  vz = self%mem(:,:,:,self%idx%pz,self%it,1)/self%mem(:,:,:,self%idx%d,self%it,1)
  w = sqrt((stagger%ddy(self%ds(2),vz)-stagger%ddz(self%ds(3),vy))**2 &
         + (stagger%ddz(self%ds(3),vx)-stagger%ddx(self%ds(1),vz))**2 &
         + (stagger%ddx(self%ds(1),vy)-stagger%ddy(self%ds(2),vx))**2)
  deallocate (vx, vy, vz)
  call trace%end (itimer)
END SUBROUTINE vorticity_magnitude

END MODULE solver_mod
