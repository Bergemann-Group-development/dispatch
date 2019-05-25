!===============================================================================
!> The solver_t data type extends the patch_t with a generic guard zone dnload
!> procedure, which covers the need for all simple, Cartesian mesh based 
!> solvers, and eliminates the need for a dnload procedure in each one of them.
!> It also adds restart functionality when initializing patches.
!===============================================================================
MODULE solver_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE kinds_mod
  USE mhd_mod
  USE vector_mod
  USE vector_ops
  USE validate_mod
  implicit none
  private
  type, public, extends(mhd_t):: solver_t
  contains
    procedure:: init
    procedure:: update
    procedure, nopass:: cast2solver
    procedure:: p2U
    procedure:: U2p
    procedure:: e2E_th
    procedure:: E_th2e
    procedure:: e2s
    procedure:: s2e
    procedure:: e2ss
    procedure:: ss2e
    procedure:: log_density => void
    procedure:: log_pressure => void
    !procedure:: gas_pressure_ngz
    !procedure:: cell_gas_pressure
    procedure:: gas_temperature
    !procedure:: gas_temperature_ngz
    !procedure:: cell_gas_temperature
    procedure:: velocity_magnitude => void
    procedure:: magnetic_field_magnitude => void
    procedure:: grav_potential => void
    procedure:: apply_heating
    procedure:: compression_magnitude
    procedure:: vorticity_magnitude
    procedure:: gas_velocity_vector
    procedure:: gas_velocity_scalar
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
  call self%mhd_t%init
  call self%extras_t%init
  call validate%init
END SUBROUTINE init

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
  call validate%check (self, d, 'before update')
  call self%mhd_t%update
  call validate%check (self, d, ' after update')
  call self%extras_t%post_update
  end associate
END SUBROUTINE update

!===============================================================================
!> Get velocities from momenta
!===============================================================================
SUBROUTINE p2U(self, U, it)
  class(solver_t) :: self
  real, dimension(:,:,:,:), pointer:: U
  integer:: it
  !.............................................................................
  real, dimension(:,:,:,:), allocatable:: dd
  !-----------------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,           it,1), &
            p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1))
  allocate(dd(size(d,1),size(d,2),size(d,3),3))
  if (self%kind(1:13) == 'stagger2e_pic') then
    dd=up(log(d))
    U = p/exp(dd)
  else
    dd=down(log(d))
    U = p/exp(dd)
  end if
  deallocate(dd)
  end associate
END SUBROUTINE p2U

!=======================================================================
!> Put velocities to momenta
!=======================================================================
SUBROUTINE U2p(self, U, it)
  class(solver_t) :: self
  real, dimension(:,:,:,:), pointer:: U
  integer:: it
  !.............................................................................
  real, dimension(:,:,:,:), allocatable:: dd
  !-----------------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d           ,it,1), &
            p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1))
  allocate(dd(size(d,1),size(d,2),size(d,3),3))
  if (self%kind(1:13) == 'stagger2e_pic') then
    dd=up(log(d))
    p = U*exp(dd)
  else
    dd=down(log(d))
    p = U*exp(dd)
  end if
  deallocate(dd)
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
            s => self%mem(:,:,:,self%idx%s,it,1))
  g1 = self%gamma-1.0
  E_th = d**g1*exp(s*g1/d)/g1
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
            s => self%mem(:,:,:,self%idx%s,it,1))
  g1 = self%gamma-1.0
  s = d*log(E_th*g1/d**g1)/g1
  end associate
END SUBROUTINE E_th2e

!===============================================================================
!> Get entropy per unit volume from energy per unit volume
!===============================================================================
SUBROUTINE e2s(self, s, it)
  class(solver_t) :: self
  integer:: it
  !.............................................................................
  real, dimension(:,:,:):: s
  real:: g1
  !-----------------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,it,1), &
            e => self%mem(:,:,:,self%idx%s,it,1))
  g1 = self%gamma-1.0
  s = d*log(e*g1/d**self%gamma)/g1
  end associate
END SUBROUTINE e2s

!===============================================================================
!> Get entropy per unit mass from energy per unit volume
!===============================================================================
SUBROUTINE e2ss(self, ss, it)
  class(solver_t) :: self
  integer:: it
  !.............................................................................
  real, dimension(:,:,:):: ss
  real:: g1
  !-----------------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,it,1), &
            e => self%mem(:,:,:,self%idx%s,it,1))
  g1 = self%gamma-1.0
  ss = log(e*g1/d**self%gamma)/g1
  end associate
END SUBROUTINE e2ss

!===============================================================================
!> Get thermal energy per unit volume from entropy per unit volume
!===============================================================================
SUBROUTINE s2e(self, s, it)
  class(solver_t) :: self
  integer:: it
  !.............................................................................
  real, dimension(:,:,:):: s
  real:: g1
  !-----------------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,it,1), &
            e => self%mem(:,:,:,self%idx%s,it,1))
  g1 = self%gamma-1.0
  e = d**self%gamma*exp(s/d*g1)/g1
  end associate
END SUBROUTINE S2e

!===============================================================================
!> Get thermal energy per unit volume from entropy per unit mass
!===============================================================================
SUBROUTINE ss2e(self, ss, it)
  class(solver_t) :: self
  integer:: it
  !.............................................................................
  real, dimension(:,:,:):: ss
  real:: g1
  !-----------------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,it,1), &
            e => self%mem(:,:,:,self%idx%s,it,1))
  g1 = self%gamma-1.0
  e = d**self%gamma*exp(ss*g1)/g1
  end associate
END SUBROUTINE ss2e

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
      call io%abort ("solver_t%gas_temperature: T<0")
    end if
  end if
  end associate
END FUNCTION gas_temperature

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
  real, dimension(:,:,:), pointer:: v
  real, dimension(:,:,:,:), pointer:: p
  !-----------------------------------------------------------------------------
  p => self%mem(:,:,:,self%idx%px:self%idx%px,self%it,1)
  if (self%kind(1:13) == 'stagger2e_pic') then
    v = norm(p/exp(up(log(self%mem(:,:,:,self%id,self%it,1)))))
  else
    v = norm(p/exp(down(log(self%mem(:,:,:,self%id,self%it,1)))))
  end if
END SUBROUTINE

!===============================================================================
SUBROUTINE magnetic_field_magnitude (self, v)
  class (solver_t):: self
  real, dimension(:,:,:), pointer:: v
  real, dimension(:,:,:,:), pointer:: b
  !-----------------------------------------------------------------------------
  b => self%mem(:,:,:,self%idx%bx:self%idx%bx,self%it,1)
  v = norm(b)
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
  USE scalar_mod
  class (solver_t):: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: d
  real(kind=KindScalarVar), dimension(:,:,:), allocatable :: lnd
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: p
  real(kind=KindScalarVar), dimension(:,:,:,:), allocatable:: ld, dd
  real(kind=KindScalarVar), dimension(self%gn(1),self%gn(2),self%gn(3),3):: v
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d,self%it,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  call allocate_vectors_a (self%gn, ld, dd)
  call allocate_scalars_a (self%gn, lnd)
  if (self%kind(1:13) == 'stagger2e_pic') then
    lnd = log(d)
    ld = up(lnd)
    dd = exp(ld)
    v = p/dd
  else
    lnd = log(d)
    ld = down(lnd)
    dd = exp(ld)
    v = p/dd
  end if
  call deallocate_vectors_a (ld, dd)
  call deallocate_scalars_a (lnd)
  nullify (d,p)
END FUNCTION gas_velocity_vector

!===============================================================================
FUNCTION gas_velocity_scalar (self, idir) RESULT (v)
  USE scalar_mod
  class (solver_t):: self
  integer:: idir
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d
  real(kind=KindScalarVar), dimension(:,:,:), allocatable:: lnd, ld, dd
  real(kind=KindScalarVar), dimension(self%gn(1),self%gn(2),self%gn(3)):: v
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d,self%it,1)
  call allocate_scalars_a (self%gn, lnd, ld, dd)
  lnd = log(d)
  if (self%kind(1:13) == 'stagger2e_pic') then
    select case (idir)
    case (1)
      ld = xup(lnd)
      dd = exp(ld)
      v = self%mem(:,:,:,self%idx%px,self%it,1) / dd
    case (2)
      ld = yup(lnd)
      dd = exp(ld)
      v = self%mem(:,:,:,self%idx%py,self%it,1) / dd
    case (3)
      ld = zup(lnd)
      dd = exp(ld)
      v = self%mem(:,:,:,self%idx%pz,self%it,1) / dd
    case default
      call io%abort ("solver_mod::gas_velocity_scalar:: invalid value of idir")
    end select
  else
    select case (idir)
    case (1)
      ld = xdn(lnd)
      dd = exp(ld)
      v = self%mem(:,:,:,self%idx%px,self%it,1) / dd
    case (2)
      ld = ydn(lnd)
      dd = exp(ld)
      v = self%mem(:,:,:,self%idx%py,self%it,1) / dd
    case (3)
      ld = zdn(lnd)
      dd = exp(ld)
      v = self%mem(:,:,:,self%idx%pz,self%it,1) / dd
    case default
      call io%abort ("solver_mod::gas_velocity_scalar:: invalid value of idir")
    end select
  end if
  call deallocate_scalars_a (lnd, ld, dd)
  nullify (d)
END FUNCTION gas_velocity_scalar

!===============================================================================
!> Compute the compression, centered in cells, given that momenta are dn-staggered
!> in stagger2, except for stagger2e_pic, where it up-staggered
!===============================================================================
SUBROUTINE compression_magnitude (self, w)
  class (solver_t):: self
  real(kind=KindScalarVar), dimension(:,:,:):: w
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: d
  real(kind=KindScalarVar), dimension(:,:,:), allocatable:: vx, vy, vz
  integer:: i
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('solver_t%compression_magnitude', itimer=itimer)
  allocate (vx(self%gn(1),self%gn(2),self%gn(3)))
  allocate (vy(self%gn(1),self%gn(2),self%gn(3)))
  allocate (vz(self%gn(1),self%gn(2),self%gn(3)))
  allocate ( d(self%gn(1),self%gn(2),self%gn(3),3))
  if (self%kind(1:13) == 'stagger2e_pic') then
    d = up(log(self%mem(:,:,:,self%idx%d,self%it,1)))
  else
    d = down(log(self%mem(:,:,:,self%idx%d,self%it,1)))
  end if
  vx = self%mem(:,:,:,self%idx%px,self%it,1)/exp(d(:,:,:,1))
  vy = self%mem(:,:,:,self%idx%py,self%it,1)/exp(d(:,:,:,2))
  vz = self%mem(:,:,:,self%idx%pz,self%it,1)/exp(d(:,:,:,3))
  if (self%kind(1:13) == 'stagger2e_pic') then
    w = max(- ddxdn(self%ds(1),vx) &
            - ddydn(self%ds(2),vy) &
            - ddzdn(self%ds(3),vz), 0.0)
  else
    w = max(- ddxup(self%ds(1),vx) &
            - ddyup(self%ds(2),vy) &
            - ddzup(self%ds(3),vz), 0.0)
  end if
  deallocate (vx, vy, vz, d)
  call trace%end (itimer)
END SUBROUTINE compression_magnitude

!===============================================================================
!> Compute the vorticity, centered in cells, given that momenta are dn-staggered
!> in stagger2, except for stagger2e_pic, where it up-staggered
!===============================================================================
SUBROUTINE vorticity_magnitude (self, w)
  class (solver_t):: self
  real(kind=KindScalarVar), dimension(:,:,:):: w
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: d
  real(kind=KindScalarVar), dimension(:,:,:), allocatable:: vx, vy, vz
  integer:: i
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('solver_t%vorticity_magnitude', itimer=itimer)
  allocate (vx(self%gn(1),self%gn(2),self%gn(3)))
  allocate (vy(self%gn(1),self%gn(2),self%gn(3)))
  allocate (vz(self%gn(1),self%gn(2),self%gn(3)))
  allocate ( d(self%gn(1),self%gn(2),self%gn(3),3))
  if (self%kind(1:13) == 'stagger2e_pic') then
    d = up(log(self%mem(:,:,:,self%idx%d,self%it,1)))
  else
    d = down(log(self%mem(:,:,:,self%idx%d,self%it,1)))
  end if
  vx = self%mem(:,:,:,self%idx%px,self%it,1)/exp(d(:,:,:,1))
  vy = self%mem(:,:,:,self%idx%py,self%it,1)/exp(d(:,:,:,2))
  vz = self%mem(:,:,:,self%idx%pz,self%it,1)/exp(d(:,:,:,3))
  if (self%kind(1:13) == 'stagger2e_pic') then
    w = sqrt(ydn(zdn(ddyup(self%ds(2),vz)-ddzup(self%ds(3),vy)))**2 &
           + zdn(xdn(ddzup(self%ds(3),vx)-ddxup(self%ds(1),vz)))**2 &
           + xdn(ydn(ddxup(self%ds(1),vy)-ddyup(self%ds(2),vx)))**2)
  else
    w = sqrt(yup(zup(ddydn(self%ds(2),vz)-ddzdn(self%ds(3),vy)))**2 &
           + zup(xup(ddzdn(self%ds(3),vx)-ddxdn(self%ds(1),vz)))**2 &
           + xup(yup(ddxdn(self%ds(1),vy)-ddydn(self%ds(2),vx)))**2)
  end if
  deallocate (vx, vy, vz)
  call trace%end (itimer)
END SUBROUTINE vorticity_magnitude

END MODULE solver_mod
