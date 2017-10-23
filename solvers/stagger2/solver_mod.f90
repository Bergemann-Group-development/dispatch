!> =============================================================================
!> The solver_t data type extends the patch_t with a generic guard zone dnload
!> procedure, which covers the need for all simple, Cartesian mesh based 
!> solvers, and eliminates the need for a dnload procedure in each one of them.
!> It also adds restart functionality when initializing patches.
!> =============================================================================
MODULE solver_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE kinds_mod
  USE mhd_mod
  USE vector_mod
  USE vector_ops
  implicit none
  private
  type, public, extends(mhd_t):: solver_t
  contains
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
    procedure:: gas_pressure_ngz
    procedure:: cell_gas_pressure
    procedure:: gas_temperature
    procedure:: gas_temperature_ngz
    procedure:: cell_gas_temperature
    procedure:: velocity_magnitude => void
    procedure:: magnetic_field_magnitude => void
    procedure:: grav_potential => void
    procedure:: apply_heating
    procedure:: gas_velocity_vector
    procedure:: gas_velocity_scalar
    generic:: gas_velocity => gas_velocity_vector, gas_velocity_scalar
  end type
  type(solver_t):: solver
CONTAINS

!> =====================================================================
!> Get velocities from momenta
!> =====================================================================
SUBROUTINE p2U(self, U, it)
  class(solver_t) :: self
  real, dimension(:,:,:,:), pointer:: U
  integer:: it
  !.....................................................................
  real, dimension(:,:,:,:), pointer:: dd
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,           it,1), &
            p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1))
  allocate(dd(size(d,1),size(d,2),size(d,3),3))
  dd=down(log(d))
  U = p/exp(dd)
  deallocate(dd)
  end associate
END SUBROUTINE p2U

!> =====================================================================
!> Put velocities to momenta
!> =====================================================================
SUBROUTINE U2p(self, U, it)
  class(solver_t) :: self
  real, dimension(:,:,:,:), pointer:: U
  integer:: it
  !.....................................................................
  real, dimension(:,:,:,:), pointer:: dd
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d           ,it,1), &
            p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1))
  allocate(dd(size(d,1),size(d,2),size(d,3),3))
  dd=down(log(d))
  p = U*exp(dd)
  deallocate(dd)
  end associate
END SUBROUTINE U2p

!> =====================================================================
!> Get thermal energy per unit mass from entropy
!> =====================================================================
SUBROUTINE e2E_th(self, E_th, it)
  class(solver_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:), pointer:: E_th
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,it,1), &
            s => self%mem(:,:,:,self%idx%s,it,1))
  g1 = self%gamma-1.0
  E_th = d**g1*exp(s*g1/d)/g1
  end associate
END SUBROUTINE e2E_th

!> =====================================================================
!> Put thermal energy to entropy per unit volume
!> =====================================================================
SUBROUTINE E_th2e(self, E_th, it)
  class(solver_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:), pointer:: E_th
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,it,1), &
            s => self%mem(:,:,:,self%idx%s,it,1))
  g1 = self%gamma-1.0
  s = d*log(E_th*g1/d**g1)/g1
  end associate
END SUBROUTINE E_th2e

!> =====================================================================
!> Get entropy per unit volume from energy per unit volume
!> =====================================================================
SUBROUTINE e2s(self, s, it)
  class(solver_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:):: s
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,it,1), &
            e => self%mem(:,:,:,self%idx%s,it,1))
  g1 = self%gamma-1.0
  s = d*log(e*g1/d**self%gamma)/g1
  end associate
END SUBROUTINE e2s

!> =====================================================================
!> Get entropy per unit mass from energy per unit volume
!> =====================================================================
SUBROUTINE e2ss(self, ss, it)
  class(solver_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:):: ss
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,it,1), &
            e => self%mem(:,:,:,self%idx%s,it,1))
  g1 = self%gamma-1.0
  ss = log(e*g1/d**self%gamma)/g1
  end associate
END SUBROUTINE e2ss

!> =====================================================================
!> Get thermal energy per unit volume from entropy per unit volume
!> =====================================================================
SUBROUTINE s2e(self, s, it)
  class(solver_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:):: s
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,it,1), &
            e => self%mem(:,:,:,self%idx%s,it,1))
  g1 = self%gamma-1.0
  e = d**self%gamma*exp(s/d*g1)/g1
  end associate
END SUBROUTINE S2e

!> =====================================================================
!> Get thermal energy per unit volume from entropy per unit mass
!> =====================================================================
SUBROUTINE ss2e(self, ss, it)
  class(solver_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:):: ss
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,it,1), &
            e => self%mem(:,:,:,self%idx%s,it,1))
  g1 = self%gamma-1.0
  e = d**self%gamma*exp(ss/d*g1)/g1
  end associate
END SUBROUTINE ss2e

!> =============================================================================
FUNCTION gas_pressure_ngz (self) RESULT (pg)
  class(solver_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: pg
  integer :: l(3), u(3)
  !-----------------------------------------------------------------------------
  l = self%li; u = self%ui
  associate (d => self%mem(l(1):u(1),l(2):u(2),l(3):u(3), 1,self%it,1), &
             s => self%mem(l(1):u(1),l(2):u(2),l(3):u(3), 2,self%it,1))
  if (self%gamma==1.0) then
    pg = d
  else
    pg = d**self%gamma*exp(s/d*(self%gamma-1))
  end if
  end associate
END FUNCTION gas_pressure_ngz

!> =============================================================================
FUNCTION cell_gas_pressure(self, c) RESULT (pg)
  class(solver_t) :: self
  real :: pg
  integer :: c(3)
  !---------------------------------------------------------------------
  associate (d => self%mem(c(1),c(2),c(3), 1,self%it,1), &
             s => self%mem(c(1),c(2),c(3), 2,self%it,1))
  if (self%gamma==1.0) then
    pg = d
  else
    pg = d**self%gamma*exp(s/d*(self%gamma-1))
  end if
  end associate
END FUNCTION cell_gas_pressure

!> =============================================================================
!> temperature from ideal gas law, taking into account that the temperature is
!> normalised by k_B and mu (so that specific gas constant = 1)
!> =============================================================================
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

!> =============================================================================
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

FUNCTION cell_gas_temperature (self, c) RESULT (tmp)
  class(solver_t):: self
  real :: tmp
  integer :: c(3)
  !-----------------------------------------------------------------------------
  associate (d => self%mem(c(1),c(2),c(3),  1,self%it,1))
  tmp = cell_gas_pressure(self,c)/d
  end associate
END FUNCTION cell_gas_temperature

!> =============================================================================
SUBROUTINE log_density (self, v)
  class(solver_t):: self
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = log(self%mem(:,:,:,self%id,self%it,1))
END SUBROUTINE

!> =============================================================================
SUBROUTINE log_pressure (self, lnd, ss, v)
  class(solver_t):: self
  real, dimension(:,:,:), pointer:: lnd, ss, v
  !-----------------------------------------------------------------------------
  v = log(self%gas_pressure())
END SUBROUTINE

!> =============================================================================
SUBROUTINE velocity_magnitude (self, v)
  class (solver_t):: self
  real, dimension(:,:,:), pointer:: v
  real, dimension(:,:,:,:), pointer:: p
  !-----------------------------------------------------------------------------
  p => self%mem(:,:,:,self%idx%px:self%idx%px,self%it,1)
  v = norm(p/exp(down(log(self%mem(:,:,:,self%id,self%it,1)))))
END SUBROUTINE

!> =============================================================================
SUBROUTINE magnetic_field_magnitude (self, v)
  class (solver_t):: self
  real, dimension(:,:,:), pointer:: v
  real, dimension(:,:,:,:), pointer:: b
  !-----------------------------------------------------------------------------
  b => self%mem(:,:,:,self%idx%bx:self%idx%bx,self%it,1)
  v = norm(b)
END SUBROUTINE

!> =============================================================================
SUBROUTINE void (self, v)
  class (solver_t):: self
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = 0.0
END SUBROUTINE

!> =============================================================================
SUBROUTINE apply_heating (self, q)
  class (solver_t):: self
  real, dimension(:,:,:):: q
END SUBROUTINE

!> =============================================================================
FUNCTION gas_velocity_vector (self) RESULT (v)
  USE scalar_mod
  class (solver_t):: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: d, lnd
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: p, ld, dd
  real(kind=KindScalarVar), dimension(self%gn(1),self%gn(2),self%gn(3),3):: v
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d,self%it,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  call allocate_vectors (self%gn, ld, dd)
  call allocate_scalars (self%gn, lnd)
  lnd = log(d)
  ld = down(lnd)
  dd = exp(ld)
  v = p/dd

  call deallocate_vectors (ld, dd)
  call deallocate_scalars (lnd)
  nullify (d,p)
END FUNCTION gas_velocity_vector

!> =============================================================================
FUNCTION gas_velocity_scalar (self, idir) RESULT (v)
  USE scalar_mod
  class (solver_t):: self
  integer:: idir
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: d, lnd, ld, dd
  real(kind=KindScalarVar), dimension(self%gn(1),self%gn(2),self%gn(3)):: v
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d,self%it,1)
  call allocate_scalars (self%gn, lnd, ld, dd)
  lnd = log(d)
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
    error stop "solver_mod::gas_velocity_scalar:: invalid value of idir"
  end select

  call deallocate_scalars (lnd, ld, dd)
  nullify (d)
END FUNCTION gas_velocity_scalar

END MODULE solver_mod
