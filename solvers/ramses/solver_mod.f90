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
    procedure:: gas_velocity
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
  integer:: i
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%i_d,           it,1), &
            p => self%mem(:,:,:,self%i_px:self%i_pz,it,1))
  do i=1,3
    U(:,:,:,i) = p(:,:,:,i)/d
  end do
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
  integer:: i
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%i_d           ,it,1), &
            p => self%mem(:,:,:,self%i_px:self%i_pz,it,1))
  do i=1,3
    p(:,:,:,i) = U(:,:,:,i)*d
  end do
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
  associate(d => self%mem(:,:,:,self%i_d,it,1), &
            s => self%mem(:,:,:,self%i_s,it,1))
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
  associate(d => self%mem(:,:,:,self%i_d,it,1), &
            s => self%mem(:,:,:,self%i_s,it,1))
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
  associate(d => self%mem(:,:,:,self%i_d,it,1), &
            e => self%mem(:,:,:,self%i_s,it,1))
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
  associate(d => self%mem(:,:,:,self%i_d,it,1), &
            e => self%mem(:,:,:,self%i_s,it,1))
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
  associate(d => self%mem(:,:,:,self%i_d,it,1), &
            e => self%mem(:,:,:,self%i_s,it,1))
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
  associate(d => self%mem(:,:,:,self%i_d,it,1), &
            e => self%mem(:,:,:,self%i_s,it,1))
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
  real, dimension(:,:,:), pointer:: v, d
  real, dimension(:,:,:,:), pointer:: p, U
  integer:: i
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,self%i_d,self%it,1)
  p => self%mem(:,:,:,self%i_px:self%i_px,self%it,1)
  v = 0.0
  do i=1,3
    v(:,:,:) = v(:,:,:) + (p(:,:,:,i)/d)**2
  end do
  v = sqrt(v)
END SUBROUTINE

!> =============================================================================
SUBROUTINE magnetic_field_magnitude (self, v)
  class (solver_t):: self
  real, dimension(:,:,:), pointer:: v
  real, dimension(:,:,:,:), pointer:: b
  integer:: i
  !-----------------------------------------------------------------------------
  b => self%mem(:,:,:,self%i_bx:self%i_bx,self%it,1)
  v = 0.0
  do i=1,3
    v(:,:,:) = v(:,:,:) + b(:,:,:,i)**2
  end do
  v = sqrt(v)
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
FUNCTION gas_velocity (self) RESULT (v)
  class (solver_t):: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: d
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: p
  real(kind=KindScalarVar), dimension(self%gn(1),self%gn(2),self%gn(3),3):: v
  integer:: i
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,self%i_d,self%it,1)
  p => self%mem(:,:,:,self%i_px:self%i_pz,self%it,1)
  do i=1,3
    v(:,:,:,i) = p(:,:,:,i)/d
  end do
END FUNCTION gas_velocity

END MODULE solver_mod
