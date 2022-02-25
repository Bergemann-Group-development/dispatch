!===============================================================================
!> The solver_t data type extends the patch_t with a generic guard zone dnload
!> procedure, which covers the need for all simple, Cartesian mesh based 
!> solvers, and eliminates the need for a dnload procedure in each one of them.
!> It also adds restart functionality when initializing patches.
!===============================================================================
MODULE solver_mod
  USE io_mod
  USE trace_mod
  USE omp_lock_mod
  USE task_mod
  USE patch_mod
  USE kinds_mod
  USE mhd_mod
  USE vector_mod
  USE vector_ops
  USE validate_mod
  USE eos_mod
  USE shared_mod
  implicit none
  private
  type, public, extends(mhd_t):: solver_t
  contains
    procedure:: pre_init
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
    procedure:: gas_pressure_sub
    procedure:: gas_pressure_ngz
    !procedure:: cell_gas_pressure
    procedure:: gas_velocity_vector
    procedure:: gas_velocity_scalar
    procedure:: compression_magnitude
    procedure:: vorticity_magnitude
    generic:: gas_velocity => gas_velocity_vector, gas_velocity_scalar
    procedure:: to_thermal
    procedure:: from_generic
    procedure:: to_generic
  end type
  type(solver_t), public:: solver
CONTAINS

!===============================================================================
!> pre-initialize the hydro solver. 
!===============================================================================
SUBROUTINE pre_init (self)
  class(solver_t), target :: self
  !.............................................................................
  call trace%begin ('solver_t%pre_init')
  if (.not. self%pre_initialized) then
    call self%mhd_t%pre_init
    self%pre_initialized = .true.
  end if
  self%patch_ptr => self
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Organize calls to the extras and the hydro solver. 
!===============================================================================
SUBROUTINE init (self)
  class(solver_t) :: self
  !.............................................................................
  call trace%begin ('solver_t%init')
  call self%pre_init
  call self%mhd_t%init
  call self%from_generic
  call validate%init
  call trace%end ()
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
  call io%abort ('solver_t%cast2solver: failed to cast a task to solver_t')
  end select
END FUNCTION cast2solver

!===============================================================================
!> Organize calls to the extras and the hydro solver.  To be compatible with
!> Riemann solvers, which initially copy the %it slot to %new, and expect
!> pre_update() methods to refer to %new, we do the same copy operation here.
!===============================================================================
SUBROUTINE update (self)
  class(solver_t) :: self
  type(lock_t), pointer:: lock
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Make sure %new refers to the next timeslot, and assign pointers
  !-----------------------------------------------------------------------------
  call trace%begin ('solver_t%update', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Proceed, referring everywhere to the %new time slot
  !-----------------------------------------------------------------------------
  lock => self%mem_lock(self%new)
  call self%mhd_t%pre_update
  call self%mhd_t%update
  call self%mhd_t%post_update
  !-----------------------------------------------------------------------------
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Get velocities from momenta
!===============================================================================
SUBROUTINE p2U(self, U, it, linear)
  class(solver_t) :: self
  real, dimension(:,:,:,:), pointer, contiguous:: U
  integer:: it
  logical, optional:: linear
  !.....................................................................
  real, pointer, contiguous:: dd(:,:,:,:), d1(:,:,:)
  integer:: i, i3, i4
  integer, save:: itimer=0
  !---------------------------------------------------------------------
  call trace%begin ('solver_t%p2U', itimer=itimer)
  associate(d => self%mem(:,:,:,self%idx%d,           it,1), &
            p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1))
  allocate (dd, mold=p)
  allocate (d1, mold=d)
  if (self%kind(1:13) == 'stagger2e_pic') then
    d1 = log(d)
    dd = up(d1)
    do i4=1,3
      do i3=1,self%gn(3)
        U(:,:,i3,i4) = p(:,:,i3,i4)/exp(dd(:,:,i3,i4))
      end do
    end do
  else if (shared%option(linear,.false.)) then
    dd = down(d)
    do i4=1,3
      do i3=1,self%gn(3)
        U(:,:,i3,i4) = p(:,:,i3,i4)/dd(:,:,i3,i4)
      end do
    end do
  else
    d1 = log(d)
    dd=down(d1)
    do i4=1,3
      do i3=1,self%gn(3)
        U(:,:,i3,i4) = p(:,:,i3,i4)/exp(dd(:,:,i3,i4))
      end do
    end do
  end if
  deallocate (dd, d1)
  end associate
  call trace%end (itimer)
END SUBROUTINE p2U

!=======================================================================
!> Put velocities to momenta
!=======================================================================
SUBROUTINE U2p(self, U, it, linear)
  class(solver_t) :: self
  real, dimension(:,:,:,:), pointer:: U
  integer:: it
  logical, optional:: linear
  !.....................................................................
  real, pointer, contiguous:: dd(:,:,:,:), d1(:,:,:)
  integer:: i, i3, i4
  integer, save:: itimer=0
  !---------------------------------------------------------------------
  call trace%begin ('solver_t%U2p', itimer=itimer)
  associate(d => self%mem(:,:,:,self%idx%d           ,it,1), &
            p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1))
  allocate (dd, mold=U)
  allocate (d1, mold=d)
  if (self%kind(1:13) == 'stagger2e_pic') then
    d1 = log(d)
    dd = up(d1)
    do i4=1,3
      do i3=1,self%gn(3)
        p(:,:,i3,i4) = U(:,:,i3,i4)*exp(dd(:,:,i3,i4))
      end do
    end do
  else if (shared%option(linear,.false.)) then
    dd = down(d)
    do i4=1,3
      do i3=1,self%gn(3)
        p(:,:,i3,i4) = U(:,:,i3,i4)*dd(:,:,i3,i4)
      end do
    end do
  else
    d1 = log(d)
    dd = down(d1)
    do i4=1,3
      do i3=1,self%gn(3)
        p(:,:,i3,i4) = U(:,:,i3,i4)*exp(dd(:,:,i3,i4))
      end do
    end do
  end if
  deallocate (dd, d1)
  end associate
  call trace%end (itimer)
END SUBROUTINE U2p

!=======================================================================
!> Get thermal energy per unit mass from entropy
!=======================================================================
SUBROUTINE e2E_th(self, E_th, it)
  class(solver_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:), pointer:: E_th
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,it,1), &
            e => self%mem(:,:,:,self%idx%e,it,1))
  E_th = e/d
  end associate
END SUBROUTINE e2E_th

!=======================================================================
!> Put thermal energy to entropy per unit volume
!=======================================================================
SUBROUTINE E_th2e(self, E_th, it)
  class(solver_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:), pointer:: E_th
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,it,1), &
            e => self%mem(:,:,:,self%idx%e,it,1))
  e = d*E_th
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
  g1 = eos%gamma-1d0
  s = d*log(e*g1/d**eos%gamma)/g1
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
  g1 = eos%gamma-1d0
  ss = log(e*g1/d**eos%gamma)/g1
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
  g1 = eos%gamma-1d0
  if (self%idx%s > 0) then
    associate(d => self%mem(:,:,:,self%idx%d,it,1), &
              e => self%mem(:,:,:,self%idx%s,it,1))
    e = d**eos%gamma*exp(s/d*g1)/g1
    end associate
  else
    associate(d => self%mem(:,:,:,self%idx%d,it,1), &
              e => self%mem(:,:,:,self%idx%e,it,1))
    e = d**eos%gamma*exp(s/d*g1)/g1
    end associate
  end if
END SUBROUTINE s2e

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
  g1 = eos%gamma-1d0
  if (self%idx%s > 0) then
    associate(d => self%mem(:,:,:,self%idx%d,it,1), &
              e => self%mem(:,:,:,self%idx%s,it,1))
    e = d**eos%gamma*exp(ss*g1)/g1
    end associate
  else
    associate(d => self%mem(:,:,:,self%idx%d,it,1), &
              e => self%mem(:,:,:,self%idx%e,it,1))
    e = d**eos%gamma*exp(ss*g1)/g1
    end associate
  end if
END SUBROUTINE ss2e

!===============================================================================
!> temperature from ideal gas law, taking into account that the temperature is
!> normalised by k_B and mu (so that specific gas constant = 1)
!===============================================================================
!> Subroutine interface to solver-independent gas pressure function
!===============================================================================
SUBROUTINE gas_pressure_sub (self, pg)
  class(solver_t):: self
  real, dimension(:,:,:), pointer, contiguous:: d, e, pg
  !-----------------------------------------------------------------------------
  d  => self%mem(:,:,:,self%idx%d ,self%it,1)
  e  => self%mem(:,:,:,self%idx%e ,self%it,1)
  if (eos%gamma==1d0) then
    pg = d*self%csound**2
  else
    call eos%lookup (d=d, e=e, pg=pg)
  end if
END SUBROUTINE gas_pressure_sub

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
             bx => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%bx,self%it,1), &
             by => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%by,self%it,1), &
             bz => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%bz,self%it,1), &
             e  => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e ,self%it,1))
  if (eos%gamma==1.0) then
    pg = d(l(1):u(1),l(2):u(2),l(3):u(3))
  else
    pg = (eos%gamma-1d0) *e(l(1):u(1),l(2):u(2),l(3):u(3))
  end if
  end associate
END FUNCTION gas_pressure_ngz

!===============================================================================
!> Subroutine interface to solver-independent gas pressure function
!===============================================================================
SUBROUTINE gas_temperature_sub (self, tt)
  class(solver_t):: self
  real, dimension(:,:,:), pointer, contiguous:: d, px, py, pz, e, eth, tt
  real, dimension(:,:,:), pointer, contiguous:: bx, by, bz
  integer:: m(3)
  !-----------------------------------------------------------------------------
  d  => self%mem(:,:,:,self%idx%d ,self%it,1)
  px => self%mem(:,:,:,self%idx%px,self%it,1)
  py => self%mem(:,:,:,self%idx%py,self%it,1)
  pz => self%mem(:,:,:,self%idx%pz,self%it,1)
  e  => self%mem(:,:,:,self%idx%e ,self%it,1)
  if (eos%gamma==1d0) then
    tt = d*self%csound**2
  else
    call eos%lookup (d=d, e=e, tt=tt)
  end if
END SUBROUTINE gas_temperature_sub

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
  contiguous:: v
  !-----------------------------------------------------------------------------
  call self%gas_pressure_sub (v)
  v = log(v)
END SUBROUTINE

!===============================================================================
SUBROUTINE velocity_magnitude (self, v)
  class (solver_t):: self
  real, dimension(:,:,:), pointer:: v, d
  real, dimension(:,:,:,:), pointer:: p, U,dd
  integer:: i
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d,self%it,1)
  dd=down(d)
  p => self%mem(:,:,:,self%idx%px:self%idx%px,self%it,1)
  v = 0.0
  do i=1,3
    v(:,:,:) = v(:,:,:) + (p(:,:,:,i)/dd(:,:,:,i))**2
  end do
  v = sqrt(v)
END SUBROUTINE

!===============================================================================
SUBROUTINE magnetic_field_magnitude (self, v)
  class (solver_t):: self
  real, dimension(:,:,:), pointer:: v
  real, dimension(:,:,:,:), pointer:: b
  !-----------------------------------------------------------------------------
  associate (bx => self%mem(:,:,:,self%idx%bx,self%it,1), &
             by => self%mem(:,:,:,self%idx%by,self%it,1), &
             bz => self%mem(:,:,:,self%idx%bz,self%it,1))
  v = xup(bx)**2 + yup(by)**2 + zup(bz)**2
  v = sqrt(v)
  end associate
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
    w = max(- ddxdn(self%ds,vx) &
            - ddydn(self%ds,vy) &
            - ddzdn(self%ds,vz), 0.0)
  else
    w = max(- ddxup(self%ds,vx) &
            - ddyup(self%ds,vy) &
            - ddzup(self%ds,vz), 0.0)
  end if
  deallocate (vx, vy, vz, d)
  call trace%end (itimer)
END SUBROUTINE compression_magnitude

!===============================================================================
!> Compute the vorticity
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
    w = sqrt(ydn(zdn(ddyup(self%ds,vz)-ddzup(self%ds,vy)))**2 &
           + zdn(xdn(ddzup(self%ds,vx)-ddxup(self%ds,vz)))**2 &
           + xdn(ydn(ddxup(self%ds,vy)-ddyup(self%ds,vx)))**2)
  else
    w = sqrt(yup(zup(ddydn(self%ds,vz)-ddzdn(self%ds,vy)))**2 &
           + zup(xup(ddzdn(self%ds,vx)-ddxdn(self%ds,vz)))**2 &
           + xup(yup(ddxdn(self%ds,vy)-ddydn(self%ds,vx)))**2)
  end if
  deallocate (vx, vy, vz)
  call trace%end (itimer)
END SUBROUTINE vorticity_magnitude

!===============================================================================
!> Compute or assign thermal energy
!===============================================================================
SUBROUTINE to_thermal (self, it, eth, kind)
  class(solver_t):: self
  integer:: it
  real, contiguous:: eth(:,:,:)
  character(len=*), optional:: kind
  !.............................................................................
  integer:: i1, i2, i3
  real, pointer, contiguous, dimension(:,:,:):: e
  !-----------------------------------------------------------------------------
  call trace%begin ('solver_t%to_thermal')
  e  => self%mem(:,:,:,self%idx%e ,it,1)
  eth = e
  call trace%end ()
END SUBROUTINE to_thermal

!===============================================================================
!> Convert from rho,T,p,B to solver-preferred variables
!===============================================================================
SUBROUTINE from_generic (self)
  class(solver_t):: self
  real(var_k), pointer, dimension(:,:,:):: dm, db, em, Tg, ss
  integer:: iv, l(3), u(3), ie
  !-----------------------------------------------------------------------------
  call trace%begin ('solver_t%from_generic')
  if (associated(self%generic)) then
    l = merge (self%mesh%lb, self%mesh%li, io%guard_zones)
    u = merge (self%mesh%ub, self%mesh%ui, io%guard_zones)
    ! Density
    db => self%generic(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%d)
    dm => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%d,self%it,1)
    dm = db
    ! Momentum and magnetic fields
    do iv=self%idx%px,self%idx%pz
      self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,self%it,1) = &
      self%generic(l(1):u(1),l(2):u(2),l(3):u(3),iv) *dm
    end do
    if (self%idx%bx > 0) then
      do iv=self%idx%bx,self%idx%bz
        self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,self%it,1) = &
        self%generic(l(1):u(1),l(2):u(2),l(3):u(3),iv)
      end do
    end if
    ! Energy or entropy
    select case(trim(self%etype))
    case('thermal')
      Tg => self%generic(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e)
      em => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e,self%it,1)
      call eos%inverse_lookup (d=dm, v=Tg, name='T', ee=em)
      em = em*dm
    case('entropy')
      Tg => self%generic(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%s)
      ss => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%s,self%it,1)
      call eos%inverse_lookup (d=dm, v=Tg, name='T', ss=ss)
      ss = ss * dm
    case default
      call io%abort("unknown etype in solver_mod::from_generic")
    end select
  end if
  call trace%end ()
END SUBROUTINE from_generic

!===============================================================================
!> Convert from rho,T,p,B to solver-preferred variables
!===============================================================================
SUBROUTINE to_generic (self)
  class(solver_t):: self
  real(var_k), pointer, contiguous, dimension(:,:,:):: dm, db, em, Tg, s
  integer:: iv, l(3), u(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('solver_t%to_generic')
  if (associated(self%generic)) then
    l = merge (self%mesh%lb, self%mesh%li, io%guard_zones)
    u = merge (self%mesh%ub, self%mesh%ui, io%guard_zones)
    ! Density
    db => self%generic(:,:,:,self%idx%d)
    dm => self%mem(:,:,:,self%idx%d,self%it,1)
    em => self%mem(:,:,:,self%idx%e,self%it,1)
    db = dm
    ! Momentum and magnetic fields
    do iv=self%idx%px,self%idx%pz
      self%generic(:,:,:,iv) = self%mem(:,:,:,iv,self%it,1)
    end do
    if (self%idx%bx > 0) then
      do iv=self%idx%bx,self%idx%bz
        self%generic(:,:,:,iv) = self%mem(:,:,:,iv,self%it,1)
      end do
    end if
    ! Energy or entropy
    select case(trim(self%etype))
    case('thermal')
      Tg => self%generic(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e)
      em => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e,self%it,1)
      call eos%lookup (d=db, e=em, tt=Tg)
    case('entropy')
      Tg => self%generic(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%s)
      s  => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%s,self%it,1)
      call eos%lookup (d=db, s=s, tt=Tg)
    case default
      call io%abort("unknown etype in solver_mod::from_generic")
    end select
  end if
  call trace%end ()
END SUBROUTINE to_generic

END MODULE solver_mod
