!===============================================================================
!> The solver_t data type extends the patch_t with a generic guard zone dnload
!> procedure, which covers the need for all simple, Cartesian mesh based 
!> solvers, and eliminates the need for a dnload procedure in each one of them.
!> It also adds restart functionality when initializing patches.
!===============================================================================
MODULE solver_mod
  USE io_mod
  USE omp_lock_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE kinds_mod
  USE mhd_mod
  USE validate_mod
  USE stagger_mod
  USE parameters_mod
  USE eos_mod
  implicit none
  private
  type, public, extends(mhd_t):: solver_t
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: update
    procedure:: dealloc
    procedure:: dnload
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
    procedure:: gas_pressure_ngz
    procedure:: gas_pressure_sub
    procedure:: gas_temperature_sub
    procedure:: gas_temperature_ngz
    !procedure:: cell_gas_temperature
    procedure:: velocity_magnitude => void
    procedure:: magnetic_field_magnitude => void
    procedure:: grav_potential => void
    procedure:: apply_heating
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
!> Organize calls to the extras and the hydro solver
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
!> Organize calls to the extras and the hydro solver
!===============================================================================
SUBROUTINE init (self)
  class(solver_t) :: self
  !.............................................................................
  call trace%begin ('solver_t%init')
  call self%pre_init
  call self%mhd_t%init
  if (io%do_generic) then
    call self%from_generic()
    call self%deallocate_generic()
  end if
  call validate%init
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Call a hiearchy of %dealloc procedures, which should deallocate everything
!> allocated in the task
!===============================================================================
SUBROUTINE dealloc (self)
  class(solver_t) :: self
  !.............................................................................
  call trace%begin ('solver_t%dealloc')
  call self%mhd_t%dealloc
  call trace%end ()
END SUBROUTINE dealloc

!===============================================================================
!> Organize calls to the extras and the MHD solver.
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
  !-----------------------------------------------------------------------------
  ! Pre- and post-update procedures: as a precaution against missing calls,
  ! extra calls to self%gpatch_t%pre_update and post_update are added here
  !-----------------------------------------------------------------------------
  call self%gpatch_t%pre_update ()
  call self%mhd_t%pre_update
  call self%mhd_t%update
  call self%gpatch_t%post_update ()
  call self%mhd_t%post_update
  !-----------------------------------------------------------------------------
  call trace%end (itimer)
END SUBROUTINE update

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
!> Organize calls to the extras and the hydro solver
!===============================================================================
SUBROUTINE dnload (self, only)
  class(solver_t) :: self
  integer, optional:: only
  !.............................................................................
  call trace%begin ('solver_t%dnload')
  call self%mhd_t%dnload (only)
  call trace%end ()
END SUBROUTINE dnload

!=======================================================================
!> Get velocities from momenta
!=======================================================================
SUBROUTINE p2U(self, U, it)
  class(solver_t) :: self
  real, dimension(:,:,:,:), pointer:: U
  integer:: it
  !.....................................................................
  real, dimension(:,:,:,:), pointer:: dd
  integer:: i
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,           it,1), &
            p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1))
  do i=1,3
    U(:,:,:,i) = p(:,:,:,i)/d
  end do
  end associate
END SUBROUTINE p2U

!=======================================================================
!> Put velocities to momenta
!=======================================================================
SUBROUTINE U2p(self, U, it)
  class(solver_t) :: self
  real, dimension(:,:,:,:), pointer:: U
  integer:: it
  !.....................................................................
  real, dimension(:,:,:,:), pointer:: dd
  integer:: i
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d           ,it,1), &
            p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1))
  do i=1,3
    p(:,:,:,i) = U(:,:,:,i)*d
  end do
  end associate
END SUBROUTINE U2p

!=======================================================================
!> Get thermal energy per unit mass from entropy
!=======================================================================
SUBROUTINE e2E_th(self, E_th, it, lb, ub)
  class(solver_t) :: self
  integer:: it, lb(3), ub(3)
  optional:: lb, ub
  !.....................................................................
  real, dimension(:,:,:), pointer:: E_th
  real:: g1
  !---------------------------------------------------------------------
  if (present(lb)) then
    associate( &
      d => self%mem(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),self%idx%d,it,1), &
      e => self%mem(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),self%idx%e,it,1))
    E_th = e/d
    end associate
  else
    associate(d => self%mem(:,:,:,self%idx%d,it,1), &
              e => self%mem(:,:,:,self%idx%e,it,1))
    E_th = e/d
    end associate
  end if
END SUBROUTINE e2E_th

!=======================================================================
!> Put thermal energy to entropy per unit volume
!=======================================================================
SUBROUTINE E_th2e(self, E_th, it, lb, ub)
  class(solver_t) :: self
  integer:: it, lb(3), ub(3)
  optional:: lb, ub
  !.....................................................................
  real, dimension(:,:,:), pointer:: E_th
  real:: g1
  !---------------------------------------------------------------------
  if (present(lb)) then
    associate( &
      d => self%mem(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),self%idx%d,it,1), &
      e => self%mem(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),self%idx%e,it,1))
    e = d*E_th
    end associate
  else
    associate(d => self%mem(:,:,:,self%idx%d,it,1), &
              e => self%mem(:,:,:,self%idx%e,it,1))
    e = d*E_th
    end associate
  end if
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
FUNCTION up(f,i) RESULT (g)
  real, dimension(:,:,:), intent(in):: f
  real, dimension(size(f,1),size(f,2),size(f,3)):: g
  integer:: i
  !.............................................................................
  g = 0.5*(cshift(f,1,i)+f)
END FUNCTION

!===============================================================================
!> Subroutine interface to solver-independent gas pressure function
!===============================================================================
SUBROUTINE gas_pressure_sub (self, pg)
  class(solver_t):: self
  real, dimension(:,:,:), pointer, contiguous:: pg
  real, dimension(:,:,:), pointer, contiguous:: d, px, py, pz, e, eth, s
  real, dimension(:,:,:), pointer, contiguous:: bx, by, bz
  integer:: m(3)
  !-----------------------------------------------------------------------------
  d  => self%mem(:,:,:,self%idx%d ,self%it,1)
  if (isothermal) then
    pg = d*self%csound**2
  else
    select case(trim(self%etype))
    case('thermal')
      e  => self%mem(:,:,:,self%idx%e ,self%it,1)
      call eos%lookup (d=d, e=e, pg=pg)
    case('total')
      m = shape(d)
      e  => self%mem(:,:,:,self%idx%e ,self%it,1)
      px => self%mem(:,:,:,self%idx%px,self%it,1)
      py => self%mem(:,:,:,self%idx%py,self%it,1)
      pz => self%mem(:,:,:,self%idx%pz,self%it,1)
      allocate (eth(m(1),m(2),m(3)))
      eth = e-0.5*(px**2+py**2+pz**2)/d
      if (self%idx%bx>0) then
        bx  => self%mem(:,:,:,self%idx%bx ,self%it,1)
        by  => self%mem(:,:,:,self%idx%by ,self%it,1)
        bz  => self%mem(:,:,:,self%idx%bz ,self%it,1)
        eth = eth-0.5*(up(bx,1)**2+up(by,2)**2+up(bz,3)**2)
      end if
      call eos%lookup (d=d, e=eth, pg=pg)
      deallocate (eth)
    case('entropy')
      s  => self%mem(:,:,:,self%idx%s ,self%it,1)
      call eos%lookup (d=d, s=s, pg=pg)
    case default
      call io%abort("unknown etype in solver_mod::gas_pressure_sub")
    end select
  end if
END SUBROUTINE gas_pressure_sub

!===============================================================================
FUNCTION gas_pressure_ngz (self) RESULT (pg)
  class(solver_t):: self
  real, dimension(:,:,:), pointer:: pg
  real, dimension(:,:,:), pointer, contiguous:: d, px, py, pz, e, eth, s
  real, dimension(:,:,:), pointer, contiguous:: bx, by, bz
  integer :: l(3), u(3), m(3)
  !-----------------------------------------------------------------------------
  l = self%li; u = self%ui
  d  => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%d ,self%it,1)
  if (eos%gamma==1.0) then
    pg = d * self%csound**2
  else
    select case(trim(self%etype))
    case('thermal')
      e  => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e ,self%it,1)
      call eos%lookup (d=d, e=e, pg=pg)
    case('total')
      m = shape(d)
      e  => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e ,self%it,1)
      px => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%px,self%it,1)
      py => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%py,self%it,1)
      pz => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%pz,self%it,1)
      allocate (eth(m(1),m(2),m(3)))
      eth = e-0.5*(px**2+py**2+pz**2)/d
      if (self%idx%bx>0) then
        bx => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%bx,self%it,1)
        by => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%by,self%it,1)
        bz => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%bz,self%it,1)
        eth = eth-0.5*(up(bx,1)**2+up(by,2)**2+up(bz,3)**2)
      end if
      call eos%lookup (d=d, e=eth, pg=pg)
      deallocate (eth)
    case('entropy')
      s  => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%s ,self%it,1)
      call eos%lookup (d=d, s=s, pg=pg)
    case default
      call io%abort("unknown etype in solver_mod::gas_pressure_ngz")
    end select
  end if
END FUNCTION gas_pressure_ngz

!===============================================================================
!> Subroutine interface to solver-independent gas temperature function
!===============================================================================
SUBROUTINE gas_temperature_sub (self, tt)
  class(solver_t):: self
  real, dimension(:,:,:), pointer, contiguous:: tt
  real, dimension(:,:,:), pointer, contiguous:: d, px, py, pz, e, eth, s
  real, dimension(:,:,:), pointer, contiguous:: bx, by, bz
  integer:: m(3)
  !-----------------------------------------------------------------------------
  d  => self%mem(:,:,:,self%idx%d ,self%it,1)
  if (isothermal) then
    tt = d*self%csound**2
  else
    select case(trim(self%etype))
    case('thermal')
      e  => self%mem(:,:,:,self%idx%e ,self%it,1)
      call eos%lookup (d=d, e=e, tt=tt)
    case('total')
      m = shape(d)
      e  => self%mem(:,:,:,self%idx%e ,self%it,1)
      px => self%mem(:,:,:,self%idx%px,self%it,1)
      py => self%mem(:,:,:,self%idx%py,self%it,1)
      pz => self%mem(:,:,:,self%idx%pz,self%it,1)
      allocate (eth(m(1),m(2),m(3)))
      eth = e-0.5*(px**2+py**2+pz**2)/d
      if (self%idx%bx>0) then
        bx  => self%mem(:,:,:,self%idx%bx ,self%it,1)
        by  => self%mem(:,:,:,self%idx%by ,self%it,1)
        bz  => self%mem(:,:,:,self%idx%bz ,self%it,1)
        eth = eth-0.5*(up(bx,1)**2+up(by,2)**2+up(bz,3)**2)
      end if
      call eos%lookup (d=d, e=eth, tt=tt)
      deallocate (eth)
    case('entropy')
      s  => self%mem(:,:,:,self%idx%s ,self%it,1)
      call eos%lookup (d=d, s=s, tt=tt)
    case default
      call io%abort("unknown etype in solver_mod::gas_temperature_sub")
    end select
  end if
END SUBROUTINE gas_temperature_sub

!===============================================================================
FUNCTION gas_temperature_ngz (self) RESULT (tt)
  class(solver_t):: self
  real, dimension(:,:,:), pointer:: tt
  real, dimension(:,:,:), pointer, contiguous:: d, px, py, pz, e, eth, s
  real, dimension(:,:,:), pointer, contiguous:: bx, by, bz
  integer :: l(3), u(3), m(3)
  !-----------------------------------------------------------------------------
  l = self%li; u = self%ui
  d  => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%d ,self%it,1)
  if (eos%gamma==1.0) then
    tt = d * self%csound**2
  else
    select case(trim(self%etype))
    case('thermal')
      e  => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e ,self%it,1)
      call eos%lookup (d=d, e=e, tt=tt)
    case('total')
      m = shape(d)
      e  => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e ,self%it,1)
      px => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%px,self%it,1)
      py => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%py,self%it,1)
      pz => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%pz,self%it,1)
      allocate (eth(m(1),m(2),m(3)))
      eth = e-0.5*(px**2+py**2+pz**2)/d
      if (self%idx%bx>0) then
        bx => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%bx,self%it,1)
        by => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%by,self%it,1)
        bz => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%bz,self%it,1)
        eth = eth-0.5*(up(bx,1)**2+up(by,2)**2+up(bz,3)**2)
      end if
      call eos%lookup (d=d, e=eth, tt=tt)
      deallocate (eth)
    case('entropy')
      s  => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%s ,self%it,1)
      call eos%lookup (d=d, s=s, tt=tt)
    case default
      call io%abort("unknown etype in solver_mod::gas_temeprature_ngz")
    end select
  end if
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
  call self%gas_pressure_sub(v)
  v = log(v)
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
SUBROUTINE magnetic_field_magnitude (self, v)
  class (solver_t):: self
  real, dimension(:,:,:), pointer:: v
  real, dimension(:,:,:,:), pointer:: b
  integer:: i
  !-----------------------------------------------------------------------------
  associate (bx => self%mem(:,:,:,self%idx%bx,self%it,1), &
             by => self%mem(:,:,:,self%idx%by,self%it,1), &
             bz => self%mem(:,:,:,self%idx%bz,self%it,1))
  v = up(bx,1)**2 + up(by,2)**2 + up(bz,3)**2
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
  w = max(- stagger%ddx(self%ds,vx) &
          - stagger%ddy(self%ds,vy) &
          - stagger%ddz(self%ds,vz), 0.0)
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
  w = sqrt((stagger%ddy(self%ds,vz)-stagger%ddz(self%ds,vy))**2 &
         + (stagger%ddz(self%ds,vx)-stagger%ddx(self%ds,vz))**2 &
         + (stagger%ddx(self%ds,vy)-stagger%ddy(self%ds,vx))**2)
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
  real, pointer, contiguous, dimension(:,:,:):: d, e, px, py, pz, bx, by, bz
  !-----------------------------------------------------------------------------
  call trace%begin ('solver_t%to_thermal')
  select case(trim(self%etype))
  case('thermal')
    e  => self%mem(:,:,:,self%idx%e ,self%it,1)
    eth = e
  case('total')
    d  => self%mem(:,:,:,self%idx%d ,it,1)
    e  => self%mem(:,:,:,self%idx%e ,it,1)
    px => self%mem(:,:,:,self%idx%px,it,1)
    py => self%mem(:,:,:,self%idx%py,it,1)
    pz => self%mem(:,:,:,self%idx%pz,it,1)
    eth = e - 0.5*(px**2+py**2+pz**2)/d
    if (self%idx%bx>0) then
      bx => self%mem(:,:,:,self%idx%bx,it,1)
      by => self%mem(:,:,:,self%idx%by,it,1)
      bz => self%mem(:,:,:,self%idx%bz,it,1)
      do i3=1,self%gn(3)-1
      do i2=1,self%gn(2)-1
      do i1=1,self%gn(1)-1
        eth(i1,i2,i3) = &
        eth(i1,i2,i3) - 0.125*((bx(i1,i2,i3)+bx(i1+1,i2,i3))**2 + &
                                    (by(i1,i2,i3)+by(i1,i2+1,i3))**2 + &
                                    (bz(i1,i2,i3)+bz(i1,i2,i3+1))**2)
      end do
      end do
      end do
    end if
  case default
    call io%abort("unknown etype in solver_mod::to_thermal")
  end select
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
      em = em * dm
    case('entropy')
      Tg => self%generic(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%s)
      ss => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%s,self%it,1)
      call eos%inverse_lookup (d=dm, v=Tg, name='T', ss=em)
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
  real(var_k), pointer, dimension(:,:,:):: dm, db, em, Tg, s
  integer:: iv, l(3), u(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('solver_t%to_generic')
  if (associated(self%generic)) then
    l = merge (self%mesh%lb, self%mesh%li, io%guard_zones)
    u = merge (self%mesh%ub, self%mesh%ui, io%guard_zones)
    ! Density
    db => self%generic(:,:,:,self%idx%d)
    dm => self%mem(:,:,:,self%idx%d,self%it,1)
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
