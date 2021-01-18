!===============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!===============================================================================
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
  real(8), parameter:: pi2=8d0*atan(1d0)
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
    procedure, nopass:: cleanup
  end type
  character(len=32), save:: type='turbulence'
  real, save, pointer, dimension(:,:,:,:):: buffer => null()
  real, save:: csound=1.
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self, solver, gamma)
  class(initial_t):: self
  character(len=64):: solver
  real(8), optional:: gamma
  logical, save:: first_time=.true.
  namelist /IC_params/ type
  !----------------------------------------------------------------------------
  call trace_begin('initial_t%init')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, IC_params)
    if (io%master) write (*, IC_params)
  end if
  call io%lock%unset
  self%solver = trim(solver)
  if (present(gamma)) then
    self%gamma = gamma
  end if
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE condition (self, m, ff, idx)
  class(initial_t) :: self
  type(index_t):: idx
  class(mesh_t)    :: m(3)
  real, dimension(:,:,:,:) :: ff
  !.............................................................................
  real(8):: g1
  real, dimension(:,:,:), pointer:: d, e, s
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t%condtion')
  select case (type)
  case ('single_solenoidal')
    call single_solenoidal (self, m, ff, idx)
  case ('single_compressive')
    call single_compressive (self, m, ff, idx)
  case ('compressive_shear')
    call compressive_shear (self, m, ff, idx)
  case ('single_advection')
    call single_advection (self, m, ff, idx)
  case ('turbulence')
    call turbulence (self, m, ff, idx)
  case ('exa256')
    call exa256 (self, m, ff, idx)
  case ('raw')
    call raw (self, m, ff, idx)
  case default
    call io%abort ('unknown IC')
  end select
  call convert (self, m, ff, idx)
  call trace_end ()
END SUBROUTINE condition

!===============================================================================
!> Convert to entropy
!===============================================================================
SUBROUTINE convert (self, m, ff, idx)
  class(initial_t) :: self
  type(index_t):: idx
  class(mesh_t)    :: m(3)
  real, dimension(:,:,:,:),target:: ff
  !.............................................................................
  real(8):: g1
  real, dimension(:,:,:), pointer:: d, e, s
  !-----------------------------------------------------------------------------
  s => ff(:,:,:,idx%s)
  if (self%gamma == 1d0) then
    s = 0.0
  else
    d => ff(:,:,:,idx%d)
    e => ff(:,:,:,idx%e)
    g1 = self%gamma-1d0
    s = d*log(e*g1/d**self%gamma)/g1
  end if
  call trace_end ()
END SUBROUTINE convert

!===============================================================================
SUBROUTINE cleanup
  if (associated(buffer)) deallocate(buffer)
END SUBROUTINE cleanup

!===============================================================================
!> Initialize a single 3-D wavenumber, with random amplitude and phase
!===============================================================================
SUBROUTINE single_solenoidal (self, m, ff, idx)
  class(initial_t) :: self
  type(index_t):: idx
  class(mesh_t)    :: m(3)
  real, dimension(:,:,:,:) :: ff
  !.............................................................................
  type(random_t) :: random
  integer        :: ix, iy, iz
  real(8)        :: kx, ky, kz, gamma, gamma_factor
  integer, save  :: seed=-33
  logical, save  :: mhd=.false.
  real, save     :: a0(3)=0.1, k(3)=1.0, d0=1.0, s0=0.0, b0(3)=0.0, u0(3)=0.0, phi(3)=0.0
  real(8), save  :: t_turn=0.1
  logical, save  :: first_time=.true.
  namelist /IC_solenoidal_params/ seed, k, a0, u0, b0, d0, s0, mhd, gamma, t_turn
  integer, save  :: itimer=0
  real           :: d, p
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:single_solenoidal', itimer=itimer)
  gamma = self%gamma
  if (first_time) then
    call io%lock%set
    if (first_time) then
      rewind (io%input)
      read (io%input, IC_solenoidal_params)
      mhd = mhd .and. any(b0 /= 0.0)
      if (io%master) write (*, IC_solenoidal_params)
      first_time = .false.
    end if
    call io%lock%unset
  end if
  self%gamma = gamma
  if (gamma==1.00) then
    gamma_factor = 1.0
  else
    gamma_factor = 1.0/(gamma*(gamma-1))
  end if
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    kz   =  pi2*k(3)*(r3(iz)+m(3)%p-0.5d0*m(3)%s)/m(3)%b + phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky   =  pi2*k(2)*(r2(iy)+m(2)%p-0.5d0*m(2)%s)/m(2)%b + phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*k(1)*(r1(ix)+m(1)%p-0.5d0*m(1)%s)/m(1)%b + phi(1)
        ff(ix,iy,iz,idx%d ) = d0
        ff(ix,iy,iz,idx%px) = a0(1)*(sin(ky)+sin(kz))
        ff(ix,iy,iz,idx%py) = a0(2)*(sin(kz)+sin(kx))
        ff(ix,iy,iz,idx%pz) = a0(3)*(sin(kx)+sin(ky))
        ff(ix,iy,iz,idx%e ) = d0*csound**2*gamma_factor
        if (self%mhd) then
          ff(ix,iy,iz,idx%bx) = b0(1)*(sin(kx)+sin(ky))
          ff(ix,iy,iz,idx%by) = b0(2)*(sin(kx)+sin(ky))
          ff(ix,iy,iz,idx%bz) = b0(3)*(sin(kx)+sin(ky))
        end if
      end do
    end do
  end do
  end associate
  call trace_end (itimer)
END SUBROUTINE single_solenoidal

!===============================================================================
!> Initialize a single 3-D wavenumber, with non-random amplitude and zero phase
!===============================================================================
SUBROUTINE single_compressive (self, m, ff, idx)
  class(initial_t) :: self
  type(index_t):: idx
  class(mesh_t)    :: m(3)
  real, dimension(:,:,:,:) :: ff
  !.............................................................................
  type(random_t) :: random
  integer        :: ix, iy, iz
  real(8)        :: kx, ky, kz, gamma, gamma_factor
  integer, save  :: seed=-33
  logical, save  :: mhd=.false.
  real, save     :: a0(3)=0.1, k(3)=1.0, d0=1.0, s0=0.0, b0(3)=0.0, u0(3)=0.0, phi(3)=0.0
  real(8), save  :: t_turn=0.1
  logical, save  :: first_time=.true.
  namelist /IC_compressive_params/ seed, k, a0, u0, b0, d0, s0, mhd, gamma, t_turn
  integer, save  :: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:single_compressive', itimer=itimer)
  gamma = self%gamma
  if (first_time) then
    call io%lock%set
    if (first_time) then
      rewind (io%input)
      read (io%input, IC_compressive_params)
      mhd = mhd .and. any(b0 /= 0.0)
      if (io%master) write (*, IC_compressive_params)
      first_time = .false.
    end if
    call io%lock%unset
  end if
  self%gamma = gamma
  if (gamma==1.00) then
    gamma_factor = 1.0
  else
    gamma_factor = 1.0/(gamma*(gamma-1))
  end if
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    kz   =  pi2*k(3)*(r3(iz)+m(3)%p-0.5d0*m(3)%s)/m(3)%b + phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky   =  pi2*k(2)*(r2(iy)+m(2)%p-0.5d0*m(2)%s)/m(2)%b + phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*k(1)*(r1(ix)+m(1)%p-0.5d0*m(1)%s)/m(1)%b + phi(1)
        ff(ix,iy,iz,idx%d ) = d0
        ff(ix,iy,iz,idx%px) = u0(1)*sin(kx)
        ff(ix,iy,iz,idx%py) = u0(2)*sin(ky)
        ff(ix,iy,iz,idx%pz) = u0(3)*sin(kz)
        ff(ix,iy,iz,idx%e ) = d0*csound**2*gamma_factor
        if (self%mhd) then
          ff(ix,iy,iz,idx%bx) = b0(1)
          ff(ix,iy,iz,idx%by) = b0(2)
          ff(ix,iy,iz,idx%bz) = b0(3)
        end if
      end do
    end do
  end do
  end associate
  call trace_end (itimer)
END SUBROUTINE single_compressive

!===============================================================================
!> Initialize a single 3-D wavenumber, with non-random amplitude and zero phase
!===============================================================================
SUBROUTINE compressive_shear (self, m, ff, idx)
  class(initial_t) :: self
  type(index_t):: idx
  class(mesh_t)    :: m(3)
  real, dimension(:,:,:,:) :: ff
  !.............................................................................
  type(random_t) :: random
  integer        :: ix, iy, iz
  real(8)        :: kx, ky, kz, gamma, gamma_factor
  integer, save  :: seed=-33, verbose=0
  logical, save  :: mhd=.false.
  real, save     :: a0(3)=0.1, k(3)=1.0, d0=1.0, s0=0.0, b0(3)=0.0, u0(3)=0.0, phi(3)=0.0
  real(8), save  :: t_turn=0.1
  logical, save  :: first_time=.true.
  namelist /IC_compressive_params/ seed, k, a0, u0, b0, d0, s0, mhd, gamma, &
                                   t_turn, verbose
  integer, save  :: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:compressive_shear', itimer=itimer)
  gamma = self%gamma
  if (first_time) then
    call io%lock%set
    if (first_time) then
      rewind (io%input)
      read (io%input, IC_compressive_params)
      mhd = mhd .and. any(b0 /= 0.0)
      if (io%master) write (*, IC_compressive_params)
      write(stdout,'(a,1p,10g12.4)') 'compressive_shear: pi2, box, stag =', &
        pi2, m(1)%b, m(1)%h
      first_time = .false.
    end if
    call io%lock%unset
  end if
  self%gamma = gamma
  if (gamma==1.00) then
    gamma_factor = 1.0
  else
    gamma_factor = 1.0/(gamma*(gamma-1))
  end if
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*k(1)*(r1(ix)+m(1)%p+m(1)%h(idx%px)*m(1)%d)/m(1)%b + phi(1)
        ff(ix,iy,iz,idx%d ) = d0
        ff(ix,iy,iz,idx%px) = u0(1)*cos(kx)
        kx   =  pi2*k(1)*(r1(ix)+m(1)%p)/m(1)%b + phi(1)
        ff(ix,iy,iz,idx%py) = u0(2)*cos(kx)
        ff(ix,iy,iz,idx%pz) = u0(3)*cos(kx)
        ff(ix,iy,iz,idx%e ) = d0*csound**2*gamma_factor
        if (self%mhd) then
          ff(ix,iy,iz,idx%bx) = b0(1)
          ff(ix,iy,iz,idx%by) = b0(2)
          ff(ix,iy,iz,idx%bz) = b0(3)
        end if
      end do
    end do
  end do
  end associate
  call trace_end (itimer)
END SUBROUTINE compressive_shear

!===============================================================================
!> Initialize a single 3-D wavenumber, with random amplitude and phase
!===============================================================================
SUBROUTINE single_advection (self, m, ff, idx)
  class(initial_t) :: self
  type(index_t):: idx
  class(mesh_t)    :: m(3)
  real, dimension(:,:,:,:) :: ff
  !.............................................................................
  type(random_t) :: random
  integer        :: ix, iy, iz
  real(8)        :: kx, ky, kz, gamma, gamma_factor
  integer, save  :: seed=-33
  logical, save  :: mhd=.false.
  logical, save  :: first_time=.true.
  real, save     :: a0(3)=0.1, k(3)=1.0, d0=1.0, s0=0.0, b0(3)=0.0, u0(3)=0.0, phi(3)=0.0
  namelist /IC_advection_params/ seed, k, a0, u0, b0, d0, s0, mhd, gamma
  integer, save  :: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:single_advection', itimer=itimer)
  gamma = self%gamma
  if (first_time) then
    call io%lock%set
    if (first_time) then
      rewind (io%input)
      read (io%input, IC_advection_params)
      mhd = mhd .and. any(b0 /= 0.0)
      if (io%master) write (*, IC_advection_params)
      first_time = .false.
    end if
    call io%lock%unset
  end if
  self%gamma = gamma
  if (gamma==1.00) then
    gamma_factor = 1.0
  else
    gamma_factor = 1.0/(gamma*(gamma-1))
  end if
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    kz   =  pi2*k(3)*(r3(iz)+m(3)%p-0.5d0*m(3)%s)/m(3)%b + phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky   =  pi2*k(2)*(r2(iy)+m(2)%p-0.5d0*m(2)%s)/m(2)%b + phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*k(1)*(r1(ix)+m(1)%p-0.5d0*m(1)%s)/m(1)%b + phi(1)
        ff(ix,iy,iz,idx%d ) = d0*(1.+a0(1)*sin(kx)+a0(2)*sin(ky)+a0(3)*sin(kz))
        ff(ix,iy,iz,idx%px) = d0*u0(1)
        ff(ix,iy,iz,idx%py) = d0*u0(2)
        ff(ix,iy,iz,idx%pz) = d0*u0(3)
        ff(ix,iy,iz,idx%e ) = d0*csound**2*gamma_factor
        if (self%mhd) then
          ff(ix,iy,iz,idx%bx) = b0(3)*sin(kz)
          ff(ix,iy,iz,idx%by) = b0(1)*sin(kx)
          ff(ix,iy,iz,idx%bz) = b0(2)*sin(ky)
        end if
      end do
    end do
  end do
  end associate
  call trace_end (itimer)
END SUBROUTINE single_advection

!===============================================================================
!> Initialize from a raw data file
!===============================================================================
SUBROUTINE exa256 (self, m, ff, idx)
  class(initial_t):: self
  type(index_t):: idx
  class(mesh_t):: m(3)
  real, dimension(:,:,:,:):: ff
  real, dimension(:), allocatable:: u, d
  real, save:: courant_factor=0.2
  !.............................................................................
  integer:: ix, iy, iz, jy, jz, rec, iv, li(3), ui(3)
  integer(8):: nx, ny, nz
  real(8):: umax, dt, gamma, gamma_factor
  logical, save:: first_time=.true.
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t::exa256', itimer=itimer)
  self%mhd = .false.
  ff(:,:,:,:) = 0.0
  allocate (u(m(1)%n))
  allocate (d(m(1)%n))
  if (first_time) then
    first_time = .false.
    open (io_unit%data, file='data/exa256/t0.02.raw', form='unformatted', &
      access='direct', recl=m(1)%n*4)
  end if
  self%time = 0.02
  ix=m(1)%p/m(1)%d+0.5-m(1)%n/2
  iy=m(2)%p/m(2)%d+0.5-m(2)%n/2
  iz=m(3)%p/m(3)%d+0.5-m(3)%n/2
  nx = m(1)%b/m(1)%d + 0.5
  ny = m(2)%b/m(2)%d + 0.5
  nz = m(3)%b/m(3)%d + 0.5
  li = m%li
  ui = m%ui
  umax = 0.0
  gamma = self%gamma
  if (gamma==1.00) then
    gamma_factor = 1.0
  else
    gamma_factor = 1.0/(gamma*(gamma-1))
  end if
  do jz=li(3),ui(3)
  do jy=li(2),ui(2)
    rec = 1   + (ix           )      /m(1)%n
    rec = rec + (iy+jy-m(2)%li)*nx   /m(1)%n
    rec = rec + (iz+jz-m(3)%li)*nx*ny/m(1)%n
    read (io_unit%data,rec=rec) d; ff(m(1)%li:m(1)%ui,jy,jz,idx%d ) = d
    rec = rec + nx*ny*nz/m(1)%n
    read (io_unit%data,rec=rec) u; ff(m(1)%li:m(1)%ui,jy,jz,idx%px) = u
    rec = rec + nx*ny*nz/m(1)%n
    read (io_unit%data,rec=rec) u; ff(m(1)%li:m(1)%ui,jy,jz,idx%py) = u
    rec = rec + nx*ny*nz/m(1)%n
    read (io_unit%data,rec=rec) u; ff(m(1)%li:m(1)%ui,jy,jz,idx%pz) = u
    umax=max(umax,3.+maxval(abs(ff(li(1):ui(1),jy,jz,idx%px)) &
                           +abs(ff(li(1):ui(1),jy,jz,idx%py)) &
                           +abs(ff(li(1):ui(1),jy,jz,idx%pz))))
  end do
  end do
  ff(:,:,:,idx%e) = ff(:,:,:,idx%d)*gamma_factor
  do iv=idx%px,idx%pz
    ff(:,:,:,iv) = ff(:,:,:,idx%d)*ff(:,:,:,iv)
  end do
  dt = courant_factor*m(1)%d/umax
  if (io%verbose > 2) write (io_unit%log,*) 'initial dt:',dt,courant_factor,self%gamma,umax
  if (io%verbose > 2) &
    write (io_unit%log,'(a,3i4,1p,2g12.3)') 'ix,iy,iz,dmin,dmax =', ix, iy, iz, &
      minval(ff(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui,1)), &
      maxval(ff(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui,1))
  call trace_end (itimer)
END SUBROUTINE exa256

!===============================================================================
!> Initialize from a raw data file, assumed to contain isothermal data, in the
!> order rho, Ux, Uy, Uz
!===============================================================================
SUBROUTINE raw (self, m, ff, idx)
  class(initial_t):: self
  type(index_t):: idx
  class(mesh_t):: m(3)
  real, dimension(:,:,:,:):: ff
  !.............................................................................
  integer, parameter:: mv=4
  integer:: n(3), i(3), iv, jv(mv)
  logical, save:: first_time=.true.
  integer, save         :: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t::raw', itimer=itimer)
  !$omp critical (raw_cr)
  self%mhd = .false.
  self%time = 0.02
  !-----------------------------------------------------------------------------
  ! At first call, allocate a buffer and read the entire file -- this is OK since
  ! we only use one thread per rank to initialize, computing also how many calls
  ! will be made, so the buffer can be deallocated when all patches are done
  !-----------------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    n = nint(m%b/m%d)
    allocate (buffer(n(1),n(2),n(3),mv))
    open (io_unit%data, file='data/exa256/t0.02.raw', form='unformatted', &
      access='direct', recl=4*product(n))
    do iv=1,mv
      read (io_unit%data, rec=iv) buffer(:,:,:,iv)
    end do
    close (io_unit%data)
  end if
  i = 1 + nint(m%p/m%d) - m%n/2                         ! index offsets
  jv(1) = idx%d                                           ! variable index translation
  jv(2) = idx%px
  jv(3) = idx%py
  jv(4) = idx%pz
  ff = 0.0
  do iv=1,mv
    ff(m(1)%li:m(1)%ui, &
       m(2)%li:m(2)%ui, &
       m(3)%li:m(3)%ui,jv(iv)) = buffer(i(1):i(1)+m(1)%n-1, &
                                        i(2):i(2)+m(2)%n-1, &
                                        i(3):i(3)+m(3)%n-1,iv)
  end do
  do iv=idx%px,idx%pz
    ff(:,:,:,iv) = ff(:,:,:,idx%d)*ff(:,:,:,iv)
  end do
  if (io%verbose > 2) &
    write (io_unit%log,'(a,3i4,1p,2g12.3)') 'ix,iy,iz,dmin,dmax =', i, &
      minval(ff(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui,1)), &
      maxval(ff(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui,1))
  !$omp end critical (raw_cr)
  call trace_end (itimer)
END SUBROUTINE raw

!===============================================================================
!> Initialize to zero, and let the forces drive from there
!===============================================================================
SUBROUTINE turbulence (self, m, ff, idx)
  class(initial_t):: self
  type(index_t):: idx
  class(mesh_t):: m(3)
  real, dimension(:,:,:,:):: ff
  real :: d, p
  !.............................................................................
  integer, save         :: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t::raw', itimer=itimer)
  d = 1.0
  ff(:,:,:,idx%d        ) = d
  if (self%gamma==1.0) then
    ff(:,:,:,idx%e        ) = d*csound**2
  else
    ff(:,:,:,idx%e        ) = d*csound**2/(self%gamma-1.)
  end if
  ff(:,:,:,idx%px:idx%pz) = 0.0
  call trace_end (itimer)
END SUBROUTINE turbulence

END MODULE
