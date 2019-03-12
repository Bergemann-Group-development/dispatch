!===============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!===============================================================================
MODULE initial_mod
  !USE const
  USE io_mod
  USE mpi_mod
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
    procedure(single_solenoidal), pointer:: condition=>null()
  contains
    procedure:: init
    procedure, nopass:: cleanup
  end type
  real, save, pointer, dimension(:,:,:,:):: buffer => null()
  real, save:: csound=1.
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self, solver, gamma)
  class(initial_t):: self
  character(len=64):: solver
  real, optional:: gamma
  character(len=32), save:: type='turbulence'
  logical, save:: first_time=.true.
  namelist /initial_params/ type
  !----------------------------------------------------------------------------
  call trace_begin('initial_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, initial_params)
    if (io%master) write (*, initial_params)
  end if
  !$omp end critical (input_cr)
  self%solver = trim(solver)
  if (present(gamma)) then
    self%gamma = gamma
  end if
  select case (type)
  case ('single_solenoidal')
    self%condition => single_solenoidal
  case ('single_compressive')
    self%condition => single_compressive
  case ('single_advection')
    self%condition => single_advection
  case ('turbulence')
    self%condition => turbulence
  case ('exa256')
    self%condition => exa256
  case ('raw')
    self%condition => raw
  case default
    print*,'UNKNOWN INITIAL CONDITION TYPE :'//trim(type)//':'
    stop
  end select
  call trace_end
END SUBROUTINE init

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
  real(8)        :: kx, ky, kz, gamma
  integer, save  :: seed=-33
  logical, save  :: mhd=.false.
  real, save     :: a0(3)=0.1, k(3)=1.0, d0=1.0, s0=0.0, b0(3)=0.0, u0(3)=0.0, phi(3)=0.0
  real(8), save  :: t_turn=0.1
  logical, save  :: first_time=.true.
  namelist /initial_solenoidal_params/ seed, k, a0, u0, b0, d0, s0, mhd, gamma, t_turn
  integer, save  :: itimer=0
  real           :: d, p
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t%single_solenoidal', itimer=itimer)
  gamma = self%gamma
  if (first_time) then
    !$omp critical (input_cr)
    if (first_time) then
      first_time = .false.
      rewind (io%input)
      read (io%input, initial_solenoidal_params)
      mhd = mhd .and. any(b0 /= 0.0)
      if (io%master) write (*, initial_solenoidal_params)
    end if
    !$omp end critical (input_cr)
  end if
  self%gamma = gamma
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
        ff(ix,iy,iz,idx%e ) = d0*csound**2*(gamma*(gamma-1.0))
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
  real(8)        :: kx, ky, kz, gamma
  integer, save  :: seed=-33
  logical, save  :: mhd=.false.
  real, save     :: a0(3)=0.1, k(3)=1.0, d0=1.0, s0=0.0, b0(3)=0.0, u0(3)=0.0, phi(3)=0.0
  real(8), save  :: t_turn=0.1
  logical, save  :: first_time=.true.
  namelist /initial_compressive_params/ seed, k, a0, u0, b0, d0, s0, mhd, gamma, t_turn
  integer, save  :: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t%single_compressive', itimer=itimer)
  gamma = self%gamma
  if (first_time) then
    !$omp critical (input_cr)
    if (first_time) then
      first_time = .false.
      rewind (io%input)
      read (io%input, initial_compressive_params)
      mhd = mhd .and. any(b0 /= 0.0)
      if (io%master) write (*, initial_compressive_params)
    end if
    !$omp end critical (input_cr)
  end if
  self%gamma = gamma
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
        if (self%gamma==1.0) then
          ff(ix,iy,iz,idx%e ) = d0*csound**2
        else
          ff(ix,iy,iz,idx%e ) = d0*csound**2/(gamma*(gamma-1))
        end if
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
  real(8)        :: kx, ky, kz, gamma
  integer, save  :: seed=-33
  logical, save  :: mhd=.false.
  logical, save  :: first_time=.true.
  real, save     :: a0(3)=0.1, k(3)=1.0, d0=1.0, s0=0.0, b0(3)=0.0, u0(3)=0.0, phi(3)=0.0
  namelist /initial_advection_params/ seed, k, a0, u0, b0, d0, s0, mhd, gamma
  integer, save  :: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t%single_advection', itimer=itimer)
  gamma = self%gamma
  if (first_time) then
    !$omp critical (input_cr)
    if (first_time) then
      first_time = .false.
      rewind (io%input)
      read (io%input, initial_advection_params)
      mhd = mhd .and. any(b0 /= 0.0)
      if (io%master) write (*, initial_advection_params)
    end if
    !$omp end critical (input_cr)
  end if
  self%gamma = gamma
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
        if (self%gamma==1.0) then
          ff(ix,iy,iz,idx%e ) = d0*csound**2
        else
          ff(ix,iy,iz,idx%e ) = d0*csound**2/(gamma*(gamma-1))
        end if
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
  real(8):: umax, dt
  logical, save:: first_time=.true.
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t%exa256', itimer=itimer)
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
  ff(:,:,:,5) = 0.5*ff(:,:,:,idx%d)*(ff(:,:,:,idx%px)**2+ &
                                     ff(:,:,:,idx%py)**2+ &
                                     ff(:,:,:,idx%pz)**2) &
              + ff(:,:,:,idx%d)/(self%gamma-1d0)
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
  real, save            :: b0(3)=[0.,0.,0.1]
  integer, save         :: itimer=0
  logical, save         :: first_time=.true.
  namelist /IC_turbulence_params/ b0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t%turbulence', itimer=itimer)
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, IC_turbulence_params)
    write (io%output, IC_turbulence_params)
  end if
  !$omp end critical (input_cr)
  if (any(b0 /= 0.0)) then
    self%mhd = .true.
  else
    self%mhd = .false.
  end if
  d = 1.0
  ff(:,:,:,idx%d        ) = d
  if (self%gamma==1.0) then
    ff(:,:,:,idx%e      ) = d*csound**2
  else
    ff(:,:,:,idx%e      ) = d*csound**2/(self%gamma-1.)
  end if
  ff(:,:,:,idx%px:idx%pz) = 0.0
  if (self%mhd .and. idx%bx > 0) then
    ff(:,:,:,idx%bx     ) = b0(1)
    ff(:,:,:,idx%by     ) = b0(2)
    ff(:,:,:,idx%bz     ) = b0(3)
  end if
  call trace_end (itimer)
END SUBROUTINE turbulence

END MODULE
