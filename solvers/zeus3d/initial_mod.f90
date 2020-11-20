!===============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE variable_mod
  USE mesh_mod
  USE random_mod
  USE index_mod
  implicit none
  private
  real(8), parameter:: pi2=8d0*atan(1d0)
  type, public:: initial_t
    integer:: seed
    logical:: mhd
    real(8):: time
    character(len=64):: solver
    procedure(single_solenoidal), pointer:: condition=>null()
  contains
    procedure:: init
    procedure, nopass:: cleanup
  end type
  real, save, pointer, dimension(:,:,:,:):: buffer => null()
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self, solver)
  class(initial_t):: self
  character(len=64):: solver
  character(len=32), save:: type='single_solenoidal'
  logical, save:: first_time=.true.
  namelist /IC_params/ type
  !----------------------------------------------------------------------------
  call trace_begin('initial_t%init')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, initial_params)
    if (io%master) write (*, initial_params)
  end if
  call io%lock%unset
  self%solver = trim(solver)
  select case (type)
  case ('single_solenoidal')
    self%condition => single_solenoidal
  case ('single_compressive')
    self%condition => single_compressive
  case ('single_advection')
    self%condition => single_advection
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
  class(initial_t):: self
  type(index_t):: idx
  class(mesh_t), pointer, dimension(:):: m
  real(kind=KindScalarVar), dimension(:,:,:,:) :: ff
  !.............................................................................
  type(random_t) :: random
  integer        :: ix, iy, iz
  real(8)        :: kx, ky, kz
  integer, save  :: seed=-33
  logical, save  :: mhd=.false.
  real, save     :: a0(3)=0.1, k(3)=1.0, d0=1.0, e0=1.0, b0(3)=0.0, u0(3)=0.0, phi(3)=0.0
  real(8), save  :: t_turn=0.1
  logical, save  :: first_time=.true.
  namelist /IC_solenoidal_params/ seed, k, a0, u0, b0, d0, e0, mhd, t_turn
  integer, save  :: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:single_solenoidal', itimer=itimer)
  if (first_time) then
    call io%lock%set
    if (first_time) then
      rewind (io%input)
      read (io%input, initial_solenoidal_params)
      mhd = mhd .and. any(b0 /= 0.0)
      if (io%master) write (*, initial_solenoidal_params)
      first_time = .false.
    end if
    call io%lock%unset
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
        ff(ix,iy,iz,idx%e) = e0
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
  class(initial_t):: self
  type(index_t):: idx
  class(mesh_t), pointer, dimension(:):: m
  real(kind=KindScalarVar), dimension(:,:,:,:) :: ff
  !.............................................................................
  type(random_t) :: random
  integer        :: ix, iy, iz
  real(8)        :: kx, ky, kz
  integer, save  :: seed=-33
  logical, save  :: mhd=.false.
  real, save     :: a0(3)=0.1, k(3)=1.0, d0=1.0, e0=1.0, b0(3)=0.0, u0(3)=0.0, phi(3)=0.0
  real(8), save  :: t_turn=0.1
  logical, save  :: first_time=.true.
  namelist /IC_compressive_params/ seed, k, a0, u0, b0, d0, e0, mhd, t_turn
  integer, save  :: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:single_compressive', itimer=itimer)
  if (first_time) then
    call io%lock%set
    if (first_time) then
      rewind (io%input)
      read (io%input, initial_compressive_params)
      mhd = mhd .and. any(b0 /= 0.0)
      if (io%master) write (*, initial_compressive_params)
      first_time = .false.
    end if
    call io%lock%unset
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
        ff(ix,iy,iz,idx%e) = e0
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
  class(initial_t):: self
  type(index_t):: idx
  class(mesh_t), pointer, dimension(:):: m
  real(kind=KindScalarVar), dimension(:,:,:,:) :: ff
  !.............................................................................
  type(random_t) :: random
  integer        :: ix, iy, iz
  real(8)        :: kx, ky, kz
  integer, save  :: seed=-33
  logical, save  :: mhd=.false.
  logical, save  :: first_time=.true.
  real, save     :: a0(3)=0.1, k(3)=1.0, d0=1.0, e0=1.0, b0(3)=0.0, u0(3)=0.0, phi(3)=0.0
  namelist /IC_advection_params/ seed, k, a0, u0, b0, d0, e0, mhd
  integer, save  :: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:single_advection', itimer=itimer)
  if (first_time) then
    call io%lock%set
    if (first_time) then
      rewind (io%input)
      read (io%input, initial_advection_params)
      mhd = mhd .and. any(b0 /= 0.0)
      if (io%master) write (*, initial_advection_params)
      first_time = .false.
    end if
    call io%lock%unset
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
        ff(ix,iy,iz,idx%e) = e0
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
  USE io_unit_mod
  class(initial_t):: self
  type(index_t):: idx
  class(mesh_t), pointer, dimension(:):: m
  real(kind=KindScalarVar), dimension(:,:,:,:) :: ff
  real, dimension(:), allocatable:: u, d
  !.............................................................................
  integer:: ix, iy, iz, jy, jz, rec, iv, li(3), ui(3)
  integer(8):: nx, ny, nz
  real(8):: umax, dt
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
    open (io_unit%data, file='data/exa256/t0.02.stag', form='unformatted', &
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
  ff(:,:,:,idx%e) = ff(:,:,:,idx%d)
  do iv=idx%px,idx%pz
    ff(:,:,:,iv) = ff(:,:,:,idx%d)*ff(:,:,:,iv)
  end do
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
  USE io_unit_mod
  class(initial_t):: self
  type(index_t):: idx
  class(mesh_t), pointer, dimension(:):: m
  real(kind=KindScalarVar), dimension(:,:,:,:):: ff
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
    open (io_unit%data, file='data/exa256/t0.02.stag', form='unformatted', &
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
  ff(:,:,:,idx%e) = ff(:,:,:,idx%d)
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

END MODULE
