!===============================================================================
!> Simple initical condition module for tests
!===============================================================================
MODULE initial_conditions_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE mesh_mod
  USE patch_mod
  USE gpatch_mod
  USE omp_mod
  USE index_mod
  USE eos_mod
  implicit none
  private
  real(8), parameter:: pi2=8d0*atan(1d0)
  type, public:: initial_conditions_t
    integer:: seed
    real:: d0
    logical:: mhd
    real(8):: time=0d0
    character(len=64):: solver
    class(gpatch_t), pointer:: patch => null()
    procedure(raw), pointer:: condition=>null()
  contains
    procedure:: pre_init
    procedure:: init
    procedure, nopass:: cleanup
  end type
  real, save, pointer, dimension(:,:,:,:):: buffer => null()
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(initial_conditions_t):: self
  class(gpatch_t), pointer:: patch
  !.............................................................................
  character(len=32), save:: type='void'
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /IC_params/ type
  !----------------------------------------------------------------------------
  call trace_begin('initial_conditions_t%pre_init')
  self%patch => patch
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, IC_params, iostat=iostat)
    call io%namelist_warning ('ramses/hydro/initial_conditions_mod.f90: IC_params', iostat)
    write (io%namelist, IC_params)
 end if
  call io%lock%unset
  select case (trim(type))
  case ('exa256')
    self%condition => exa256
  case ('raw')
    self%condition => raw
  case ('void')
    self%condition => void
  case default
    print*,'UNKNOWN INITIAL CONDITION TYPE :'//trim(type)//':'
  end select
  call trace_end
END SUBROUTINE pre_init

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self)
  class(initial_conditions_t):: self
  !.............................................................................
  character(len=32), save:: type='void'
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /IC_params/ type
  !----------------------------------------------------------------------------
  call trace_begin('initial_conditions_t%init')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, IC_params, iostat=iostat)
    call io%namelist_warning ('ramses/hydro/initial_conditions_mod.f90: IC_params', iostat)
    write (io%namelist, IC_params)
 end if
  call io%lock%unset
  select case (trim(type))
  case ('exa256')
    self%condition => exa256
  case ('raw')
    self%condition => raw
  case ('void')
    self%condition => void
  case default
    print*,'UNKNOWN INITIAL CONDITION TYPE :'//trim(type)//':'
  end select
  call trace_end
END SUBROUTINE init

!===============================================================================
SUBROUTINE cleanup
  if (associated(buffer)) deallocate(buffer)
END SUBROUTINE cleanup

!===============================================================================
!> Void initializer
!===============================================================================
SUBROUTINE void (self)
  class(initial_conditions_t) :: self
END SUBROUTINE void

!===============================================================================
!> Initialize from a raw data file
!===============================================================================
SUBROUTINE exa256 (self)
  class(initial_conditions_t)       :: self
  !.............................................................................
  type(index_t),            pointer :: idx
  class(mesh_t),            pointer :: m(:)
  real, dimension(:,:,:,:), pointer :: ff
  real, dimension(:),   allocatable :: u, d
  integer                           :: ix, iy, iz, jy, jz, rec, iv, li(3), ui(3)
  integer(8)                        :: nx, ny, nz
  real(8)                           :: umax, dt
  logical,                     save :: first_time=.true.
  integer,                     save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_conditions_t::exa256', itimer=itimer)
  m => self%patch%mesh
  ff => self%patch%mem(:,:,:,:,1,1)
  idx => self%patch%idx
  !.............................................................................
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
  li = m%lb
  ui = m%ub
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
    if (idx%bx > 0) then
      ff(li(1):ui(1),jy,jz,idx%bx) = 0.0
      ff(li(1):ui(1),jy,jz,idx%by) = 0.0
      ff(li(1):ui(1),jy,jz,idx%bz) = 0.0
    end if
  end do
  end do
  do iv=idx%px,idx%pz
    ff(:,:,:,iv) = ff(:,:,:,idx%d)*ff(:,:,:,iv)
  end do
  call energy (self)
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
SUBROUTINE raw (self)
  class(initial_conditions_t)       :: self
  !.............................................................................
  type(index_t),            pointer :: idx
  class(mesh_t),            pointer :: m(:)
  real, dimension(:,:,:,:), pointer :: ff
  integer,                parameter :: mv=4
  integer                           :: n(3), i(3), iv, jv(mv)
  logical,                     save :: first_time=.true.
  integer,                     save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t::raw', itimer=itimer)
  m => self%patch%mesh
  ff => self%patch%mem(:,:,:,:,1,1)
  idx => self%patch%idx
  !.............................................................................
  !$omp critical (raw_cr)
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
  call energy (self)
  if (io%verbose > 2) &
    write (io_unit%log,'(a,3i4,1p,2g12.3)') 'ix,iy,iz,dmin,dmax =', i, &
      minval(ff(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui,1)), &
      maxval(ff(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui,1))
  !$omp end critical (raw_cr)
  call trace_end (itimer)
END SUBROUTINE raw

!===============================================================================
!> Calculate total energy
!===============================================================================
SUBROUTINE energy (self)
  class(initial_conditions_t)       :: self
  !.............................................................................
  type(index_t),            pointer :: idx
  class(mesh_t),            pointer :: m(:)
  real, dimension(:,:,:,:), pointer :: ff
  real(8)                           :: fgamma
  !-----------------------------------------------------------------------------
  call trace%begin ('initial_conditions_t%energy')
  if (eos%gamma==1d0) then
    fgamma = 1d0
  else
    fgamma = eos%gamma*(eos%gamma-1d0)
  end if
  ff(:,:,:,idx%e) = 0.5*(ff(:,:,:,idx%px)**2+ &
                         ff(:,:,:,idx%py)**2+ &
                         ff(:,:,:,idx%pz)**2)/ff(:,:,:,idx%d) &
                  + ff(:,:,:,idx%e)
  if (idx%bx > 0) then
    ff(:,:,:,idx%e) = ff(:,:,:,idx%e) + &
      0.5*(ff(:,:,:,idx%bx)**2 + &
           ff(:,:,:,idx%by)**2 + &
           ff(:,:,:,idx%bz)**2)
  end if
  !print *, 'IIII', idx%e, fgamma, &
    !minval(ff(m(1)%lb:m(1)%ub,m(2)%lb:m(2)%ub,m(3)%lb:m(3)%ub,idx%d)), &
    !maxval(ff(m(1)%lb:m(1)%ub,m(2)%lb:m(2)%ub,m(3)%lb:m(3)%ub,idx%d)), &
    !minval(ff(m(1)%lb:m(1)%ub,m(2)%lb:m(2)%ub,m(3)%lb:m(3)%ub,idx%e)), &
    !maxval(ff(m(1)%lb:m(1)%ub,m(2)%lb:m(2)%ub,m(3)%lb:m(3)%ub,idx%e))
  call trace%end()
END SUBROUTINE energy

END MODULE
