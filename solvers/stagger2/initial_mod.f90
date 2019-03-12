!===============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE vector_mod
  USE mesh_mod
  USE random_mod
  USE stagger_mod
  USE index_mod
  implicit none
  private
  real(8), parameter:: pi2=8d0*atan(1d0)
  type, public:: initial_t
    integer:: id=0
    integer:: seed
    logical:: mhd
    real(8):: gamma=0d0
    real(8):: a(3), b(3), phi(3), d0, p0, e0
    real(8):: a0(3), k(3), b0(3), u0(3)
    real(8):: time
    character(len=64):: solver
    procedure(single_solenoidal), pointer:: condition=>null()
  contains
    procedure:: init
    procedure, nopass:: cleanup
  end type
  logical, save:: first_time=.true.
  real, save, pointer, dimension(:,:,:,:):: buffer => null()
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self, solver, gamma)
   class(initial_t):: self
   character(len=64):: solver
   real, optional:: gamma
   !............................................................................
   integer:: iostat
   integer, save:: seed=-33
   logical, save:: mhd=.true.
   real, save:: a0(3)=0.1, k(3)=1.0, d0=1.0, p0=1.0, b0(3)=0.0, u0(3)=0.0
   character(len=32), save:: type='void'
   namelist /initial_params/ type, seed, k, a0, u0, b0, d0, p0, mhd
   !----------------------------------------------------------------------------
   call trace_begin('initial_t%init')
   if (present(gamma)) self%gamma = gamma
   !$omp critical (input_cr)
   if (first_time) then
     first_time = .false.
     mhd = self%mhd
     rewind (io%input)
     read (io%input, initial_params, iostat=iostat)
     mhd = mhd .or. any(b0 /= 0.0)
     if (io%master) write (*, initial_params)
   end if
   self%solver = trim(solver)
   !$omp end critical (input_cr)
   self%seed   = seed
   self%k      = k
   self%a0     = a0
   self%b0     = b0
   self%u0     = u0
   self%d0     = d0
   self%p0     = p0
   self%gamma  = max(self%gamma,1d0+1d-5)
   self%e0     = p0/(self%gamma-1d0)
   self%mhd    = mhd
   self%time   = 0d0
   select case (type)
   case ('void')
     self%condition => void
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
     print*,'UNKNOWN INITIAL CONDITION TYPE'
   end select
   call trace_end
END SUBROUTINE init

!===============================================================================
!> Do-nothing IC; when initializing in experiment_t%init()
!===============================================================================
SUBROUTINE void (self, m, ff, idx)
  class(initial_t)      :: self
  type(index_t):: idx
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(:,:,:,:) :: ff
END SUBROUTINE void

!===============================================================================
!> Initialize a single 3-D wavenumber, with random amplitude and phase
!===============================================================================
SUBROUTINE single_solenoidal (self, m, ff, idx)
  class(initial_t)      :: self
  type(index_t):: idx
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(:,:,:,:) :: ff
  !.............................................................................
  type(random_t)        :: random
  integer               :: ix, iy, iz
  real(8)               :: kx, ky, kz
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:single_solenoidal')
  self%phi = 0.0
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    kz   =  pi2*self%k(3)*(r3(iz)+m(3)%p-0.5d0*m(3)%s)/m(3)%b + self%phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky   =  pi2*self%k(2)*(r2(iy)+m(2)%p-0.5d0*m(2)%s)/m(2)%b + self%phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*self%k(1)*(r1(ix)+m(1)%p-0.5d0*m(1)%s)/m(1)%b + self%phi(1)
        ff(ix,iy,iz,idx%d ) = self%d0*(1.0 + self%a0(1)*sin(kx) + self%a0(2)*sin(ky) + self%a0(3)*sin(kz))
        ff(ix,iy,iz,idx%e ) = self%e0*(1.0 + self%a0(1)*sin(kx) + self%a0(2)*sin(ky) + self%a0(3)*sin(kz))
        ff(ix,iy,iz,idx%px) = self%u0(1)*(sin(ky)+sin(kz))
        ff(ix,iy,iz,idx%py) = self%u0(2)*(sin(kz)+sin(kx))
        ff(ix,iy,iz,idx%pz) = self%u0(3)*(sin(kx)+sin(ky))
        if (self%mhd) then
          ff(ix,iy,iz,idx%bx) = self%b0(1)*(sin(ky)+sin(kz))
          ff(ix,iy,iz,idx%by) = self%b0(2)*(sin(kz)+sin(kx))
          ff(ix,iy,iz,idx%bz) = self%b0(3)*(sin(kx)+sin(ky))
        end if
      end do
    end do
  end do
  end associate

  if (self%solver(1:9) == 'stagger2_') then
    associate (e=>ff(:,:,:,idx%e), d=>ff(:,:,:,idx%d))
    call e2s(e, d, self%gamma)
    end associate
  end if
  call trace_end
END SUBROUTINE single_solenoidal

!===============================================================================
!> Initialize a single 3-D wavenumber, with random amplitude and phase
!===============================================================================
SUBROUTINE single_compressive (self, m, ff, idx)
  class(initial_t)        :: self
  type(index_t):: idx
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(:,:,:,:) :: ff
  !.............................................................................
  type(random_t)        :: random
  integer               :: ix, iy, iz
  real(8)               :: kx, ky, kz
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:single_compressive')
  self%phi = 0.0
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    kz   =  pi2*self%k(3)*(r3(iz)+m(3)%p-0.5d0*m(3)%s)/m(3)%b + self%phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky   =  pi2*self%k(2)*(r2(iy)+m(2)%p-0.5d0*m(2)%s)/m(2)%b + self%phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*self%k(1)*(r1(ix)+m(1)%p-0.5d0*m(1)%s)/m(1)%b + self%phi(1)
        ff(ix,iy,iz,idx%d ) = self%d0
        ff(ix,iy,iz,idx%e ) = self%e0
        ff(ix,iy,iz,idx%px) = self%u0(1)*sin(kx)
        ff(ix,iy,iz,idx%py) = self%u0(2)*sin(ky)
        ff(ix,iy,iz,idx%pz) = self%u0(3)*sin(kz)
        if (self%mhd) then
          ff(ix,iy,iz,idx%bx) = self%b0(1)
          ff(ix,iy,iz,idx%by) = self%b0(2)
          ff(ix,iy,iz,idx%bz) = self%b0(3)
        end if
      end do
    end do
  end do
  end associate

  if (self%solver(1:9) == 'stagger2_') then
    associate (e=>ff(:,:,:,idx%e), d=>ff(:,:,:,idx%d))
    call e2s(e, d, self%gamma)
    end associate
  end if
  call trace_end
END SUBROUTINE single_compressive

!===============================================================================
!> Initialize a single 3-D wavenumber, with random amplitude and phase
!===============================================================================
SUBROUTINE single_advection (self, m, ff, idx)
  class(initial_t)        :: self
  type(index_t):: idx
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(:,:,:,:) :: ff
  !.............................................................................
  type(random_t)        :: random
  integer               :: ix, iy, iz
  real(8)               :: kx, ky, kz
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:single_compressive')
  self%phi = 0.0
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    kz   =  pi2*self%k(3)*(r3(iz)+m(3)%p-0.5d0*m(3)%s)/m(3)%b + self%phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky   =  pi2*self%k(2)*(r2(iy)+m(2)%p-0.5d0*m(2)%s)/m(2)%b + self%phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*self%k(1)*(r1(ix)+m(1)%p-0.5d0*m(1)%s)/m(1)%b + self%phi(1)
        ff(ix,iy,iz,idx%d ) = self%d0*(1.+self%a0(1)*sin(kx)+self%a0(2)*sin(ky)+self%a0(3)*sin(kz))
        ff(ix,iy,iz,idx%e ) = self%e0
        ff(ix,iy,iz,idx%px) = self%d0*self%u0(1)
        ff(ix,iy,iz,idx%py) = self%d0*self%u0(2)
        ff(ix,iy,iz,idx%pz) = self%d0*self%u0(3)
        if (self%mhd) then
          ff(ix,iy,iz,idx%bx) = self%b0(3)*sin(kz)
          ff(ix,iy,iz,idx%by) = self%b0(1)*sin(kx)
          ff(ix,iy,iz,idx%bz) = self%b0(2)*sin(ky)
        end if
      end do
    end do
  end do
  end associate

  if (self%solver(1:9) == 'stagger2_') then
    associate (e=>ff(:,:,:,idx%e), d=>ff(:,:,:,idx%d))
    call e2s(e, d, self%gamma)
    end associate
  end if
  call trace_end
END SUBROUTINE single_advection

!===============================================================================
!> Initialize from a raw data file
!===============================================================================
SUBROUTINE exa256 (self, m, ff, idx)
  USE io_unit_mod
  class(initial_t)      :: self
  type(index_t):: idx
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(:,:,:,:) :: ff
  real, dimension(:), allocatable:: u, d
  !.............................................................................
  integer               :: ix, iy, iz, jy, jz, rec
  integer(8)            :: nx, ny, nz
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t::exa256')
  self%mhd = .false.
  ff(:,:,:,:) = 0.0
  allocate (u(m(1)%n))
  allocate (d(m(1)%n))
  open (io_unit%data, file='data/exa256/t0.02.stag', form='unformatted', &
    access='direct', recl=m(1)%n*4)
  self%time = 0.02
  ix=m(1)%p/m(1)%d+0.5-m(1)%n/2
  iy=m(2)%p/m(2)%d+0.5-m(2)%n/2
  iz=m(3)%p/m(3)%d+0.5-m(3)%n/2
  nx = m(1)%b/m(1)%d + 0.5
  ny = m(2)%b/m(2)%d + 0.5
  nz = m(3)%b/m(3)%d + 0.5
  do jz=m(3)%li,m(3)%ui
  do jy=m(2)%li,m(2)%ui
    rec = 1   + (ix           )      /m(1)%n
    rec = rec + (iy+jy-m(2)%li)*nx   /m(1)%n
    rec = rec + (iz+jz-m(3)%li)*nx*ny/m(1)%n
    !print*,ix,iy,iz,rec
    read (io_unit%data,rec=rec) d; ff(m(1)%li:m(1)%ui,jy,jz,idx%d) = d
                                   ff(m(1)%li:m(1)%ui,jy,jz,idx%e) = d
    rec = rec + nx*ny*nz/m(1)%n
    read (io_unit%data,rec=rec) u; ff(m(1)%li:m(1)%ui,jy,jz,idx%px) = u
    rec = rec + nx*ny*nz/m(1)%n
    read (io_unit%data,rec=rec) u; ff(m(1)%li:m(1)%ui,jy,jz,idx%py) = u
    rec = rec + nx*ny*nz/m(1)%n
    read (io_unit%data,rec=rec) u; ff(m(1)%li:m(1)%ui,jy,jz,idx%pz) = u
  end do
  end do
  if (io%verbose > 0) &
    print *,'ix,iy,iz,dmin,dmax =', ix, iy, iz, &
      minval(ff(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui,1)), &
      maxval(ff(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui,1))

  if (self%solver(1:9) == 'stagger2_') then
    associate (e=>ff(:,:,:,idx%e), d=>ff(:,:,:,idx%d))
    call e2s(e, d, self%gamma)
    end associate
  end if
  call trace_end
END SUBROUTINE exa256

!===============================================================================
!> Initialize from a raw data file
!===============================================================================
SUBROUTINE raw (self, m, ff, idx)
  USE io_unit_mod
  class(initial_t):: self
  type(index_t):: idx
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(:,:,:,:):: ff
  !.............................................................................
  integer:: n(3), i(3), iv, iw, iostat
  integer, save:: mv=4
  logical, save:: first_time=.true.
  character(len=64), save:: filename = 'data/exa256/t0.02.stag'
  namelist /raw_params/ filename
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t::raw')
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
    rewind (io%input); read (io%input,raw_params,iostat=iostat)
    if (io%master) write (*,raw_params)
    n = nint(m%b/m%d)
    allocate (buffer(n(1),n(2),n(3),mv))
    open (io_unit%data, file=filename, form='unformatted', &
      access='direct', recl=4*product(n))
    do iv=1,mv
      read (io_unit%data, rec=iv) buffer(:,:,:,iv)
    end do
    close (io_unit%data)
  end if
  i = 1 + nint(m%p/m%d) - m%n/2                         ! index offsets
  ff = 0.0
  do iw=1,mv
    iv = merge(1,iw+1,iw==1)
    ff(m(1)%li:m(1)%ui, &
       m(2)%li:m(2)%ui, &
       m(3)%li:m(3)%ui,iv) = buffer(i(1):i(1)+m(1)%n-1, &
                                    i(2):i(2)+m(2)%n-1, &
                                    i(3):i(3)+m(3)%n-1,iw)
  end do
  ff(:,:,:,idx%e) = ff(:,:,:,idx%d)
  if (io%verbose > 2) &
    write (io_unit%log,'(a,3i4,1p,2g12.3)') 'ix,iy,iz,dmin,dmax =', i, &
      minval(ff(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui,1)), &
      maxval(ff(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui,1))
  !$omp end critical (raw_cr)

  if (self%solver(1:9) == 'stagger2_') then
    associate (e=>ff(:,:,:,idx%e), d=>ff(:,:,:,idx%d))
    call e2s(e, d, self%gamma)
    end associate
  end if
  call trace_end
END SUBROUTINE raw

!===============================================================================
SUBROUTINE cleanup
  if (associated(buffer)) deallocate(buffer)
END SUBROUTINE cleanup

!=======================================================================
!> Get entropy per unit volume from energy per unit volume
!> Note that energy and entropy are assumed to occupy the same memory
!=======================================================================
SUBROUTINE e2s(s, d, gamma)
  real(8):: gamma
  !.....................................................................
  real, dimension(:,:,:):: s, d
  real:: g1
  !---------------------------------------------------------------------
  g1 = gamma-1.0
  s = d*log(s*g1/d**gamma)/g1
END SUBROUTINE e2s

END MODULE
