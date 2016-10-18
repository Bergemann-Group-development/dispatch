!> =============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!> =============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE random_mod
  USE hydro_parameters
  USE omp_mod
  implicit none
  private
  real(8), parameter:: pi2=8d0*atan(1d0)
  type, public:: initial_t
    integer:: seed
    logical:: mhd
    real(8):: a(3), b(3), phi(3), d0, s0
    real(8):: a0(3), k(3), b0(3), u0(3)
    real(8):: time
    procedure(single_solenoidal), pointer:: selected=>null()
  contains
    procedure:: init
    procedure, nopass:: cleanup
  end type
  logical, save:: first_time=.true.
  real, save, pointer, dimension(:,:,:,:):: buffer => null()
CONTAINS

!> =============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!> =============================================================================
SUBROUTINE init (self)
   class(initial_t):: self
   integer, save:: seed=-33
   logical, save:: mhd=.true.
   real, save:: a0(3)=0.1, k(3)=1.0, d0=1.0, s0=0.0, b0(3)=0.0, u0(3)=0.0
   character(len=32), save:: type='single_solenoidal'
   namelist /initial_params/ seed, k, a0, u0, b0, d0, s0, mhd, type
   !----------------------------------------------------------------------------
   call trace_begin('initial_t%init')
   !$omp critical (input_cr)
   if (first_time) then
     first_time = .false.
     rewind (io%input)
     read (io%input, initial_params)
     mhd = any(b0 /= 0.0)
     if (io%master) write (*, initial_params)
   end if
   !$omp end critical (input_cr)
   self%seed   = seed
   self%k      = k
   self%a0     = a0
   self%b0     = b0
   self%u0     = u0
   self%d0     = d0
   self%s0     = s0
   self%mhd    = mhd
   self%time   = 0d0
   select case (type)
   case ('single_solenoidal')
     self%selected => single_solenoidal
   case ('single_compressive')
     self%selected => single_compressive
   case ('single_advection')
     self%selected => single_advection
   case ('exa256')
     self%selected => exa256
   case ('raw')
     self%selected => raw
   case default
     print*,'UNKNOWN INITIAL CONDITION TYPE :'//trim(type)//':'
     stop
   end select
   call trace_end
END SUBROUTINE init

!> =============================================================================
SUBROUTINE cleanup
  if (associated(buffer)) deallocate(buffer)
END SUBROUTINE cleanup

!> =============================================================================
!> Initialize a single 3-D wavenumber, with random amplitude and phase
!> =============================================================================
SUBROUTINE single_solenoidal (self, m, ff)
  class(initial_t)      :: self
  class(mesh_t)          :: m(3)
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
        ff(ix,iy,iz,1) = self%d0
        ff(ix,iy,iz,2) = self%a0(1)*(sin(ky)+sin(kz))
        ff(ix,iy,iz,3) = self%a0(2)*(sin(kz)+sin(kx))
        ff(ix,iy,iz,4) = self%a0(3)*(sin(kx)+sin(ky))
        ff(ix,iy,iz,5) = self%d0*csound**2/(gamma*(gamma-1)) + &
          half*ff(ix,iy,iz,1)*(ff(ix,iy,iz,2)**2+ff(ix,iy,iz,3)**2+ff(ix,iy,iz,4)**2)
        if (self%mhd) then
          ff(ix,iy,iz,6) = self%b0(1)*(sin(kx)+sin(ky))
          ff(ix,iy,iz,7) = self%b0(2)*(sin(kx)+sin(ky))
          ff(ix,iy,iz,8) = self%b0(3)*(sin(kx)+sin(ky))
        end if
      end do
    end do
  end do
  end associate
  call trace_end
END SUBROUTINE single_solenoidal

!> =============================================================================
!> Initialize a single 3-D wavenumber, with non-random amplitude and zero phase
!> =============================================================================
SUBROUTINE single_compressive (self, m, ff)
  class(initial_t)      :: self
  class(mesh_t)          :: m(3)
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
        ff(ix,iy,iz,1) = self%d0
        ff(ix,iy,iz,2) = self%u0(1)*sin(kx)
        ff(ix,iy,iz,3) = self%u0(2)*sin(ky)
        ff(ix,iy,iz,4) = self%u0(3)*sin(kz)
        ff(ix,iy,iz,5) = self%d0*csound**2/(gamma*(gamma-1)) + &
          half*ff(ix,iy,iz,1)*(ff(ix,iy,iz,2)**2+ff(ix,iy,iz,3)**2+ff(ix,iy,iz,4)**2)
        if (self%mhd) then
          ff(ix,iy,iz,6) = self%b0(1)
          ff(ix,iy,iz,7) = self%b0(2)
          ff(ix,iy,iz,8) = self%b0(3)
        end if
      end do
    end do
  end do
  end associate
  call trace_end
END SUBROUTINE single_compressive

!> =============================================================================
!> Initialize a single 3-D wavenumber, with random amplitude and phase
!> =============================================================================
SUBROUTINE single_advection (self, m, ff)
  class(initial_t)        :: self
  class(mesh_t)          :: m(3)
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
        ff(ix,iy,iz,1) = self%d0*(1.+self%a0(1)*sin(kx)+self%a0(2)*sin(ky)+self%a0(3)*sin(kz))
        ff(ix,iy,iz,2) = self%d0*self%u0(1)
        ff(ix,iy,iz,3) = self%d0*self%u0(2)
        ff(ix,iy,iz,4) = self%d0*self%u0(3)
        ff(ix,iy,iz,5) = self%d0*csound**2/(gamma*(gamma-1)) + &
          half*ff(ix,iy,iz,1)*(ff(ix,iy,iz,2)**2+ff(ix,iy,iz,3)**2+ff(ix,iy,iz,4)**2)
        if (self%mhd) then
          ff(ix,iy,iz,6) = self%b0(3)*sin(kz)
          ff(ix,iy,iz,7) = self%b0(1)*sin(kx)
          ff(ix,iy,iz,8) = self%b0(2)*sin(ky)
        end if
      end do
    end do
  end do
  end associate
  call trace_end
END SUBROUTINE single_advection

!> =============================================================================
!> Initialize from a raw data file
!> =============================================================================
SUBROUTINE exa256 (self, m, ff)
  USE io_unit_mod
  class(initial_t):: self
  class(mesh_t):: m(3)
  real, dimension(:,:,:,:):: ff
  real, dimension(:), allocatable:: u, d
  !.............................................................................
  integer:: ix, iy, iz, jy, jz, rec, i, li(3), ui(3)
  integer(8):: nx, ny, nz
  real(dp):: umax, dt
  logical, save:: first_time=.true.
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t::exa256')
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
    read (io_unit%data,rec=rec) d; ff(m(1)%li:m(1)%ui,jy,jz,1) = d
    rec = rec + nx*ny*nz/m(1)%n
    read (io_unit%data,rec=rec) u; ff(m(1)%li:m(1)%ui,jy,jz,2) = u
    rec = rec + nx*ny*nz/m(1)%n
    read (io_unit%data,rec=rec) u; ff(m(1)%li:m(1)%ui,jy,jz,3) = u
    rec = rec + nx*ny*nz/m(1)%n
    read (io_unit%data,rec=rec) u; ff(m(1)%li:m(1)%ui,jy,jz,4) = u
    umax=max(umax,3.+maxval(abs(ff(li(1):ui(1),jy,jz,2)) &
                           +abs(ff(li(1):ui(1),jy,jz,3)) &
                           +abs(ff(li(1):ui(1),jy,jz,4))))
  end do
  end do
  ff(:,:,:,5) = 0.5*ff(:,:,:,1)*(ff(:,:,:,2)**2+ &
                                 ff(:,:,:,3)**2+ &
                                 ff(:,:,:,4)**2) &
              + ff(:,:,:,1)/(gamma-one)
  do i=2,4
    ff(:,:,:,i) = ff(:,:,:,1)*ff(:,:,:,i)
  end do
  dt = courant_factor*m(1)%d/umax
  if (io%verbose > 2) write (io_unit%log,*) 'initial dt:',dt,courant_factor,gamma,umax
  if (io%verbose > 2) &
    write (io_unit%log,'(a,3i4,1p,2g12.3)') 'ix,iy,iz,dmin,dmax =', ix, iy, iz, &
      minval(ff(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui,1)), &
      maxval(ff(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui,1))
  call trace_end
END SUBROUTINE exa256

!> =============================================================================
!> Initialize from a raw data file
!> =============================================================================
SUBROUTINE raw (self, m, ff)
  USE io_unit_mod
  class(initial_t):: self
  class(mesh_t):: m(3)
  real, dimension(:,:,:,:):: ff
  !.............................................................................
  integer:: n(3), i(3), iv
  integer, save:: mv=4
  logical, save:: first_time=.true.
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
    n = nint(m%b/m%d)
    allocate (buffer(n(1),n(2),n(3),mv))
    open (io_unit%data, file='data/exa256/t0.02.raw', form='unformatted', &
      access='direct', recl=4*product(n))
    do iv=1,mv
      read (io_unit%data, rec=iv) buffer(:,:,:,iv)
    end do
    close (io_unit%data)
  end if
  i = 1 + nint(m%p/m%d) - m%n/2                          ! index offsets
  ff = 0.0
  do iv=1,mv
    ff(m(1)%li:m(1)%ui, &
       m(2)%li:m(2)%ui, &
       m(3)%li:m(3)%ui,iv) = buffer(i(1):i(1)+m(1)%n-1, &
                                    i(2):i(2)+m(2)%n-1, &
                                    i(3):i(3)+m(3)%n-1,iv)
  end do
  do iv=2,4
    ff(:,:,:,iv) = ff(:,:,:,1)*ff(:,:,:,iv)
  end do
  if (io%verbose > 2) &
    write (io_unit%log,'(a,3i4,1p,2g12.3)') 'ix,iy,iz,dmin,dmax =', i, &
      minval(ff(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui,1)), &
      maxval(ff(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui,1))
  !$omp end critical (raw_cr)
  call trace_end
END SUBROUTINE raw

END MODULE
