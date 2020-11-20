!===============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE random_mod
  USE index_mod
  implicit none
  private
  real(8), parameter:: pi2=8d0*atan(1d0)
  type, public:: initial_t
    integer:: seed
    logical:: mhd
    real(8):: a(3), b(3), phi(3), d0, e0
    real(8):: a0(3), k(3), b0(3), u0(3)
    real(8):: time
    procedure(raw), pointer:: condition=>null()
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
SUBROUTINE init (self)
   class(initial_t):: self
   integer, save:: seed=-33
   logical, save:: mhd=.true.
   real, save:: a0(3)=0.1, k(3)=1.0, d0=1.0, e0=1.0, b0(3)=0.0, u0(3)=0.0
   character(len=32), save:: type='single_solenoidal'
   namelist /IC_params/ seed, k, a0, u0, b0, d0, e0, mhd, type
   !----------------------------------------------------------------------------
   call trace_begin('initial_t%init')
   call io%lock%set
   if (first_time) then
     first_time = .false.
     rewind (io%input)
     read (io%input, initial_params)
     mhd = any(b0 /= 0.0)
     if (io%master) write (*, initial_params)
   end if
   call io%lock%unset
   self%seed   = seed
   self%k      = k
   self%a0     = a0
   self%b0     = b0
   self%u0     = u0
   self%d0     = d0
   self%e0     = e0
   self%mhd    = mhd
   self%time   = 0d0
   select case (type)
   case ('single_solenoidal')
     self%condition => single_solenoidal
   case ('single_compressive')
     self%condition => single_compressive
   case ('single_advection')
     self%condition => single_advection
   case ('raw')
     self%condition => raw
   case default
     print*,'UNKNOWN INITIAL CONDITION TYPE'
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
FUNCTION single_solenoidal (self, m, idx) RESULT (ff)
  class(initial_t)      :: self
  type(index_t)         :: idx
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,8) :: ff
  !.............................................................................
  type(random_t)        :: random
  integer               :: ix, iy, iz
  real(8)               :: kx, ky, kz
  type(index_t)         :: idx
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:single_solenoidal')
  self%phi = 0.0
  do iz=m(3)%lb,m(3)%ub
    kz   =  pi2*self%k(3)*(m(3)%r(iz)+m(3)%p)/m(3)%b + self%phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky   =  pi2*self%k(2)*(m(2)%r(iy)+m(2)%p)/m(2)%b + self%phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*self%k(1)*(m(1)%r(ix)+m(1)%p)/m(1)%b + self%phi(1)
        ff(ix,iy,iz,idx%d) = self%d0
        ff(ix,iy,iz,idx%e) = self%e0
        ff(ix,iy,iz,idx%px) = self%a0(1)*(sin(ky)+sin(kz))
        ff(ix,iy,iz,idx%py) = self%a0(2)*(sin(kz)+sin(kx))
        ff(ix,iy,iz,idx%pz) = self%a0(3)*(sin(kx)+sin(ky))
        ff(ix,iy,iz,idx%bx) = self%b0(1)*(sin(kx)+sin(ky))
        ff(ix,iy,iz,idx%by) = self%b0(2)*(sin(kx)+sin(ky))
        ff(ix,iy,iz,idx%bz) = self%b0(3)*(sin(kx)+sin(ky))
      end do
    end do
  end do
  call trace_end
END FUNCTION single_solenoidal

!===============================================================================
!> Initialize a single 3-D wavenumber, with random amplitude and phase
!===============================================================================
FUNCTION single_compressive (self, m, idx) RESULT (ff)
  class(initial_t)        :: self
  type(index_t):: idx
  class(mesh_t), pointer, dimension(:):: m
  type(index_t)           :: idx
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,8) :: ff
  !.............................................................................
  type(random_t)        :: random
  integer               :: ix, iy, iz
  real(8)               :: kx, ky, kz
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:single_compressive')
  self%phi = 0.0
  do iz=m(3)%lb,m(3)%ub
    kz   =  pi2*self%k(3)*(m(3)%r(iz)+m(3)%p)/m(3)%b + self%phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky   =  pi2*self%k(2)*(m(2)%r(iy)+m(2)%p)/m(2)%b + self%phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*self%k(1)*(m(1)%r(ix)+m(1)%p)/m(1)%b + self%phi(1)
        ff(ix,iy,iz,idx%d) = self%d0
        ff(ix,iy,iz,idx%e) = self%e0
        ff(ix,iy,iz,idx%px) = self%u0(1)*sin(kx)
        ff(ix,iy,iz,idx%py) = self%u0(2)*sin(ky)
        ff(ix,iy,iz,idx%pz) = self%u0(3)*sin(kz)
        ff(ix,iy,iz,idx%bx) = self%b0(1)
        ff(ix,iy,iz,idx%by) = self%b0(2)
        ff(ix,iy,iz,idx%bz) = self%b0(3)
      end do
    end do
  end do
  call trace_end
END FUNCTION single_compressive

!===============================================================================
!> Initialize a single 3-D wavenumber, with random amplitude and phase
!===============================================================================
FUNCTION single_advection (self, m, idx) RESULT (ff)
  class(initial_t)        :: self
  type(index_t)           :: idx
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,8) :: ff
  !.............................................................................
  type(random_t)        :: random
  integer               :: ix, iy, iz
  real(8)               :: kx, ky, kz
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:single_compressive')
  self%phi = 0.0
  do iz=m(3)%lb,m(3)%ub
    kz   =  pi2*self%k(3)*(m(3)%r(iz)+m(3)%p)/m(3)%b + self%phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky   =  pi2*self%k(2)*(m(2)%r(iy)+m(2)%p)/m(2)%b + self%phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*self%k(1)*(m(1)%r(ix)+m(1)%p)/m(1)%b + self%phi(1)
        ff(ix,iy,iz,idx%d) = self%d0*(1.+self%a0(1)*sin(kx)+self%a0(2)*sin(ky)+self%a0(3)*sin(kz))
        ff(ix,iy,iz,idx%e) = self%e0
        ff(ix,iy,iz,idx%px) = self%d0*self%u0(1)
        ff(ix,iy,iz,idx%py) = self%d0*self%u0(2)
        ff(ix,iy,iz,idx%pz) = self%d0*self%u0(3)
        ff(ix,iy,iz,idx%bx) = self%b0(3)*sin(kz)
        ff(ix,iy,iz,idx%by) = self%b0(1)*sin(kx)
        ff(ix,iy,iz,idx%bz) = self%b0(2)*sin(ky)
      end do
    end do
  end do
  call trace_end
END FUNCTION single_advection

!===============================================================================
!> Initialize from a raw data file
!===============================================================================
FUNCTION raw (self, m, idx) RESULT (ff)
  USE io_unit_mod
  class(initial_t):: self
  type(index_t):: idx
  class(mesh_t), pointer, dimension(:) :: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,8):: ff
  !.............................................................................
  integer:: n(3), i(3), iv
  integer, save:: mv=4
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t::raw')
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
    open (io_unit%data, file='data/ramses256/t0.02.raw', form='unformatted', &
      access='direct', recl=4*product(n))
    do iv=1,mv
      read (io_unit%data, rec=iv) buffer(:,:,:,iv)
    end do
    close (io_unit%data)
  end if
  i = 1 + nint(m%p/m%d)-m%n/2                            ! index offsets
  ff = 0.0
  do iv=1,mv
    ff(m(1)%li:m(1)%ui, &
       m(2)%li:m(2)%ui, &
       m(3)%li:m(3)%ui,iv) = buffer(i(1):i(1)+m(1)%n-1, &
                                    i(2):i(2)+m(2)%n-1, &
                                    i(3):i(3)+m(3)%n-1,iv)
  end do
  do iv=2,4
    ff(:,:,:,iv) = ff(:,:,:,idx%d)*ff(:,:,:,iv)
  end do
  if (io%verbose > 2) &
    write (io_unit%log,'(a,3i4,1p,2g12.3)') 'ix,iy,iz,dmin,dmax =', i, &
      minval(ff(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui,idx%d)), &
      maxval(ff(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui,idx%d))
  call trace_end
END FUNCTION raw

END MODULE
