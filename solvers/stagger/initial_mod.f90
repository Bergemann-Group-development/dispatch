!> =============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!> =============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE vector_mod
  USE mesh_mod
  USE random_mod
  implicit none
  private
  real(8), parameter:: pi2=8d0*atan(1d0)
  type, public:: initial_t
    integer:: seed
    logical:: mhd
    real(8):: a(3), b(3), phi(3), d0, s0
    real(8):: a0(3), k(3), b0(3), u0(3)
    procedure(single_solenoidal), pointer:: selected=>null()
  contains
    procedure init
  end type
  logical, save:: first_time=.true.
CONTAINS

!> =============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!> =============================================================================
SUBROUTINE init (self)
   class(initial_t):: self
   integer, save:: seed=-33
   logical, save:: mhd=.true.
   real, save:: a0(3)=0.1, k(3)=1.0, d0=1.0, s0=0.0, b0(3)=0.0, u0(3)=0.0
   character(len=32):: type='single_solenoidal'
   namelist /initial_params/ seed, k, a0, u0, b0, d0, s0, mhd, type
   !----------------------------------------------------------------------------
   call trace_begin('initial_t%init')
   if (first_time) then
    !$omp critical (experiment_cr)
     if (first_time) then
       first_time = .false.
       rewind (io%input)
       read (io%input, initial_params)
       mhd = any(self%b0 /= 0.0)
       if (io%master) write (*, initial_params)
     end if
    !$omp end critical (experiment_cr)
   end if
   self%seed   = seed
   self%k      = k
   self%a0     = a0
   self%b0     = b0
   self%u0     = u0
   self%d0     = d0
   self%s0     = s0
   self%mhd    = mhd
   select case (type)
   case ('single_solenoidal')
     self%selected => single_solenoidal
   case ('single_compressive')
     self%selected => single_compressive
   case ('single_advection')
     self%selected => single_advection
   case default
     print*,'UNKNOWN INITIAL CONDITION TYPE'
     stop
   end select
   call trace_end
END SUBROUTINE init

!> =============================================================================
!> Initialize a single 3-D wavenumber, with random amplitude and phase
!> =============================================================================
FUNCTION single_solenoidal (self, m) RESULT (ff)
  class(initial_t)      :: self
  type(mesh_t)          :: m(3)
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,8) :: ff
  !.............................................................................
  type(random_t)        :: random
  integer               :: ix, iy, iz
  real(8)               :: kx, ky, kz
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:single_solenoidal')
  if (first_time) call self%init
  self%phi = 0.0
  do iz=m(3)%lb,m(3)%ub
    kz   =  pi2*self%k(3)*(m(3)%r(iz)+m(3)%p)/m(3)%b + self%phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky   =  pi2*self%k(2)*(m(2)%r(iy)+m(2)%p)/m(2)%b + self%phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*self%k(1)*(m(1)%r(ix)+m(1)%p)/m(1)%b + self%phi(1)
        ff(ix,iy,iz,1) = self%d0
        ff(ix,iy,iz,2) = self%s0
        ff(ix,iy,iz,3) = self%a0(1)*(sin(ky)+sin(kz))
        ff(ix,iy,iz,4) = self%a0(2)*(sin(kz)+sin(kx))
        ff(ix,iy,iz,5) = self%a0(3)*(sin(kx)+sin(ky))
        ff(ix,iy,iz,6) = self%b0(1)*(sin(kx)+sin(ky))
        ff(ix,iy,iz,7) = self%b0(2)*(sin(kx)+sin(ky))
        ff(ix,iy,iz,8) = self%b0(3)*(sin(kx)+sin(ky))
      end do
    end do
  end do
  do ix=m(1)%lb,m(1)%ub
    print'(a,i4,f22.16)', 'initial:', ix, ff(ix,1,1,3)
  end do
  call trace_end
END FUNCTION single_solenoidal

!> =============================================================================
!> Initialize a single 3-D wavenumber, with random amplitude and phase
!> =============================================================================
FUNCTION single_compressive (self, m) RESULT (ff)
  class(initial_t)        :: self
  type(mesh_t)          :: m(3)
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,8) :: ff
  !.............................................................................
  type(random_t)        :: random
  integer               :: ix, iy, iz
  real(8)               :: kx, ky, kz
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:single_compressive')
  if (first_time) call self%init
  print*,self%k,self%a0
  self%phi = 0.0
  do iz=m(3)%lb,m(3)%ub
    kz   =  pi2*self%k(3)*(m(3)%r(iz)+m(3)%p)/m(3)%b + self%phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky   =  pi2*self%k(2)*(m(2)%r(iy)+m(2)%p)/m(2)%b + self%phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*self%k(1)*(m(1)%r(ix)+m(1)%p)/m(1)%b + self%phi(1)
        ff(ix,iy,iz,1) = self%d0
        ff(ix,iy,iz,2) = self%s0
        ff(ix,iy,iz,3) = self%u0(1)*sin(kx)
        ff(ix,iy,iz,4) = self%u0(2)*sin(ky)
        ff(ix,iy,iz,5) = self%u0(3)*sin(kz)
        ff(ix,iy,iz,6) = self%b0(1)
        ff(ix,iy,iz,7) = self%b0(2)
        ff(ix,iy,iz,8) = self%b0(3)
      end do
    end do
  end do
  call trace_end
END FUNCTION single_compressive

!> =============================================================================
!> Initialize a single 3-D wavenumber, with random amplitude and phase
!> =============================================================================
FUNCTION single_advection (self, m) RESULT (ff)
  class(initial_t)        :: self
  type(mesh_t)          :: m(3)
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,8) :: ff
  !.............................................................................
  type(random_t)        :: random
  integer               :: ix, iy, iz
  real(8)               :: kx, ky, kz
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:single_compressive')
  if (first_time) call self%init
  print*,'adv:',self%k,m%p
  self%phi = 0.0
  do iz=m(3)%lb,m(3)%ub
    kz   =  pi2*self%k(3)*(m(3)%r(iz)+m(3)%p)/m(3)%b + self%phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky   =  pi2*self%k(2)*(m(2)%r(iy)+m(2)%p)/m(2)%b + self%phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*self%k(1)*(m(1)%r(ix)+m(1)%p)/m(1)%b + self%phi(1)
        ff(ix,iy,iz,1) = self%d0*(1.+self%a0(1)*sin(kx)+self%a0(2)*sin(ky)+self%a0(3)*sin(kz))
        ff(ix,iy,iz,2) = self%s0
        ff(ix,iy,iz,3) = self%d0*self%u0(1)
        ff(ix,iy,iz,4) = self%d0*self%u0(2)
        ff(ix,iy,iz,5) = self%d0*self%u0(3)
        ff(ix,iy,iz,6) = self%b0(3)*sin(kz)
        ff(ix,iy,iz,7) = self%b0(1)*sin(kx)
        ff(ix,iy,iz,8) = self%b0(2)*sin(ky)
      end do
    end do
  end do
  print*,'MK',maxval(abs(ff(:,:,:,3:5)))
  call trace_end
END FUNCTION single_advection

END MODULE
