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
  real, parameter:: pi2=8.0*atan(1.0)
  type, public:: initial_t
    integer:: seed
    real:: a(3), b(3), phi(3), d0, s0
    real:: a0(3), k(3), b0(3)
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
   real, save:: a0(3)=0.1, k(3)=1.0, d0=1.0, s0=0.0, b0(3)=0.0
   character(len=32):: type='single_solenoidal'
   namelist /initial_params/ seed, k, a0, b0, d0, s0, type
   !----------------------------------------------------------------------------
   call trace_begin('initial_t%init')
   if (first_time) then
    !$omp critical (experiment_cr)
     if (first_time) then
       first_time = .false.
       rewind (io%input)
       read (io%input, initial_params)
       if (io%master) write (*, initial_params)
     end if
    !$omp end critical (experiment_cr)
   end if
   self%seed   = seed
   self%k      = k
   self%a0     = a0
   self%b0     = b0
   self%d0     = d0
   self%s0     = s0
   select case (type)
   case ('single_solenoidal')
     self%selected => single_solenoidal
   case ('single_compressive')
     self%selected => single_compressive
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
  class(initial_t)        :: self
  type(mesh_t)          :: m(3)
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,8) :: ff
  !.............................................................................
  type(random_t)        :: random
  integer               :: ix, iy, iz
  real                  :: kx, ky, kz
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:single_solenoidal')
  if (first_time) call self%init
  self%a = self%a0*random%ran3()
  self%b = self%b0*random%ran3()
  self%phi = pi2*random%ran3()
  do iz=m(3)%lb,m(3)%ub
    kz   =  pi2*self%k(3)*(m(3)%r(iz)+m(3)%p) + self%phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky   =  pi2*self%k(2)*(m(2)%r(iy)+m(2)%p) + self%phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*self%k(1)*(m(1)%r(ix)+m(1)%p) + self%phi(1)
        ff(ix,iy,iz,1) = self%d0
        ff(ix,iy,iz,2) = self%s0
        ff(ix,iy,iz,3) = self%a(1)*(sin(ky)+sin(kz))
        ff(ix,iy,iz,4) = self%a(2)*(sin(kz)+sin(kx))
        ff(ix,iy,iz,5) = self%a(3)*(sin(kx)+sin(ky))
        ff(ix,iy,iz,6) = self%b(1)*(sin(kx)+sin(ky))
        ff(ix,iy,iz,7) = self%b(2)*(sin(kx)+sin(ky))
        ff(ix,iy,iz,8) = self%b(3)*(sin(kx)+sin(ky))
      end do
    end do
  end do
  print*,'MK1',maxval(abs(ff(:,:,:,3:5)))
  call trace_end
END FUNCTION single_solenoidal

!> =============================================================================
!> Forcing on a single 3-D wavenumber, changing amplitude and phase periodically
!> =============================================================================
FUNCTION single_compressive (self, m) RESULT (ff)
  class(initial_t)        :: self
  type(mesh_t)          :: m(3)
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,8) :: ff
  !.............................................................................
  type(random_t)        :: random
  integer               :: ix, iy, iz
  real                  :: kx, ky, kz
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:single_compressive')
  if (first_time) call self%init
  self%a = self%a0*random%ran3()
  self%b = self%b0*random%ran3()
  self%phi = pi2*random%ran3()
  do iz=m(3)%lb,m(3)%ub
    kz   =  pi2*self%k(3)*(m(3)%r(iz)+m(3)%p) + self%phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky   =  pi2*self%k(2)*(m(2)%r(iy)+m(2)%p) + self%phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*self%k(1)*(m(1)%r(ix)+m(1)%p) + self%phi(1)
        ff(ix,iy,iz,1) = self%d0
        ff(ix,iy,iz,2) = self%s0
        ff(ix,iy,iz,3) = self%a(1)*sin(kx)
        ff(ix,iy,iz,4) = self%a(2)*sin(ky)
        ff(ix,iy,iz,5) = self%a(3)*sin(kz)
        ff(ix,iy,iz,6) = self%b(1)*sin(kx)
        ff(ix,iy,iz,7) = self%b(2)*sin(ky)
        ff(ix,iy,iz,8) = self%b(3)*sin(kz)
      end do
    end do
  end do
  print*,'MK',maxval(abs(ff(:,:,:,3:5)))
  call trace_end
END FUNCTION single_compressive

END MODULE
