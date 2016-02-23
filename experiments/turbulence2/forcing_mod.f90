!> =============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!> =============================================================================
MODULE forcing_mod
  USE io_mod
  USE trace_mod
  USE vector_mod
  USE mesh_mod
  USE random_mod
  implicit none
  private
  real, parameter:: pi2=8.0*atan(1.0)
  type, public:: force_t
    integer:: seed=-11
    integer:: gn(3)
    real:: a(3), phi(3)
    real:: t_turn, t_save=0.0
    real:: ampl(3), k(3)
    procedure(single_solenoidal), pointer:: selected=>null()
  end type
 logical, save:: first_time=.true.
CONTAINS

!> =============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!> =============================================================================
SUBROUTINE init (self)
   class(force_t):: self
   real:: t_turn=0.3, ampl(3)=0.0, k(3)=1.0
   namelist /force_params/ k, ampl, t_turn
   !----------------------------------------------------------------------------
   !$omp critical (experiment_cr)
   if (first_time) then
     first_time = .false.
     rewind (io%input)
     k      = self%k
     ampl   = self%ampl
     t_turn = self%t_turn
     read (io%input, force_params)
     if (io%master) write (*, force_params)
     self%t_turn = t_turn
     self%ampl   = ampl
     self%k      = k
     self%selected => single_solenoidal 
   end if
   !$omp end critical (experiment_cr)
END SUBROUTINE init

!> =============================================================================
!> Forcing on a single 3-D wavenumber, changing amplitude and phase periodically
!> =============================================================================
FUNCTION single_solenoidal (self, time, m) RESULT (ff)
  class(force_t)        :: self
  real(8)               :: time
  real, dimension(self%gn(1),self%gn(2),self%gn(3),3) :: ff
  type(mesh_t)          :: m(3)
  !.............................................................................
  type(random_t)        :: random
  integer               :: ix, iy, iz
  real                  :: kx, ky, kz
  !-----------------------------------------------------------------------------
  if (time >= self%t_save) then
    self%a = self%ampl*random%ran3()
    self%phi = pi2*random%ran3()
    self%t_save = self%t_save + self%t_turn
  end if

  do iz=m(3)%lb,m(3)%ub
    kz = pi2*self%k(3)*(m(3)%r(iz)+m(3)%p) + self%phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky = pi2*self%k(2)*(m(2)%r(iy)+m(2)%p) + self%phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx = pi2*self%k(1)*(m(1)%r(ix)+m(1)%p) + self%phi(1)
        ff(ix,iy,iz,1) = self%a(1)*(sin(ky)+sin(kz))
        ff(ix,iy,iz,2) = self%a(2)*(sin(kz)+sin(kx))
        ff(ix,iy,iz,3) = self%a(3)*(sin(kx)+sin(ky))
      end do
    end do
  end do
END FUNCTION single_solenoidal

END MODULE
