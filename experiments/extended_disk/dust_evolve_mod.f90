!===============================================================================
!> Port and test dust evolution from Podolak_1.py
!===============================================================================
MODULE dust_evolve_mod
  USE omp_lib
  USE units_mod
  USE math_mod
  USE trace_mod
  USE io_mod
  USE scaling_mod
  USE particles_mod
  USE omp_timer_mod
  implicit none
  private
  type, public:: dust_evolve_t
    real(KIND=8):: start_time=0d0                           ! start time
    logical:: on=.false.
  contains
    procedure:: init
    procedure:: update
  end type
  type(dust_evolve_t), public:: dust_evolve
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self, ln_bin, lm_cgs, rho_s)
  class(dust_evolve_t) :: self
  integer              :: ln_bin               ! number of bins
  real                 :: lm_cgs(:)            ! particle masses
  real                 :: rho_s               ! particle solid density
  !.............................................................................
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE update (self, dt_dispatch, vrms_in, dd_in, vd_in, dd_out, vd_out)
  class(dust_evolve_t):: self
  real(kind=8), intent(in)  :: dt_dispatch, vrms_in(:)
  real,         intent(in)  :: dd_in(:), vd_in(:,:)
  real,         intent(out) :: dd_out(:), vd_out(:,:)
  dd_out = 0.0
  vd_out = 0.0
END SUBROUTINE update

END MODULE dust_evolve_mod
