!> =====================================================================
!> Simple opacity, expressed as
!> kappa = (T/T0)^-3.5*(rho/rho0)^2
!> T0 = 290K, rho0 = 1 g cm-3
!> =====================================================================
MODULE opacities_mod
  USE io_mod
  USE trace_mod
  USE scaling_mod
  implicit none
  private
  type, public ::  opacity_t
    real :: t0    = 290
    real :: rho0 = 1
  contains
    procedure :: init   => init_shadow_opacity
    procedure :: solve  => rt_shadow_opacity
  end type
CONTAINS

SUBROUTINE init_shadow_opacity(self)
  class(opacity_t) :: self
  ! --------------------------------------------------------------------
  call trace_begin('opacity_mod::rt_shadow_opacity',3)
  self%t0   = 1./(290./scaling%temp)
  self%rho0 = 1./(1.  /scaling%d)
  call trace_end
END SUBROUTINE init_shadow_opacity

FUNCTION rt_shadow_opacity (self, rho, temp, nwb) RESULT (kappa)
  class(opacity_t) :: self
  real(4), dimension(:,:,:), intent(in) :: rho, temp
  real(4), dimension(size(rho,1),size(rho,2),size(rho,3), nwb) :: kappa
  integer :: nwb, wl
  !.....................................................................
  call trace_begin('opacity_mod::rt_shadow_opacity',3)
  
  kappa(:,:,:,1) = ((temp(:,:,:)*self%t0)**(-3.5)*(rho(:,:,:)*self%rho0)**2)
  
  if (nwb>1) then
   do wl = 2, nwb
     kappa(:,:,:,wl) = kappa(:,:,:,1)
   end do
  end if
  
  if (io%verbose>2) then
    if (any(isnan(rho)).or.any(isnan(temp))) then
      error stop "rt_shadow_opac%solve: rho or temp has NaN"
    end if
    if (any(isnan(kappa))) then
      print *,'=================== kappa ================'
      print *,kappa
      error stop "rt_shadow_opacity%check: kappa has NaN"
    end if
  end if
  
  call trace_end
END FUNCTION rt_shadow_opacity

!> =============================================================================
!> Set the opacity to be constant per unit mass and independent of temperature
!> assuming kappa_grey is 1 gm cm-3
!> =============================================================================
FUNCTION opac_gr (rho, temp, wl) result(kappa)
  real(4), dimension(:,:,:), intent(in) :: rho, temp
  integer :: wl, iwl
  real(4), dimension(size(rho,1),size(rho,2),size(rho,3)) :: rho0
  real(4), dimension(size(rho,1),size(rho,2),size(rho,3), wl) :: kappa
  real(4) :: op_sc, sc_rho
  !.............................................................................
  op_sc = scaling%m/scaling%l**2
  sc_rho = 1_8/scaling%d
  rho0 = rho * sc_rho
  kappa(:,:,:,1) = rho0 * scaling%m/scaling%l**2
  if (wl>1) then
    do iwl =2, wl
      kappa(:,:,:,iwl) = kappa(:,:,:,1)
    end do
  end if
  
  if (io%verbose > 2) then
    if (any(isnan(kappa))) then
      print *,'=================== kappa ================'
      print *,kappa
      error stop "rt_shadow_opacity%check: kappa has NaN"
    end if
  end if
  
END FUNCTION opac_gr

END MODULE opacities_mod
