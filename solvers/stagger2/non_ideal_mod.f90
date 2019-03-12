!===============================================================================
!> non-ideal MHD module
!===============================================================================
MODULE non_ideal_mod
  USE io_mod
  USE trace_mod
  USE io_unit_mod
  USE scalar_mod
  USE stagger_mod
  implicit none
  private
  type, public:: non_ideal_t
    real:: eta_Ohm, gamma_AD, rho_ions
    logical:: mhd_Ohm, mhd_AD, is_used
  contains
    procedure:: init
    procedure:: update
  end type
  type(non_ideal_t), public:: non_ideal
CONTAINS

!===============================================================================
!> Read parameters
!===============================================================================
SUBROUTINE init (self)
  class (non_ideal_t):: self
  integer:: iostat
  real:: gamma_AD=75., rho_ions=1.0, eta_Ohm=0.0
  logical:: mhd_AD=.false., mhd_Ohm=.false.
  namelist /non_ideal_params/mhd_AD, mhd_Ohm, gamma_AD, eta_Ohm, rho_ions
  !---------------------------------------------------------------------------
  rewind (io_unit%input)
  read (io_unit%input, non_ideal_params,iostat=iostat)
  if (io%master) write (*, non_ideal_params)
  self%mhd_Ohm  = mhd_Ohm
  self%eta_Ohm  = eta_Ohm
  self%mhd_AD   = mhd_AD
  self%gamma_AD = gamma_AD
  self%rho_ions = rho_ions
  self%is_used  = mhd_Ohm .or. mhd_AD
END SUBROUTINE init

!===============================================================================
!> Update heating Q, electric field E, and Courant velocity u_max
!===============================================================================
SUBROUTINE update (self, gn, ds, d, Jx, Jy, Jz, Bx, By, Bz, Q, E, u_max)
  class (non_ideal_t):: self
  integer:: gn(3)
  real(8):: ds
  real:: u_max, gamma_AD
  logical::mhd_AD, mhd_Ohm
  real, dimension(:,:,:,:):: E
  real, dimension(:,:,:):: d, Jx, Jy, Jz, Bx, By, Bz, Q
  real, dimension(:,:,:), allocatable:: Ud_x, Ud_y, Ud_z
  real, dimension(:,:,:), allocatable:: ambi
  integer, save:: itimer=0
  !---------------------------------------------------------------------------
  ! Below is the actual ambipolar drift parameter as given in equation 20 in
  ! Padoan, Zweibel and Nordlund (2000).
  ! To only write ambi, the definition in Padoan et al. is already multiplied
  ! by density**(-3/2) as it appears in the calculation of the drift velocity.
  ! The term d**1.5 is a cosnequence of the single fluid approximation and the
  ! relation n_i = K(n_n/10^5cm^-3)^k + K' (n_n/10^3 cm^-3)^-2
  ! K = 3*10^-3cm^-3, k=0.5 and K'=4.64*10^-4cm^-3.
  ! based on assuming that dust grains do NOT affect ionization balance significantly
  ! => second term falls off quickly => n_i~n_n^0.5
  ! In order to compute this parameter correctly, we need to know TEMPERATURE,
  ! RESOLUTION (not directly applicable here because we use AMR), BOX LENGTH,
  ! AVERAGE DENSITY, COSMIC-RAY IONISATION (K_i parameter:
  ! K_i = n_i*n_H**(-0.5) = x_i*n_H**(-0.5) ) and AVERAGE VOLUMETRIC FLOW RATE av_sigv:
  !---------------------------------------------------------------------------
  ! ambi = 0.28 * (0.1*T)**0.5 * (N/128d0) * (L/10d0)**(-1) * (av_n/200d0)**(-0.5)
  !      * (K_i/1e-5)**(-1) * (av_sigv/2e-9)**(-1) * d**(-1.5)
  !
  ! Alternative definition in Duffin and Pudritz (2008) assumes constant parameter
  ! for gamma_AD
  ! gamma_AD = av_sigv / (m_i + m_n) = 3.28e13 g^-1cm^3s^-1
  ! ambi = 1.4 / (gamma_AD*d**1.5)
  ! rho_i = 1.
  !---------------------------------------------------------------------------
  if (.not.self%is_used) return
  call trace%begin ('non_ideal_t%update', itimer=itimer)
  if (self%mhd_AD) then
    call allocate_scalars_a (gn, ambi, Ud_x, Ud_y, Ud_z)
    ambi = 1./(self%gamma_AD*d*self%rho_ions)
    Ud_x = xdn(ambi) * (zup(Jy)*xdn(zup(Bz)) - yup(Jz)*xdn(yup(By)))
    Ud_y = ydn(ambi) * (xup(Jz)*ydn(xup(Bx)) - zup(Jx)*ydn(zup(Bz)))
    Ud_z = zdn(ambi) * (yup(Jx)*zdn(yup(By)) - xup(Jy)*zdn(xup(Bx)))
    E(:,:,:,1) = E(:,:,:,1) - (Ud_y*Bz - Ud_z*By)
    E(:,:,:,2) = E(:,:,:,2) - (Ud_z*Bx - Ud_x*Bz)
    E(:,:,:,3) = E(:,:,:,3) - (Ud_x*By - Ud_y*Bx)
    Q = Q + self%rho_ions*d*self%gamma_AD*(xup(Ud_x)**2+yup(Ud_y)**2+zup(Ud_z)**2)
    call deallocate_scalars_a (ambi, Ud_x, Ud_y, Ud_z)
    u_max = max(u_max,maxval((xup(Bx)**2+yup(By)**2+zup(Bz)**2)/(self%gamma_AD*self%rho_ions*d*ds)))
  end if
  if (self%mhd_Ohm) then
    E(:,:,:,1) = E(:,:,:,1) + self%eta_Ohm*Jx
    E(:,:,:,2) = E(:,:,:,2) + self%eta_Ohm*Jy
    E(:,:,:,3) = E(:,:,:,3) + self%eta_Ohm*Jz
    Q = Q + yup(zup(self%eta_Ohm*Jx**2)) + zup(xup(self%eta_Ohm*Jy**2)) + xup(yup(self%eta_Ohm*Jz**2))
    u_max = max(u_max,self%eta_Ohm/ds)
  end if
  call trace%end (itimer)
END SUBROUTINE update

END MODULE non_ideal_mod
