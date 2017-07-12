module poisson_commons 
  use amr_commons
  use poisson_parameters

  real(dp),allocatable,dimension(:),  target :: phi_gs, phi_gas         ! Potential
  real(dp),allocatable,dimension(:)          :: phi_old, phi_prv        ! Potential
  real(dp),allocatable,dimension(:,:),target :: f_gs, f_gas             ! 3-force

  real(dp),allocatable,dimension(:), target  :: rho_gas, rho_gs         ! Density
  real(dp),allocatable,dimension(:), target  :: rho_top_gas, rho_top_gs ! Density at last CIC level                                 
  ! Multigrid lookup table for amr -> mg index mapping
  integer, allocatable, dimension(:) :: lookup_mg   ! Lookup table

  ! Communicator arrays for multigrid levels
  type(communicator), allocatable, dimension(:,:) :: active_mg
  type(communicator), allocatable, dimension(:,:) :: emission_mg

  ! Minimum MG level
  integer :: levelmin_mg

  ! Multigrid safety switch
  logical, allocatable, dimension(:) :: safe_mode

  ! Multipole coefficients
  real(dp),dimension(1:ndim+1), target ::multipole_gas, multipole_gs

  real(dp),dimension(:),   pointer :: phi, rho, rho_top, multipole
  real(dp),dimension(:,:), pointer :: f

  ! N-body switch
  logical :: do_nbody

  ! MG statistics
  integer(kind=8), save:: iter_mg=0, n_mg=0
end module poisson_commons
