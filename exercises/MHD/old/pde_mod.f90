!*******************************************************************************
! Partial Differential Equations corresponding to Magneto-Hydro-Dynamics
!*******************************************************************************
MODULE pde_mod
  USE mhd_mod
  USE vectors
  USE operators
  type(vector):: U, E, f                        ! velocity, electric field, flux
  real, dimension(:,:,:), pointer::     pg, pm  ! gas and magnetic pressure
  character(len=16):: eqn='mhd'
  real:: eta=0.001                              ! diffusion speed
  real:: nu =0.001                              ! diffusion speed
PUBLIC
  real:: gamma=5./3.                            ! ideal gas adiabatic index
CONTAINS

!*******************************************************************************
SUBROUTINE init_pde
  implicit none
  namelist /pde/ eqn, gamma, eta, nu
!...............................................................................
  call allocate_vector (U)
  call allocate_vector (E)
  call allocate_vector (f)
  call allocate_scalar (pg)
  call allocate_scalar (pm)
  rewind(input); read(input,pde); if(master) write(*,pde)
  eta = eta*sum(m%d**2)
  nu  = nu *sum(m%d**2)
END SUBROUTINE init_pde

!*******************************************************************************
SUBROUTINE pde (v, d)
  implicit none
  type(mhd_t):: v, d
  !-----------------------------------------------------------------------------
  ! Compute velocity U=momentum/density, and electric field -E = U x B
  !-----------------------------------------------------------------------------
  call trace_begin ('pde')
  U%x = v%p%x/v%d;                         call scalar_stats (U%x, 'Ux')        ! U = p/d
  U%y = v%p%y/v%d;                         call scalar_stats (U%y, 'Uy')
  U%z = v%p%z/v%d;                         call scalar_stats (U%z, 'Uz')
  call cross_sub (v%B, U, E);              call scalar_stats (E%x, 'Ex')        ! E = - U x B
  !-----------------------------------------------------------------------------
  ! Conservation of magnetic flux
  !-----------------------------------------------------------------------------
  call curl_minus (E, d%B);                call scalar_stats (d%B%z, 'Bz')      ! dB/dt = -curl(E)
  call resistivity (v, d)                                                       ! add resistive effects
  call boundary_vector (d%B)
  if (eqn=='induction') return                                                   ! bail out if enough
  !-----------------------------------------------------------------------------
  ! Conservation of mass
  !-----------------------------------------------------------------------------
  call div_minus (v%p, d%d);               call scalar_stats (d%d, 'dd')        ! d(rho)/dt = -div(p) = -div(d*U)
  !-----------------------------------------------------------------------------
  ! Conservation of entropy
  !-----------------------------------------------------------------------------
  f%x = v%s*U%x
  f%y = v%s*U%y
  f%z = v%s*U%z
  call div_minus (f, d%s);                 call scalar_stats (d%s, 'ds')        ! d(s)/dt = -div(s*U)
  !-----------------------------------------------------------------------------
  ! Conservation of momentum
  !-----------------------------------------------------------------------------
  pg =  v%d**gamma*exp(v%s/v%d*(gamma-1))                                       ! ideal gas
  pm =  0.5*(v%B%x**2+v%B%y**2+v%B%z**2)
  f%x = v%p%x*U%x - v%B%x*v%B%x + pg + pm
  f%y = v%p%x*U%y - v%B%x*v%B%y
  f%z = v%p%x*U%z - v%B%x*v%B%z
  call div_minus (f, d%p%x)                                                     ! d(p%x)/dt = -div(p*U_x + P x_hat)
  f%x = v%p%y*U%x - v%B%y*v%B%x
  f%y = v%p%y*U%y - v%B%y*v%B%y + pg + pm
  f%z = v%p%y*U%z - v%B%y*v%B%z
  call div_minus (f, d%p%y)                                                     ! d(p%y)/dt = -div(p*U_y + P y_hat)
  f%x = v%p%z*U%x - v%B%z*v%B%x
  f%y = v%p%z*U%y - v%B%z*v%B%y
  f%z = v%p%z*U%z - v%B%z*v%B%z + pg + pm
  call div_minus (f, d%p%z)                                                     ! d(p%z)/dt = -div(p*U_z + P z_hat)
  call viscosity (v, d)                                                         ! add viscous effects
  !-----------------------------------------------------------------------------
  ! Ghost zone fill
  !-----------------------------------------------------------------------------
  call boundary        (d%d)
  call boundary        (d%s)
  call boundary_vector (d%p)
  !
  call trace_end
END SUBROUTINE pde

!*******************************************************************************
! Simple dissipative effects
!*******************************************************************************
SUBROUTINE resistivity (v, d)
  implicit none
  type(mhd_t):: v, d
  d%B%x = d%B%x + eta*laplace(v%B%x)
  d%B%y = d%B%y + eta*laplace(v%B%y)
  d%B%z = d%B%z + eta*laplace(v%B%z)
END SUBROUTINE resistivity
SUBROUTINE viscosity (v, d)
  implicit none
  type(mhd_t):: v, d
  d%p%x = d%p%x + nu*laplace(v%p%x)
  d%p%y = d%p%y + nu*laplace(v%p%y)
  d%p%z = d%p%z + nu*laplace(v%p%z)
END SUBROUTINE viscosity

!*******************************************************************************
SUBROUTINE courant (v, c, dt)
  USE io, only: do_trace
  implicit none
  real(8):: c, dt
  type(mhd_t):: v
!-------------------------------------------------------------------------------
! Evaluate Courant condition
!-------------------------------------------------------------------------------
  call trace_begin ('courant')
  if (eqn=='induction') then
    dt = c*minval(m%d)/maxval(sqrt(dot(U,U)))                                  ! velocity magnitue
  else
    dt = c*minval(m%d)/maxval(sqrt(dot(U,U)) + sqrt((2.*pm+gamma*pg)/v%d))     ! velocity magnitue + fast mode speed
  end if
  if (do_trace) print*,'dt =',dt
  call trace_end
END SUBROUTINE courant
END MODULE pde_mod
