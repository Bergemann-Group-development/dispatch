!===============================================================================
!> HD module for the Sedov blast wave exercise. This implements the hydrodynamic
!> equations in conservative form on a staggered mesh, where density and entropy
!> are centered at cell-centers, linear momenta are centered on cell faces (walls),
!> and momentum fluxes are centered on cell edges.
!===============================================================================
MODULE mhd_mod
  USE io_mod
  USE mpi_mod
  USE trace_mod
  USE scalar_mod
  USE vector_mod
  USE tensor_mod
  USE gpatch_mod
  USE mesh_mod
  USE initial_mod
  USE index_mod
  USE timestep_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  type, public:: gauss_t                                                        ! Gaussian
    real:: w                                                                    ! width
    real:: p                                                                    ! power exponent
    real:: o(3)                                                                 ! offset
  end type
  !-----------------------------------------------------------------------------
  type, public, extends(gpatch_t):: mhd_t
    real, dimension(:,:,:),   pointer:: d, dd_dt                                ! density
    real, dimension(:,:,:),   pointer:: s, ds_dt                                ! entropy
    real, dimension(:,:,:,:), pointer:: p, dp_dt                                ! momenta
    real, dimension(:,:,:),   pointer:: log_d => null()                         ! scratch (permanent)
    real, dimension(:,:,:),   pointer:: p_gas => null()                         ! scratch (permanent)
    real, dimension(:,:,:,:), pointer:: vel => null()                           ! scratch (permanent)
    logical:: mhd=.false.
    real:: nu(2)                                                                ! diffusion coeffs
    real:: amplitude_d                                                          ! density amplitude
    real:: amplitude_s                                                          ! entropy amplitude
    type(gauss_t):: gauss                                                       ! Gaussian
    type(scalar_t):: scalar
    type(vector_t):: vector
    type(tensor_t):: tensor
  contains
    procedure:: init                                                            ! initialize
    procedure:: pde                                                             ! time derivs
    procedure:: update                                                          ! update
    procedure:: deallocate                                                      ! deallocate scratch
    procedure:: log_density                                                     ! for refinement
    procedure:: log_pressure                                                    ! for refinement
    procedure:: velocity_magnitude                                              ! for refinement
    procedure:: magnetic_field_magnitude                                        ! for refinement
  end type
  !-----------------------------------------------------------------------------
CONTAINS

!===============================================================================
!> Map variables to the current memory slot and initialize the deriv operators.
!> mhd_t%init is called from timestep_t%init, again called from task_list_t%init.
!===============================================================================
SUBROUTINE init(self)
  class(mhd_t):: self
  real, dimension(:,:,:,:), pointer:: mem
  integer:: ix, iy, iz
  real:: x, y, z, f
  logical, save:: first_time=.true.
  real, save:: nu(2)=[0.5,2.0], gamma=5./3.
  namelist /hd_params/ nu, gamma
  !.............................................................................
  call trace%begin ('hd_t%init')
  call self%gpatch_t%init
  !-----------------------------------------------------------------------------
  ! Read and print the namelist only once
  !-----------------------------------------------------------------------------
  call io%lock%set
  if (first_time) then
    first_time=.false.
    rewind (io%input)
    read (io%input, hd_params)
    if (io%master) write (*,hd_params)
    call timestep%init
    ! fool-proofing
    if (timestep%time_order <= 0) call mpi%abort("The Stagger solvers require `time_order` > 0. Abort!")
  end if
  call io%lock%unset
  !-----------------------------------------------------------------------------
  ! Copy parameter values to private values for all patches (note ", save" above)
  !-----------------------------------------------------------------------------
  self%gamma = gamma
  self%nu = nu
  !-----------------------------------------------------------------------------
  ! Physical initial conditions: Sedov blast or advection
  !-----------------------------------------------------------------------------
  self%d => self%mem(:,:,:,self%idx%d ,self%it,1)
  self%s => self%mem(:,:,:,self%idx%s ,self%it,1)
  self%p => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  !-----------------------------------------------------------------------------
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Evaluates the time derivatives of the hydrodynamic variables
!===============================================================================
SUBROUTINE pde (self)
  class(mhd_t):: self
  real:: ds, void
  real, dimension(:,:,:), pointer:: &
    c_sound, log_d, viscosity, p_artificial, &
    diffusion, entropy_per_mass, dissipation
  real, dimension(:,:,:,:), pointer:: &
    delta_u, d_flux, s_flux
  real, dimension(:,:,:,:,:), pointer:: &
    strain_tensor, viscous_stress, p_flux
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('mhd_t%pde', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Allocate scratch scalars (7 scalar chunks)
  !-----------------------------------------------------------------------------
  associate (mx=>self%mesh(1)%gn, my=>self%mesh(2)%gn, mz=>self%mesh(3)%gn)
  if (.not.associated(self%p_gas)) &
  allocate (self%p_gas          (mx,my,mz))
  if (.not.associated(self%log_d)) &
  allocate (self%log_d          (mx,my,mz))
  allocate (c_sound             (mx,my,mz))
  allocate (entropy_per_mass    (mx,my,mz))
  allocate (diffusion           (mx,my,mz))
  allocate (p_artificial        (mx,my,mz))
  allocate (dissipation         (mx,my,mz))
  !-----------------------------------------------------------------------------
  ! Allocate scratch vectors (15 scalar chunks)
  !-----------------------------------------------------------------------------
  if (.not.associated(self%vel)) &
  allocate (self%vel            (mx,my,mz,3))
  allocate (d_flux              (mx,my,mz,3))
  allocate (s_flux              (mx,my,mz,3))
  !-----------------------------------------------------------------------------
  ! Allocate scratch tensors (27 scalar chunks)
  !-----------------------------------------------------------------------------
  allocate (strain_tensor       (mx,my,mz,3,3))
  allocate (viscous_stress      (mx,my,mz,3,3))
  allocate (p_flux              (mx,my,mz,3,3))
  end associate
  !-----------------------------------------------------------------------------
  ! Memory mapping (changes with time steps)
  !-----------------------------------------------------------------------------
  self%d             => self%mem(:,:,:,self%idx%d ,self%it,1)
  self%s             => self%mem(:,:,:,self%idx%s ,self%it,1)
  self%dd_dt         => self%mem(:,:,:,self%idx%d ,self%it,2)
  self%ds_dt         => self%mem(:,:,:,self%idx%s ,self%it,2)
  self%p             => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  self%dp_dt         => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,2)
  !-----------------------------------------------------------------------------
  ! Mesh pointers in operators
  !-----------------------------------------------------------------------------
  self%scalar%mesh        => self%mesh
  self%vector%mesh        => self%mesh
  self%tensor%mesh        => self%mesh
  self%vector%scalar%mesh => self%mesh
  !-----------------------------------------------------------------------------
  ! Association of names, for readability
  !-----------------------------------------------------------------------------
  associate (d=>self%d, s=>self%s, p=>self%p, dd_dt=>self%dd_dt, ds_dt=>self%ds_dt, &
             dp_dt=>self%dp_dt, gamma=>self%gamma, nu=>self%nu, &
             log_d=>self%log_d, p_gas=>self%p_gas, velocity=>self%vel, &
             scalar=>self%scalar, vector=>self%vector, tensor=>self%tensor)
  !-----------------------------------------------------------------------------
  ! Sanity control: is the mass density positive?
  !-----------------------------------------------------------------------------
  if (minval(d) < 0.0) then
    print*,'density negative: id, loc=', self%id, minloc(d)
    call self%output
    call mpi%abort
  end if
  !-----------------------------------------------------------------------------
  ! Equation of state for ideal gas
  !-----------------------------------------------------------------------------
  entropy_per_mass = s/d
  if (gamma==1.0) then
    p_gas        = d
  else
    p_gas        = d**gamma*exp(entropy_per_mass*(gamma-1.0))
  end if
  !-----------------------------------------------------------------------------
  ! Velocity field from linear momenta, log-interpolated to cell faces
  !-----------------------------------------------------------------------------
  log_d          = log(d)
  velocity       = p/exp(vector%face(log_d))
  !-----------------------------------------------------------------------------
  ! Articifical pressure and diffusion added for stability
  !-----------------------------------------------------------------------------
  ds             = maxval(self%mesh%d,self%n > 1)
  p_artificial   = (nu(2)*ds)**2*d*min(0.0,scalar%div(velocity))**2
  c_sound        = sqrt(gamma*((p_gas+p_artificial)/d))
  diffusion      = (nu(1)*ds)*d*c_sound
  !-----------------------------------------------------------------------------
  ! Mass flux, with diffusive term proportional to -grad(ln(d))
  !-----------------------------------------------------------------------------
  d_flux         = p - vector%face(diffusion)*vector%grad(log_d)
  !-----------------------------------------------------------------------------
  ! Entropy flux is mass flux * face-centered entropy per unit mass
  !-----------------------------------------------------------------------------
  s_flux         = d_flux*vector%face(entropy_per_mass) &
                 - vector%face(diffusion)*vector%grad(entropy_per_mass)
  !-----------------------------------------------------------------------------
  ! Momentum flux is the sum of the Reynolds stress d*u*u, the viscous stress,
  ! and a diagonal pressure contribution
  !-----------------------------------------------------------------------------
  strain_tensor  = tensor%strain (velocity)
  viscous_stress = tensor%edge(diffusion) * strain_tensor
  p_flux         = tensor%edge(d,log_d) * tensor%reynolds(velocity) &
                 - viscous_stress
  call tensor%add_diagonal (p_flux, p_gas+p_artificial)
  !-----------------------------------------------------------------------------
  ! Viscous dissipation per unit volume (notice its prop. to strain_tensor**2)
  !-----------------------------------------------------------------------------
  dissipation    = scalar%magnitude (viscous_stress*strain_tensor)
  !-----------------------------------------------------------------------------
  ! Time derivatives
  !-----------------------------------------------------------------------------
  dd_dt          = - scalar%div (d_flux)
  ds_dt          = - scalar%div (s_flux) + dissipation*(d/p_gas)
  dp_dt          = - vector%div (p_flux)
  call scalar%stats (d    , 'd    ')
  call scalar%stats (s    , 's    ')
  call scalar%stats (dd_dt, 'dd_dt')
  call scalar%stats (ds_dt, 'ds_dt')
  call vector%stats (dp_dt, 'dp_dt')
  !-----------------------------------------------------------------------------
  ! Courant condition
  !-----------------------------------------------------------------------------
  self%u_max     = max(scalar%maxval(scalar%magnitude(velocity)+c_sound), &
                    ds*scalar%maxval(abs(dd_dt/d)))
  self%dt        = self%courant*ds/self%u_max
  !-----------------------------------------------------------------------------
  ! Deallocate scratch not needed in refine
  !-----------------------------------------------------------------------------
  deallocate (c_sound, entropy_per_mass, diffusion, p_artificial, dissipation, &
              d_flux, s_flux, strain_tensor, viscous_stress, p_flux)
  !-----------------------------------------------------------------------------
  call trace%end (itimer)
  end associate
END SUBROUTINE pde

!===============================================================================
!===============================================================================
SUBROUTINE update (self)
  class(mhd_t):: self
  !-----------------------------------------------------------------------------
  ! The call to output below catches the state after both guard zone downloads
  ! and boundary conditions.  The call here is a fall back; if another call is
  ! placed at an earlier point in the task handling the call here will have no
  ! effect, since out_next has already been updated.
  !-----------------------------------------------------------------------------
  call self%pde
  call self%courant_condition
  if (io%verbose>2) &
    call self%check_density('before timestep')
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call self%check_density(' after timestep')
END SUBROUTINE update

!===============================================================================
!> Deallocate scratch held in object, after use in refine procedure
!===============================================================================
SUBROUTINE deallocate (self)
  class(mhd_t):: self
 !-----------------------------------------------------------------------------
  !deallocate (self%p_gas, &
  !            self%log_d, &
  !            self%velocity)
END SUBROUTINE deallocate

!===============================================================================
!> Compute the log density and other similar things -- for use in refinement
!===============================================================================
SUBROUTINE log_density (self, out)
  class(mhd_t):: self
  real, dimension(:,:,:):: out
  !...............................................................................
  call trace%begin ('mhd_t%log_density')
  out = log(self%d)
  call trace%end()
END SUBROUTINE log_density
!===============================================================================
SUBROUTINE log_pressure (self, out)
  class(mhd_t):: self
  real, dimension(:,:,:):: out
  !...............................................................................
  out = log(self%p_gas)
END SUBROUTINE log_pressure
!===============================================================================
SUBROUTINE velocity_magnitude (self, out)
  class(mhd_t):: self
  real, dimension(:,:,:):: out
  !...............................................................................
  out = self%scalar%magnitude(self%vel)
END SUBROUTINE velocity_magnitude
!===============================================================================
SUBROUTINE magnetic_field_magnitude (self, out)
  class(mhd_t):: self
  real, dimension(:,:,:):: out
  !...............................................................................
  out = 0.0
END SUBROUTINE magnetic_field_magnitude
!===============================================================================
END MODULE mhd_mod
