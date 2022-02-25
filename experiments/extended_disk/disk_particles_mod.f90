MODULE disk_particles_mod
  USE io_mod
  USE scaling_mod
  USE task_mod
  USE particle_solver_mod
  USE particles_mod
  USE gpatch_mod
  USE disk_model_mod
  USE vector_ops
  USE link_mod
  USE math_mod
  implicit none
  private
  type, public, extends(particle_solver_t):: disk_particles_t
    real :: m_star
    real :: omega
    real :: gm
    logical :: ICs_done = .false.
    real(8), pointer, contiguous:: v_patch(:,:,:,:) => null()
  contains
    procedure:: init
    procedure:: initial_condition
    procedure:: set_attributes
    procedure:: particle_force
    !procedure:: import
    procedure:: refine
  end type
CONTAINS

!===============================================================================
!> Initialize -- a chance to add parameters, or just pass the call on
!===============================================================================
SUBROUTINE init (self)
  class(disk_particles_t):: self
  !.............................................................................
  if (.not.self%on) return
  call self%set_attributes
  call self%particle_solver_t%init
END SUBROUTINE init

!===============================================================================
!> Set attributes by call rather then accessing the attributes directly
!===============================================================================
SUBROUTINE set_attributes (self)
  class(disk_particles_t):: self
  !.............................................................................
  self%m_star = disk_model%m_star
  !FIXME: self%omega = self%patch%omega(3)
  self%gm = scaling%grav*self%m_star
END SUBROUTINE set_attributes

!===============================================================================
!> Particle initial conditions: power-law mid-plane density profile, Gaussian
!> vertical density profile, Keplerian velocities
!===============================================================================
SUBROUTINE initial_condition (self)
  class(disk_particles_t)     :: self
  !.............................................................................
END SUBROUTINE initial_condition

!===============================================================================
!> External force on a particle
!===============================================================================
FUNCTION particle_force (self, p, vp, grad_phi, v) RESULT (a)
  class(disk_particles_t)                  :: self
  class(gpatch_t), pointer                 :: patch
  real, dimension(3)                       :: p, v, a, grad_phi
  real, dimension(2)                       :: vp
  !.............................................................................
  a = 0.0
END FUNCTION particle_force

!===============================================================================
!> Import particles from nbors export lists
!> Verbatim copy from particle_solver_mod, with coordinate and velocity
!> transformations added
!===============================================================================
SUBROUTINE import (self, count)
  class(disk_particles_t)    :: self
  logical, optional          :: count
  !.............................................................................
END SUBROUTINE import

!===============================================================================
!> Adoption of the refinement criterion of the FLASH simulations presented in
!> Schäfer et al. (2020). In their simulations, the resolution is increased if
!> the number of dust particles in a cell exceeds 10 and decreased if it is less
!> than 1.
!> We consider their simulations spanning 10 au to 100 au in the radial
!> dimension, but assume that the vertical gas density profile is given by
!> rho_g(z) = rho_g(z=0)*exp(-z^2/2H_g^2) rather than their Eq. 6.
!> Under these conditions, the mass of each dust particle (their Eq. 15) is
!> equal to m_d = 1.359e26*dust_to_gas g, where dust_to_gas is the dust-to-gas
!> surface density ratio.
!> Their simulations are two-dimensional and do not include the azimuthal
!> dimension. Each particle therefore represents the entire mass in a ring with
!> the circumference 2*pi*r, where r is the radius. To adopt their refinement
!> criterion for three-dimensional DISPATCH simulations, we thus multiply the
!> above mass threshold by the ratio of azimuthal cell size to 2*pi*r.
!===============================================================================
FUNCTION refine (self)
  class(disk_particles_t)               :: self
  integer                               :: refine
  !.............................................................................
  class(gpatch_t), pointer              :: patch
  class(mesh_t), pointer, dimension(:)  :: mesh
  integer                               :: l(3), u(3), ix, iy, iz
  real(8)                               :: r, mass_thres, mass
  real                                  :: sigma_ratio=1.0
  !-----------------------------------------------------------------------------
  call trace%begin ('disk_particles_t%refine_Flash')
  patch => self%patch
  mesh => self%patch%mesh
  !call self%density (self%dd)
  l = mesh%li
  u = mesh%ui
  do iz = l(3), u(3)
    do iy = l(2), u(2)
      do ix = l(1), u(1)
        r = sqrt((patch%radius+mesh(1)%r(ix))**2 + mesh(2)%r(iy)**2)
        mass_thres = 1.359e26*sigma_ratio*mesh(2)%d / &
          (2d0*math%pi*r*scaling%m)
        mass = sum (self%dd(ix,iy,iz,:))*mesh(1)%d*mesh(2)%d*mesh(3)%d
        if (mass >= 10d0 * mass_thres) then
          patch%frefine(ix,iy,iz) = 1.
        else if (mass < mass_thres) then
          patch%frefine(ix,iy,iz) = -1.
        end if
      end do
    end do
  end do
  refine = 0
  if (maxval(patch%frefine(l(1):u(1),l(2):u(2),l(3):u(3))) == 1.) then
    refine = 1
  else if (any(patch%frefine(l(1):u(1),l(2):u(2),l(3):u(3)) == -1.)) then
    refine = -1
  end if
  call trace%end ()
END FUNCTION refine

END MODULE disk_particles_mod
