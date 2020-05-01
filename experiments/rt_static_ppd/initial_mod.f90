!===============================================================================
!> Radiative transfer static protoplanetary disk benchmark
!> http://adsabs.harvard.edu/abs/2015A&A...574A..81R
!> See also: Kuiper & Klessen, 2013, A&A, 555, A7;
!>           Pinte et al., 2009, A&A, 498, 967;
!>           Pascucci et al., 2004, A&A, 417, 793.
!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE units_mod
  USE scaling_mod
  USE mesh_mod
  USE index_mod
  implicit none
  private
  type, public:: initial_t
    logical:: mhd=.false.
    real:: gamma=1.4
    real:: temp0=14.7, temp_exp=-1.0/3.0
    real:: rho0, r_max, z_scale, r_scale, rho_floor
    procedure(initial_conditions), pointer:: condition => initial_conditions
    procedure(disk_density_interface), pointer:: disk_density => pascucci_disk
    class(mesh_t), pointer:: mesh(:)
  contains
    procedure :: init
    procedure :: initial_conditions
  end type
  abstract interface
    FUNCTION disk_density_interface (self, z, r, phi) RESULT (rho)
      import initial_t
      class(initial_t):: self
      real, intent(in):: z, r, phi
      real:: rho
    END FUNCTION disk_density_interface
  end interface

  logical, save:: first_time = .true.
  real, parameter:: tiny = 1.0e-37

CONTAINS
!===============================================================================
!> Read namelist and perform simple initialisations.
!===============================================================================
SUBROUTINE init (self, solver, gamma)
  class(initial_t):: self
  character(len=*):: solver
  real, optional:: gamma
  character(len=8):: disk_type='pascucci'
  real, save:: rho0=8.321e-21, z_scale=0.125, r_scale=0.5, &
               rho_floor=1.0e-36 ! the floor value is determined by the available precision
  namelist /disk_params/ disk_type, rho0, z_scale, r_scale, rho_floor
  !----------------------------------------------------------------------------
  call trace_begin('initial_t%init')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, disk_params)
    if (io%master) write (*, disk_params)
  end if
  call io%lock%unset

  self%rho0 = rho0
  self%z_scale = z_scale
  self%r_scale = r_scale
  self%rho_floor = rho_floor

  if (trim(disk_type) == 'pascucci') then
    ! set up a Pascucci et al. disk
    self%disk_density => pascucci_disk
  else if (trim(disk_type) == 'pinte') then
    ! set up a Pinte et al.-type disk
    self%disk_density => pinte_disk
  else
    call mpi%abort('Invalid disk type selected. Abort!')
  end if

  self%r_max = self%mesh(2)%p + self%mesh(2)%r(self%mesh(2)%ui)
  ! Convert from AU (implied) to cm (needed).
  self%r_max = self%r_max * cgs%AU
  ! The scaling parameters are expressed in terms of `r_max`.
  self%z_scale = self%z_scale * self%r_max
  self%r_scale = self%r_scale * self%r_max

  call trace_end
END SUBROUTINE init

!===============================================================================
!> Set the initial conditions
!===============================================================================
SUBROUTINE initial_conditions (self, mesh, mem, idx)
  class(initial_t) :: self
  type(index_t):: idx
  class(mesh_t), pointer :: mesh(:)
  real, dimension(:,:,:,:), pointer :: mem
  class(mesh_t), pointer :: m1 => null(), m2 => null(), m3 => null()
  real, dimension(:,:,:,:), pointer :: q => null()
  integer :: i1, i2, i3
  real :: z, r, phi, rideal, t0, p0
  !----------------------------------------------------------------------------
  call trace_begin('initial_t::initial_conditions')

  m1 => mesh(1)
  m2 => mesh(2)
  m3 => mesh(3)
  q => mem(:,:,:,:)

  ! fool-proofing
  if (.not. associated(self%disk_density)) &
    call mpi%abort('You have not provided a function to calculate the density. Abort!')

  rideal = cgs%k_b / scaling%mu / cgs%m_u
  do i3=m3%lb,m3%ub
    phi = m3%p + m3%r(i3)
    do i2=m2%lb,m2%ub
      r = abs(m2%p  + m2%r(i2))
      do i1=m1%lb,m1%ub
        z = m1%p  + m1%r(i1)

        ! density
        q(i1,i2,i3,1) = self%disk_density(z,r,phi)
        q(i1,i2,i3,1) = max(self%rho_floor, q(i1,i2,i3,1))
        ! entropy
        ! assume the temperature is initially locally isothermal
        t0 = self%temp0 * (r / self%r_max)**self%temp_exp
        p0 = q(i1,i2,i3,1) * rideal * t0
        q(i1,i2,i3,2) = (log(p0) - log(q(i1,i2,i3,1)) * self%gamma) &
                      * q(i1,i2,i3,1) / (self%gamma - 1.0)
        
        ! momentum components
        q(i1,i2,i3,3) = 0.0
        q(i1,i2,i3,4) = 0.0
        q(i1,i2,i3,5) = 0.0
        ! magnetic field components
        q(i1,i2,i3,6) = 0.0
        q(i1,i2,i3,7) = 0.0
        q(i1,i2,i3,8) = 0.0
      end do
    end do
  end do

  nullify (m1, m2, m3, q)
  call trace_end
END SUBROUTINE initial_conditions

!===============================================================================
!> Calculate the PPD gas density following Pascucci et al. (2004).
!===============================================================================
FUNCTION pascucci_disk (self, z, r, phi) RESULT (rho)
  class(initial_t):: self
  real, intent(in):: z, r, phi
  real:: rho, h
  !----------------------------------------------------------------------------
  h = self%z_scale * (r / self%r_scale)**1.125
  rho = self%rho0 * ((r + tiny) / self%r_scale)**(-1.0) & 
      * exp(-0.25 * cgs%pi * (z / (h + tiny))**2)

END FUNCTION pascucci_disk

!===============================================================================
!> Calculate the PPD gas density following Pinte et al. (2009).
!===============================================================================
FUNCTION pinte_disk (self, z, r, phi) RESULT (rho)
  class(initial_t):: self
  real, intent(in):: z, r, phi
  real:: rho, h
  !----------------------------------------------------------------------------
  h = self%z_scale * (r / self%r_scale)**1.125
  rho = self%rho0 * ((r + tiny) / self%r_scale)**(-2.625) & 
      * exp(-0.5 * (z / (h + tiny))**2)

END FUNCTION pinte_disk

END MODULE
