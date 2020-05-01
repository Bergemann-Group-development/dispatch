!===============================================================================
!> A dummy `initial_mod` to override the default one.
!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE mesh_mod
  USE index_mod
  USE patch_mod
  implicit none
  private

  type, public:: initial_t
    character(len=64):: solver
  contains
    procedure:: init
    procedure:: condition
  end type

  real(8), parameter:: pi = asin(1.0) * 2.0
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self, solver)
  class(initial_t):: self
  character(len=64):: solver
  !----------------------------------------------------------------------------
  self%solver = solver
END SUBROUTINE init

!===============================================================================
!> Initialise a Rayleigh-Taylor experiment.
!> Initial conditions are given in Abel (2011) or Hopkins (2015).
!>
!> Namelist parameters:
!>
!===============================================================================
SUBROUTINE condition (self, m, f, idx)
  class(initial_t):: self
  type(index_t):: idx
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: f
  class(mesh_t), pointer:: m(:), m1, m2, m3
  integer:: i1, i2, i3
  real(4):: y0=0.5, d1=1.0, d2=2.0, deltay=0.025, deltavy=0.025, dymin=0.3, dymax=0.7, y_grav_acc=-0.5, gamma=1.4
  real(4):: d0, p0
  real(8):: x, y, z
  logical, save:: first_time=.true.
  namelist / experiment_params / y0, d1, d2, deltay, deltavy, dymin, dymax, y_grav_acc, gamma
  !----------------------------------------------------------------------------
  call trace%begin('initial_t%init')

  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*, experiment_params)
  end if
  call io%lock%unset

  m1 => m(1)
  m2 => m(2)
  m3 => m(3)
  ! initialise field variables
  do i3=m3%lb,m3%ub
    z = m3%centre_nat + m3%r(i3)
    do i2=m2%lb,m2%ub
      y = m2%centre_nat + m2%r(i2)
      do i1=m1%lb,m1%ub
        x = m1%centre_nat + m1%r(i1)

        ! initial density
        d0 = d1 + (d2 - d1) / (1.0 + exp((y0 - y)/deltay))
        f(i1,i2,i3,idx%d) = d0
        ! add y-velocity perturbation
        if (y >= dymin .and. y <= dymax) then
          f(i1,i2,i3,idx%py) = deltavy * (1.0 + cos(8.0 * pi * (x + 0.25))) &
                                           * (1.0 + cos(5.0 * pi * (y - 0.50))) & 
                                     * d0
        else
          f(i1,i2,i3,idx%py) = 0.0
        end if
        f(i1,i2,i3,idx%px) = 0.0
        f(i1,i2,i3,idx%pz) = 0.0

        ! assign energy/entropy/whatever
        ! the sound speed is set to 1.0 in the upper region
        p0 = d2 / gamma + y_grav_acc * d0 * (y - y0)
        if (trim(self%solver) == 'stagger2_mhd_patch') then
          f(i1,i2,i3,idx%s) = (log(p0) - log(d0) * gamma) &
                                * d0 / (gamma - 1.0)
        else if (trim(self%solver) == 'stagger2e_mhd_patch') then
          f(i1,i2,i3,idx%e) = p0 / (gamma - 1.0)
        else if (trim(self%solver) == 'ramses_mhd_patch') then
          f(i1,i2,i3,idx%e) = p0 + 0.5 * f(i1,i2,i3,idx%py)**2 / d0
        else if (trim(self%solver) == 'zeus_mhd_patch') then
          f(i1,i2,i3,idx%e) = p0 / (gamma - 1.0)
        end if

        ! this experiment is hydro only.
        if (size(f,4) > idx%bx) then
          f(i1,i2,i3,idx%bx) = 0.0
          f(i1,i2,i3,idx%by) = 0.0
          f(i1,i2,i3,idx%bz) = 0.0
        end if
      end do
    end do
  end do

  call trace%end()
END SUBROUTINE condition

END MODULE initial_mod
