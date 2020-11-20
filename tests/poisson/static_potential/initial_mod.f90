!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE index_mod
  USE kinds_mod
  implicit none
  private

  type, public:: initial_t
    character(len=64):: solver
    real:: r0
  contains
    procedure:: init
    procedure:: condition
  end type

  real:: global_r0

public func_d0, func_phi0
CONTAINS

!===============================================================================
SUBROUTINE init (self, solver, gamma)
  class(initial_t):: self
  character(len=64):: solver
  real, optional:: gamma
  real:: r0 = 1.0
  logical:: first_time = .true.
  namelist / IC_params / r0
  !----------------------------------------------------------------------------
  call trace%begin('initial_t%init')

  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, initial_params)
    if (io%master) write (*, initial_params)
  end if
  call io%lock%unset
  self%r0 = r0
  ! save global copy for use by `func_d0`.
  global_r0 = r0
  self%solver = solver
  
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Initialise a static potential and correspondiong density given by Guillet &
!> Teyssier (2011, JCoPh, 230, p. 4756; eqn. 4). There are no dynamics in this
!> setup, rather it is used to test the accuracy of the Poisson solver.
!>
!> The potential is given by:
!>   Phi(r) = ln( (r/r0)^2 + 1),
!> where r0 controls the concentration of the profile. The density is meanwhile
!> given by:
!>   \rho(r) = 4r0^2 / (r^2 + r0^2)^2.
!===============================================================================
SUBROUTINE condition (self, m, f, idx)
  class(initial_t):: self
  type(index_t):: idx
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: f
  class(mesh_t), pointer:: m(:), m1, m2, m3
  integer:: i1, i2, i3
  real(kind=KindScalarVar):: x, y ,z
  !----------------------------------------------------------------------------
  call trace%begin('initial_t%condition')

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

        f(i1,i2,i3,idx%d) = func_d0(x,y)
        ! intentionally intialise the grav. pot. to a poor estimate of the true value
        f(i1,i2,i3,idx%phi) = 0.0

      end do
    end do
  end do

  call trace%end()
END SUBROUTINE condition

!===============================================================================
!> Function to calculate the exact potential (eq. 4 of Guillet & Teyssier)
!===============================================================================
FUNCTION func_phi0 (x, y) result (phi0)
  real(kind=KindScalarVar), intent(in):: x, y
  real(kind=KindScalarVar):: r, phi0
  !----------------------------------------------------------------------------
  r = sqrt(x**2 + y**2)
  phi0 = log((r / global_r0)**2 + 1.0) + 1.0

END FUNCTION func_phi0

!===============================================================================
!> Function to calculate the density that matches the exact potential
!===============================================================================
FUNCTION func_d0 (x, y) result (d0)
  real(kind=KindScalarVar), intent(in):: x, y
  real(kind=KindScalarVar):: r, d0
  !----------------------------------------------------------------------------
  r = sqrt(x**2 + y**2)
  d0 = 4.0 * global_r0**2 / (r**2 + global_r0**2)**2

END FUNCTION func_d0

END MODULE initial_mod
