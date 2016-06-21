!*******************************************************************************
! Partial Differential Equations corresponding to Magneto-Hydro-Dynamics
!*******************************************************************************
MODULE mhd_mod
  !.............................................................................
  USE patch_mod
  implicit none
PRIVATE
  type, public, extends(patch_t):: mhd_t
    logical:: mhd=.true.
  contains
    procedure:: init
    procedure:: update
    procedure:: log_density
    procedure:: log_pressure
    procedure:: grav_potential
    procedure:: velocity_magnitude
    procedure:: magnetic_field_magnitude
  end type
  logical:: first_time=.true.
  character(len=16):: eqn='mhd'
  real, public :: eta=0.001                           ! diffusion speed
  real, public :: nu(6)=(/0.05,0.7,3.0,2.0,0.1,1.0/)  ! diffusion speeds
  real(kind=8):: courant_number=0.2
CONTAINS

!*******************************************************************************
SUBROUTINE init (self)
  class(mhd_t):: self
!...............................................................................
  !$omp critical (pde)
  !print'(a,i4,3(2x,6f8.3))','MK',m(1)%id,(m(i)%h(4:9),i=1,3)
  if (first_time) then
    first_time = .false.
  end if
  !$omp end critical (pde)
END SUBROUTINE init

!*******************************************************************************
SUBROUTINE update (self)
  class(mhd_t):: self
  !-----------------------------------------------------------------------------
END SUBROUTINE update

!*******************************************************************************
SUBROUTINE log_density (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:):: v
  !-----------------------------------------------------------------------------
END SUBROUTINE

!*******************************************************************************
SUBROUTINE log_pressure (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:):: v
  !-----------------------------------------------------------------------------
END SUBROUTINE

!*******************************************************************************
SUBROUTINE grav_potential (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:):: v
  !-----------------------------------------------------------------------------
END SUBROUTINE

!*******************************************************************************
SUBROUTINE velocity_magnitude (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:):: v
  !-----------------------------------------------------------------------------
END SUBROUTINE

!*******************************************************************************
SUBROUTINE magnetic_field_magnitude (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:):: v
  !-----------------------------------------------------------------------------
END SUBROUTINE

END MODULE
