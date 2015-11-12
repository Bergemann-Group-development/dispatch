!> =============================================================================
!> =============================================================================
MODULE experiment_mod
  USE mhd_mod
  USE dust_mod
  USE planet_mod
  USE list_mod
  implicit none
  type, extends(mhd_t):: experiment_t
  contains
    procedure:: init => init_experiment
  end type
  type(mhd_t), pointer:: mhd
  type(dust_t), pointer:: dust
  type(planet_t), pointer:: planet
PRIVATE
PUBLIC experiment_t
CONTAINS

!> =============================================================================
SUBROUTINE init_experiment (self)
  class(experiment_t):: self
  allocate (mhd, dust, planet)
  call mhd%init
  call dust%init
  call planet%init
END SUBROUTINE init_experiment

END MODULE experiment_mod

!> =============================================================================
!> =============================================================================
PROGRAM test
  USE experiment_mod
  implicit none
  type(experiment_t), pointer:: ex

  allocate (ex)
  call ex%init
END
