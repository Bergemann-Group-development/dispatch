!===============================================================================
!> Experiment module for static RT disk benchmarks (following
!> Kuiper & Klessen, 2013, A&A, 555m A7).
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE mhd_mod
  USE initial_mod
  USE timestep_mod
  USE eos_mod
  implicit none
  private
  type, public, extends(mhd_t):: experiment_t
  contains
    procedure init
    procedure update
  end type
  logical, save:: first_time=.true.
CONTAINS

!===============================================================================
!> Set up the experiment.
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  real, pointer:: ff(:,:,:,:)
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init')
  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  self%nv = 8
  self%mhd = .false.
  call self%mhd_t%init()
  call self%set(bits%static)
  !----------------------------------------------------------------------------
  ! Initial conditions
  !----------------------------------------------------------------------------
  call self%initial%init (self%kind, eos%gamma)
  ff => self%mem(:,:,:,:,self%it,1)
  call self%initial%initial_conditions(self%mesh, ff)
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Update procedure for the experiment.
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  type(timestep_t):: timestep
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')

  call self%mhd_t%update
  call self%courant_condition
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)

  call trace_end
END SUBROUTINE update

END MODULE
