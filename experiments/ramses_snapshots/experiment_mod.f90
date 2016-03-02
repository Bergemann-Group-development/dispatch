!> =============================================================================
!> Using a RAMSES snapshot, set up a DISPATCH simulation.
!>
!> In this first implementation, a patch is created for every 2x2x2 cells ("grid").
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE mhd_mod
  USE scaling_mod
  implicit none
  private
  type, public, inherits(mhd_t):: experiment_t
    integer:: nout, minlvl, maxlvl,
    real:: xmin, xmax, ymin, ymax, zmin, zmax, rank_load
    character(len=64):: path='.'
  contains
    procedure init
  end type
CONTAINS

!> =============================================================================
!> Facilitate call to `read_ramses_snapshot`.  All heavy lifting is done in that
!> subroutine.
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  logical, save:: first_time=.true
  integer, save:: nout, minlvl, maxlvl, 
  real, save:: xmin, xmax, ymin, ymax, zmin, zmax, rank_load
  character(len=64), save:: path='.'
  namelist / init_params / nout, xmin, xmax, ymin, ymax, zmin, zmax, &
                           minlvl, maxlvl, rank_load, path
  !.............................................................................
  call trace_begin('experiment_t%init')
  !-----------------------------------------------------------------------------
  ! The first thread reads the namelist parameters and save them
  !-----------------------------------------------------------------------------
  !$omp critical (init_cr)
  if (first_time) then
    first_time=.false.
    ! default values
    nout = -1
    xmin = 0.0; xmax = 1.0
    ymin = 0.0; ymax = 1.0
    zmin = 0.0; zmax = 1.0
    minlvl = 1
    maxlvl = 1
    rank_load = 1.0
    path = './'
    ! read namelist
    rewind (io%input); read (io%input,init_params)
    if (io%master) write(*, init_params)
  end if
  !$omp end critical (init_cr)
  !-----------------------------------------------------------------------------
  ! All threads save the parameters in their private parts, including the 
  ! scaling information.  The selection of RAMSES data we read is based on
  ! the parameters defined here.
  !-----------------------------------------------------------------------------
  self%xmin = xmin
  self%xmax = xmax
  self%ymin = ymin
  self%ymax = ymax
  self%zmin = zmin
  self%zmax = zmax
  self%minlvl = minlvl
  self%maxlvl = maxlvl
  self%rank_load = rank_load
  self%path = path
  call scaling%init
  !-----------------------------------------------------------------------------
  ! Each thread reads the RAMSES grids belonging to it, filling in the data,
  ! initializing the patch structure of every grid
  !-----------------------------------------------------------------------------
  call self%read_ramses
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  type(timestep_t):: timestep
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call self%mhd_t%update
  call self%courant_condition
  call timestep%update (self%id, self%iit, self%dt, self%mem)
  call trace_end
END SUBROUTINE update

END MODULE initial_mod
