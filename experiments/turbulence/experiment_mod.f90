!> =============================================================================
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE scalar_mod
  USE vector_mod
  USE mhd_mod
  USE forcing_mod
  USE timestep_mod
  implicit none
  private
  real, parameter:: pi2=8.0*atan(1.0)
  type, public, extends(mhd_t):: experiment_t
    integer:: i_d=1
    real:: d0, B0(3)
    real:: G
    type(force_t):: force
  contains
    procedure init
    procedure update
  end type
CONTAINS

!> =============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  real, save:: B0(3)=[0.0,0.0,1.0], d0=1.0
  logical, save:: first_time=.true.
  real, dimension(:,:,:), pointer:: s
  real, dimension(:,:,:,:), pointer:: v
  integer:: n(3)
  namelist /experiment_params/ d0, B0
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init')
  !$omp critical (experiment_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*, experiment_params)
  end if
  !$omp end critical (experiment_cr)
  self%d0 = d0
  self%B0 = B0
  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  self%nv = 8
  self%nt = 7
  call self%patch_t%init
  !----------------------------------------------------------------------------
  ! Map variables to memory
  !----------------------------------------------------------------------------
  n = self%mesh%gn
  s=>self%mem(:,:,:,  1,1,1); call self%d%scalar_map    (s, n, self%ds)
  s=>self%mem(:,:,:,  2,1,1); call self%s%scalar_map    (s, n, self%ds)
  v=>self%mem(:,:,:,3:5,1,1); call self%p%vector_map    (v, n, self%ds)
  v=>self%mem(:,:,:,6:8,1,1); call self%B%vector_map    (v, n, self%ds)
  s=>self%mem(:,:,:,  1,1,2); call self%dddt%scalar_map (s, n, self%ds)
  s=>self%mem(:,:,:,  2,1,2); call self%dsdt%scalar_map (s, n, self%ds)
  v=>self%mem(:,:,:,3:5,1,2); call self%dpdt%vector_map (v, n, self%ds)
  v=>self%mem(:,:,:,6:8,1,2); call self%dBdt%vector_map (v, n, self%ds)
  !----------------------------------------------------------------------------
  ! Initial conditions
  !----------------------------------------------------------------------------
  self%d = self%d0
  self%p = 0.0
  self%B = self%B0
END SUBROUTINE init

!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  type(timestep_t):: timestep
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call self%mhd_t%update
  !self%dpdt = self%dpdt + self%force%selected (self%time, self%mesh)
  call self%courant_condition
  call timestep%update (self%id, self%iit, self%dt, self%mem)
  call self%patch_t%update
  call trace_end
END SUBROUTINE update

END MODULE
