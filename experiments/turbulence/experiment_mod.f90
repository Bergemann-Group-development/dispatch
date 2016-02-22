!> =============================================================================
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE scalar_mod
  USE vector_mod
  USE mhd_mod
  USE forcing_mod
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
   real:: B0(3)=[0.0,0.0,1.0], d0=1.0
   logical, save:: first_time=.true.
   real, dimension(:,:,:), pointer:: s
   real, dimension(:,:,:,:), pointer:: v
   namelist /experiment_params/ d0, B0
   !----------------------------------------------------------------------------
   !$omp critical (experiment_cr)
   if (first_time) then
     first_time = .false.
     rewind (io%input)
     read (io%input, experiment_params)
     if (io%master) write (*, experiment_params)
     self%d0 = d0
     self%B0 = B0
   end if
   !$omp end critical (experiment_cr)
   s=>self%mem(:,:,:,  1,1,1); call self%d%scalar_map (s)
   s=>self%mem(:,:,:,  2,1,1); call self%s%scalar_map (s)
   v=>self%mem(:,:,:,3:5,1,1); call self%p%vector_map (v)
   v=>self%mem(:,:,:,6:8,1,1); call self%B%vector_map (v)
   s=>self%mem(:,:,:,  1,1,2); call self%dddt%scalar_map (s)
   s=>self%mem(:,:,:,  2,1,2); call self%dsdt%scalar_map (s)
   v=>self%mem(:,:,:,3:5,1,2); call self%dpdt%vector_map (v)
   v=>self%mem(:,:,:,6:8,1,2); call self%dBdt%vector_map (v)
   self%d = 1.0
   self%p = 0.0
   self%B = self%B0
END SUBROUTINE init

!> =============================================================================
SUBROUTINE update (self)
   class(experiment_t):: self
   !----------------------------------------------------------------------------
   call self%mhd_t%update
   self%dpdt = self%dpdt + self%force%selected (self%time, self%mesh)
END SUBROUTINE update

END MODULE
