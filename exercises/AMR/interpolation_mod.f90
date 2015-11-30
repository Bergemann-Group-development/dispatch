!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE interpolation_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  type, public :: interpolation_t
    real   :: factor=1.5
    logical:: do_upwind=.false.
  contains
    procedure:: init
    procedure:: prolongation
    procedure:: restriction
  end type
  logical:: first_time = .true.
CONTAINS
!> =============================================================================
SUBROUTINE init (self, dims, ds)
  class(derivatives_t):: self
  integer:: dims(3)
  real:: ds(3)
  logical:: do_upwind = .false.
  real:: factor = 1.5
  namelist /interpol_params/ do_upwind, factor
  !.............................................................................
  call trace_begin('interpolation_t%init', 1)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, interpol_params)
    if (io%master) write (*, interpol_params)
  end if
  self%do_upwind = do_upwind
  self%factor = factor
  self%n = dims
  self%ds = ds
  call trace_end
END SUBROUTINE init
!> =============================================================================
!> Prolongate data from a coarse patch to a fine patch.
!> init is a switch to indicate we need to prolongate all the patch (initial refinement)
!> guard is a switch that signals updating of the guard cells
!> =============================================================================
SUBROUTINE prolongate (self, coarse, fine, init, guard)
  class(interpolation_t):: self
  class(patch_t):: coarse, fine
  logical :: init, guard
  !.............................................................................
  call trace_begin ('interpolation_t%prolongate', 3)
  if (first_time) then
    first_time = .false.
  end if
  call trace_end
END FUNCTION prolongate
!> =============================================================================
!> Restrict data from a fine patch to a coarse patch.
!> upload is a switch to indicate we need to restrict in interior of the fine patch (upload of data)
!> guard is a switch that signals updating of the guard cells
!> =============================================================================
SUBROUTINE restrict (self, coarse, fine, upload, guard)
  class(interpolation_t):: self
  class(patch_t):: coarse, fine
  logical :: upload, guard
  !.............................................................................
  call trace_begin ('interpolation_t%prolongation', 3)
  if (first_time) then
    first_time = .false.
  end if
  call trace_end
END FUNCTION restrict
!
END MODULE interpolation_mod
