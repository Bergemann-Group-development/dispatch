!===============================================================================
!> Interface module to the Bifrost code, with global initialization and
!> finalization (if needed) of the code.
!===============================================================================
MODULE extern_mod
  USE params_mod
  USE parallel_mod
  USE units_mod
  USE stagger_mod
  USE bifrost_timestep_mod
  USE bifrost_io_mod
  USE math_mod
  USE quench_mod
  USE timing_mod
  USE adaption_mod
  implicit none
  private
  type, public:: extern_t
  contains
    procedure:: init
    procedure:: end
  end type
  character(len=80), save:: id = &
  '$Id$'
  type(extern_t), public:: extern
CONTAINS

!===============================================================================
!> This procedure is called once, at the start of the run
!===============================================================================
SUBROUTINE init (self)
  class(extern_t):: self
  integer :: tl_tp, tl_cr, tl_cm
  !-----------------------------------------------------------------------------
  reread = 0
  call init_params                                                      ! reads in the basic simulation parameters. Could and should be replaced with a simple subroutine inside DISPATCH
  call init_parallel
  dispatch = .true.
!
  if (myid.eq.R_MASTER) then
    write (*,*)
    write (*,FMT=ilo) 'main    '
    write (*,*) id
    if (time_lim > 0.0) then
      call system_clock(tl_tp,COUNT_RATE=tl_cr,COUNT_MAX=tl_cm)
    endif
  endif
!
  call init_units
  call init_io
  call init_timestep
  call init_math
  call init_quench
  call init_timing
  call init_adaption
  call init_standard_boundaries
!> The latter should be initialized per patch basis until they will be reworked
!> into globally accessible ones or fully incorporated into other modules.
  ! working on this one:           
  !call init_stagger

  clean=nstepstart 
END SUBROUTINE init



!===============================================================================
!> End of run
!===============================================================================
SUBROUTINE end (self)
  class(extern_t):: self
  !-----------------------------------------------------------------------------
  call finalize_timing
  !***call finalize_io
  call finalize_parallel
END SUBROUTINE end

END MODULE extern_mod
