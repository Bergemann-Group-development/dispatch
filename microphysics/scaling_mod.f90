!===============================================================================
!> Define code units, in terms of (here) CGS units.  Demonstrates how to pass
!> on the three standard scaling values (l,t,m), and benefit from the standard
!> metadata produced and passed on to Python by self%set (...)
!===============================================================================
MODULE scaling_mod
  USE units_mod
  USE io_mod
  USE trace_mod
  USE standard_scaling_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Extend the units_t data typem which handles scaling metadata output
  !-----------------------------------------------------------------------------
  type, public, extends(standard_scaling_t):: scaling_t
  contains
    procedure:: init
  end type
  type(scaling_t), public:: scaling
CONTAINS

!===============================================================================
!> Default scaling module, optionally overloaded with local version.
!===============================================================================
SUBROUTINE init (self)
  class(scaling_t):: self
  real(8), save:: l=1d0, d=1d0, v=1d0
  character(4), save:: system='CGS'
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /scaling_params/ system, l, v, d
  !----------------------------------------------------------------------------
  call trace_begin('scaling_t%init')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, scaling_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, scaling_params)
    end if
    call io%namelist_warning ('scaling_params', iostat)
    write (io%namelist, scaling_params)
    call self%units_t%init ()
    call self%set (system, l=l, t=l/v, m=d*l**3)
  end if
  call io%lock%unset
END SUBROUTINE init

END MODULE scaling_mod
