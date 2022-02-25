!===============================================================================
!===============================================================================
MODULE deriv_mod
  USE io_mod
  USE derivatives_mod
  USE central_mod
  USE upwind_mod
  implicit none
  type, public:: deriv_t
  contains
    procedure, nopass:: init
  end type
  private
  type(upwind_t), pointer:: upwind
  type(central_t), pointer:: central
  class(derivatives_t), public, pointer:: d
  type(deriv_t), public:: deriv
  logical:: first_time = .true.
CONTAINS

!===============================================================================
SUBROUTINE init (ds,n)
  real(8), optional:: ds(3)
  integer, optional:: n(3)
  real(8):: ds1(3)
  integer:: n1(3)
  logical:: do_upwind = .false.
  namelist /deriv_params/ do_upwind
  integer:: iostat
  !.............................................................................
  if (.not.first_time) return
  first_time = .false.
  if (present(ds)) then
    ds1 = ds
  else
    ds1 = 1.0
  end if
  if (present(n)) then
    n1 = n
  else
    n1 = 1
  end if
  rewind (io%input)
  read (io%input, deriv_params,iostat=iostat)
  if (io%master) write (*, deriv_params)
  !-----------------------------------------------------------------------------
  if (do_upwind) then
    allocate (upwind)
    d => upwind
  else
    allocate (central)
    d => central
  end if
  d%do_upwind = do_upwind
  d%ds = ds1
  d%n  = n1
END SUBROUTINE init
!===============================================================================
END MODULE deriv_mod
