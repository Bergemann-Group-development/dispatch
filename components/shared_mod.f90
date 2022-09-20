!===============================================================================
!> This module is a placeholder for shared information, solving the problem of
!> making a piece of information available btw any two points in the code, without
!> having to worry about module dependency.  Since this module does not depend
!> on anything else, we can put shared information here
!===============================================================================
MODULE shared_mod
  implicit none
  public
  type, public:: shared_t
    integer:: levelmin      = 99
    integer:: levelmax      = 0
    integer:: levelroot     = 0
    integer:: n_repeat(3)   = 1
    integer:: n_patch       = 0
    integer:: n_initial     = 0
    integer:: n_add         = 0
    integer:: n_patch_total = 0
    integer:: patch_n(3)    = 0
    integer:: n_queue       = 0
    integer:: mpi_size      = 0
    integer:: mpi_rank      = 0
    integer:: benchmark     = 0
    logical:: mpi_master    = .true.
    logical:: all_src       = .false.
    logical:: refine_on     = .false.
    logical:: debug         = .false.
    real(8):: refine_ratio  = 2d0
    real(8):: tlist_size(3) = 1d0
    real:: f_zoom(3)        = 1.0
    real:: extra_grace      = 0.0
    real(8):: mean_density  = 0d0
    real(8):: sec_per_report= 1d0
  contains
    procedure, nopass:: abort
    procedure, nopass:: option_l
    procedure, nopass:: option_i
    procedure, nopass:: option_i_3d
    procedure, nopass:: option_i8
    procedure, nopass:: option_r
    procedure, nopass:: option_r8
    procedure, nopass:: option_r8_3d
    procedure, nopass:: option_s
    generic:: option => option_l, option_i, option_i8, option_r, option_r8, &
                        option_r8_3d, option_s, option_i_3d
  end type
  type(shared_t):: shared
CONTAINS

SUBROUTINE abort (label)
  character(len=*):: label
  write (*,*) label
  stop 'shared%abort'
END SUBROUTINE abort

!===============================================================================
!> Convenience functions, to handle optional arguments
!===============================================================================
FUNCTION option_l (in, def) RESULT (out)
  logical:: in, def, out
  optional:: in, def
  if (present(in)) then
    out = in
  else if (present(def)) then
    out = def
  else
    out = .false.
  end if
END FUNCTION
!===============================================================================
FUNCTION option_i (in, def) RESULT (out)
  integer:: in, def, out
  optional:: in
  if (present(in)) then
    out = in
  else
    out = def
  end if
END FUNCTION
!===============================================================================
FUNCTION option_i_3d (in, def) RESULT (out)
  integer, dimension(3):: in, def, out
  optional:: in
  if (present(in)) then
    out = in
  else
    out = def
  end if
END FUNCTION
!===============================================================================
FUNCTION option_i8 (in, def) RESULT (out)
  integer(8):: in, def, out
  optional:: in
  if (present(in)) then
    out = in
  else
    out = def
  end if
END FUNCTION
!===============================================================================
FUNCTION option_r (in, def) RESULT (out)
  real(4):: in, def, out
  optional:: in
  if (present(in)) then
    out = in
  else
    out = def
  end if
END FUNCTION
!===============================================================================
FUNCTION option_r8 (in, def) RESULT (out)
  real(8):: in, out, def
  optional:: in
  if (present(in)) then
    out = in
  else
    out = def
  end if
END FUNCTION
!===============================================================================
FUNCTION option_r8_3d (in, def) RESULT (out)
  real(8), dimension(3):: in, out, def
  optional:: in
  if (present(in)) then
    out = in
  else
    out = def
  end if
END FUNCTION
!===============================================================================
FUNCTION option_s (in, def) RESULT (out)
  character(*):: in, def
  character(64):: out
  optional:: in
  if (present(in)) then
    out = in
  else
    out = def
  end if
END FUNCTION

END MODULE shared_mod
