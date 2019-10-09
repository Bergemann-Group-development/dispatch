!===============================================================================
!> Read the kappa data from Tomida & Hori 2016, which comes in two tables, with
!> Rosseland and Planck means
!===============================================================================
MODULE opacity_reader_mod
  USE io_unit_mod
  USE trace_mod
  implicit none
  private
  type, public:: opac_reader_t
    integer:: nT=0, nd=0, nw=0
    real, dimension(:),     pointer:: T     => null()
    real, dimension(:),     pointer:: d     => null()
    real, dimension(:,:,:), pointer:: table => null()
  contains
    procedure:: init
  end type
  type (opac_reader_t), public:: opac_reader
CONTAINS

!===============================================================================
!> Read in the data from a specific table file, which we might as well 
!> hardwire, since it is specific
!===============================================================================
SUBROUTINE init (self)
  class (opac_reader_t):: self
  real:: row(3)
  integer:: iT, id, iw
  !-----------------------------------------------------------------------------
  call trace%begin ('opac_reader_t%init (void)')
  call trace%end()
END SUBROUTINE init

END MODULE opacity_reader_mod
