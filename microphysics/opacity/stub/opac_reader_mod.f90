!> =============================================================================
!> =============================================================================
MODULE opac_reader_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  type, public:: opac_reader_t
    integer:: skip=1
    integer:: nT=0, nd=0
    real(8), dimension(:), pointer:: T, d
    real(8), dimension(:,:,:), pointer:: table
    character(len=64):: dir='data/opacity'
  contains
    procedure:: init
  end type
  type (opac_reader_t), public:: opac_reader
CONTAINS

!> =============================================================================
!> Read in the data from a specific table file, which we might as well 
!> hardwire, since it is specific
!> =============================================================================
SUBROUTINE init (self)
  class (opac_reader_t):: self
END SUBROUTINE init

END MODULE opac_reader_mod
