!===============================================================================
!> Data type that holds information about point sources
!===============================================================================
MODULE point_sources_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  !.............................................................................
  type, public:: point_source_t
    real(8):: position(3)
    real(8):: radius
    real:: temperature
  contains
    procedure:: init => init_point_source
  end type
  !.............................................................................
  type, public:: point_sources_t
    integer:: n = 0
    type(point_source_t), pointer:: point_source(:) => null()
  contains
    procedure:: init
  end type
  !.............................................................................
  type(point_sources_t), public:: point_sources
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init_point_source (self)
  class(point_source_t):: self
END SUBROUTINE init_point_source

!===============================================================================
!===============================================================================
SUBROUTINE init (self)
  class(point_sources_t):: self
  integer:: iostat, n
  real(8):: position(3)=0d0, radius=1d0, temperature=300
  namelist /point_source_params/ position, radius, temperature
  !-----------------------------------------------------------------------------
  call io%lock%set
  rewind (io%input)
  n = 0
  do while (iostat==0)
    read (io%input, point_source_params, iostat=iostat)
    n = n + 1
  end do
  self%n  = n
  rewind (io%input)
  allocate (self%point_source(n))
  do while (iostat==0)
    read (io%input, point_source_params, iostat=iostat)
    self%point_source(n)%position = position
    self%point_source(n)%radius = radius
    self%point_source(n)%temperature = temperature
  end do
  call io%lock%unset
END SUBROUTINE

END MODULE point_sources_mod