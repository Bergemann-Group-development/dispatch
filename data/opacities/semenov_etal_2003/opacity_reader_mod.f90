!> =============================================================================
!> =============================================================================
MODULE opacity_reader_mod
  USE io_mod
  USE trace_mod
  USE opacity_mod
  USE scaling_mod
  implicit none
  private
  type, public:: opac_reader_t
    integer:: skip=1
    integer:: nT=0, nd=0, nw=1
    real, dimension(:), pointer:: T, d
    real, dimension(:,:,:), pointer:: table
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
  integer, save:: nT=100, nd=100, nw=1
  real, save:: Tmin=5., Tmax=3e4, dmin=1e-16, dmax=1e-9
  real:: tmp
  integer:: i, j
  namelist /semenov_params/ nT, nd, Tmin, Tmax, dmin, dmax
  allocate (self%T(nT), self%d(nd), self%table(nT,nd,nw))
  do i=1,nT
    self%T(i) = exp(alog(Tmin)+alog(Tmax/Tmin)*(i-1.0))/(nT-1.0)
  end do
  do j=1,nd
    self%d    = exp(alog(dmin)+alog(dmax/dmin)*(j-1.0))/(nd-1.0)
    self%table(:,j,1) = opacity%mean(self%d,self%T)
  end do
  do j=1,nd
    self%d(j) = exp(alog(dmin)+alog(dmax/dmin)*(j-1.0))/(nd-1.0)
  end do
END SUBROUTINE init

END MODULE opacity_reader_mod
