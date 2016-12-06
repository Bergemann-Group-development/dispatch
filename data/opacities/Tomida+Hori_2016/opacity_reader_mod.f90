!> =============================================================================
!> Read the kappa data from Tomida & Hori 2016, which comes in two tables, with
!> Rosseland and Planck means
!> =============================================================================
MODULE opacity_reader_mod
  USE io_unit_mod
  USE trace_mod
  implicit none
  private
  type, public:: opacity_reader_t
    integer:: nT=450, nd=450, nw=2
    real(8), dimension(:), pointer:: T, d
    real(8), dimension(:,:,:), pointer:: table
  contains
    procedure:: init
  end type
  type (opacity_reader_t), public:: opacity_reader
CONTAINS

!> =============================================================================
!> Read in the data from a specific table file, which we might as well 
!> hardwire, since it is specific
!> =============================================================================
SUBROUTINE init (self)
  class (opacity_reader_t):: self
  real:: row(3)
  integer:: iT, id, iw
  !-----------------------------------------------------------------------------
  call trace%begin ('opacity_reader_t%init (Tomida+Hori_2016)')
  allocate (self%T(self%nT), self%d(self%nd), self%table(self%nT,self%nd,self%nw))
  open (io_unit%tmp+1, form='formatted', status='old', &
        file=trim(io_unit%top)//'/data/opacities/Tomida+Hori_2016/kr.dat')
  open (io_unit%tmp+2, form='formatted', status='old', &
        file=trim(io_unit%top)//'/data/opacities/Tomida+Hori_2016/kp.dat')
  do iT=1,self%nT
    do id=1,self%nd
      do iw=1,self%nw
        read (io_unit%tmp+iw,*) row
        self%table(iT,id,iw) = 10.**row(3)
      end do
      self%d(id) = 10.**row(1)
    end do
    self%T(iT) = 10.**row(2)
  end do
  close (io_unit%tmp+1)
  close (io_unit%tmp+2)
  call trace%end
END SUBROUTINE init

END MODULE opacity_reader_mod
