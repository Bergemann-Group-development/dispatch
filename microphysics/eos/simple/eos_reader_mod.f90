!> =============================================================================
!> =============================================================================
MODULE eos_reader_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  type, public:: eos_reader_t
    integer:: skip=1
    integer:: nT, nd
    real(8), dimension(:), pointer:: T, d
    real(8), dimension(:,:,:), pointer:: table
    character(len=64):: dir='data/eos'
  contains
    procedure:: init
  end type
  type (eos_reader_t), public:: eos_reader
CONTAINS

!> =============================================================================
!> Read in the data from a specific table file, which we might as well 
!> hardwire, since it is specific
!> =============================================================================
SUBROUTINE init (self)
  class (eos_reader_t):: self
  integer:: u, nv, nT, nd, iT, id, i
  integer:: skip=1
  character(len=16):: file='EOS.cgs'
  namelist /eos_params/ file, skip
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_reader_t%init')
  open (io_unit%data, form='formatted', status='old', &
    file=trim(io%top)//'/data/eos/simple/'//trim(file))
  rewind (io_unit%input)
  read (io_unit%input, eos_params)
  write (io_unit%output, eos_params)
  !-----------------------------------------------------------------------------
  u = io_unit%data
  do i=1,self%skip
    read (u,'(a)')
  end do
  !-----------------------------------------------------------------------------
  read (u,*) nT
  self%nT = nT
  allocate (self%T(nT))
  do iT=1,nT
    read (u,*) self%T(iT)
  end do
  !-----------------------------------------------------------------------------
  read (u,*) nd
  self%nd = nd
  allocate (self%d(nd))
  do id=1,nd
    read (u,*) self%d(id)
  end do
  !-----------------------------------------------------------------------------
  read (u,*) nv
  print*,'MK', self%d(nd), self%T(nt), nv
  allocate (self%table(nT,nd,nv))
  do id=1,nd
    do iT=1,nT
      read (u,*) self%table(iT,id,:)
    end do
  end do
  close (u)
  call trace%end
END SUBROUTINE init

END MODULE eos_reader_mod
