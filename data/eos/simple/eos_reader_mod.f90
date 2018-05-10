!> =============================================================================
!> =============================================================================
MODULE eos_reader_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  type, public:: eos_reader_t
    integer:: nT, nd
    integer:: nv=5
    real(8), dimension(:)    , pointer:: log_T, log_d
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
  character(len=16):: file='EOS.cgs'
  real(8):: d, T
  namelist /eos_params/ file
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_reader_t%init')
  rewind (io_unit%input)
  read (io_unit%input, eos_params)
  write (io_unit%output, eos_params)
  !-----------------------------------------------------------------------------
  u = io_unit%data
  open (u, file=trim(io%top)//'/data/eos/simple/'//trim(file), &
    form='formatted', status='old')
  !-----------------------------------------------------------------------------
  read (u,'(a)')
  read (u,*) nT
  self%nT = nT
  allocate (self%log_T(nT))
  do iT=1,nT
    read (u,*) T
    self%log_T(iT) = log(T)
  end do
  print 1, 'min(T), max(T) =', exp(minval(self%log_T)), exp(maxval(self%log_T))
  1 format(1x,a,1p,2e15.3)
  !-----------------------------------------------------------------------------
  read (u,*) nd
  self%nd = nd
  allocate (self%log_d(nd))
  do id=1,nd
    read (u,*) d
    self%log_d(id) = log(d)
  end do
  print 1, 'min(d), max(d) =', exp(minval(self%log_d)), exp(maxval(self%log_d))
  !-----------------------------------------------------------------------------
  read (u,*) nv
  self%nv = nv
  allocate (self%table(nT,nd,nv))
  do id=1,nd
    do iT=1,nT
      read (u,*) self%table(iT,id,:)
      self%table(iT,id,:) = log(self%table(iT,id,:))
    end do
  end do
  close (u)
  print 1, 'min(P), max(P) =', exp(minval(self%table(:,:,1))), exp(maxval(self%table(:,:,1)))
  print 1, 'min(E), max(E) =', exp(minval(self%table(:,:,2))), exp(maxval(self%table(:,:,2)))
  call trace%end
END SUBROUTINE init

END MODULE eos_reader_mod
