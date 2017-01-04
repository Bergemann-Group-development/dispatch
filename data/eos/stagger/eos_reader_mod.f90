!> =============================================================================
!> EOS module, reading tabular EOS compatible with STAGGER code
!> =============================================================================
MODULE eos_reader_mod
  USE io_mod
  USE io_unit_mod
  USE mpi_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Table related parameters
  !-----------------------------------------------------------------------------
  integer iupdte,nvar,mrho,mtable,njon
  real(kind=4):: eosxmin,ur,ul,ut,eps,tff,grv,abnd,dbox
  real(kind=4), allocatable, dimension(:):: tmean,tamp,rhm,xcorr,thmin,thmax, &
    dth,eemin,eemax,deetab,tab
  integer, allocatable, dimension(:):: itab,mtab
  integer, save:: nbelow=0, ncall=0, mbox
  character(len=64):: tablefile
  logical, save:: do_eos
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  ! Public part
  !-----------------------------------------------------------------------------
  type, public:: eos_reader_t
    integer:: nd, nT
    logical:: on = .false.
    real, dimension(:,:,:), pointer:: table
  contains
    procedure:: init
    procedure:: lookup
  end type
  type (eos_reader_t), public:: eos_reader
CONTAINS

!> =============================================================================
!> Initialize EOS module (only once)
!> =============================================================================
SUBROUTINE init (self)
  class (eos_reader_t):: self
  namelist /eos_params/ do_eos, tablefile
  !-----------------------------------------------------------------------------
  !$omp critical (eos_cr)
  if (first_time) then
    do_eos = .true.
    tablefile = '../../../data/eos/table.dat'
    rewind (io%input); read (io%input,eos_params)
    if (io%master) write (*,eos_params)
    !-----------------------------------------------------------------------------
    ! Read the table file
    !-----------------------------------------------------------------------------
    open (io_unit%data, file=tablefile, form='unformatted', status='old', action='read')
    read (io_unit%data) mrho,iupdte,nvar,mbox,eosxmin,dbox,ul,ut,ur,eps,tff,grv,abnd
    njon = nvar-mbox
    allocate (tmean(mrho),tamp(mrho),rhm(mrho),xcorr(mrho),thmin(mrho),thmax(mrho))
    allocate (dth(mrho),eemin(mrho),eemax(mrho),deetab(mrho),itab(mrho),mtab(mrho))
    if (io%do_trace.and.io%master) print *,mpi%rank, 'reading 12 times', mrho, ' values'
    read (io_unit%data) tmean,tamp,rhm,xcorr,thmin,thmax,dth,eemin,eemax,deetab,itab,mtab
    read (io_unit%data) mtable
    allocate (tab(mtable))
    read (io_unit%data) tab
    close (io_unit%data)
    first_time = .false.
  end if
  self%on = do_eos
  !$omp end critical (eos_cr)
  !$omp barrier
END SUBROUTINE init

END MODULE eos_reader_mod
