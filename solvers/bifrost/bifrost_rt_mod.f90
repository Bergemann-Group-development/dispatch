!> =============================================================================
!> Interface module, with a data type that can store all required BIFROST RT
!> details as part of the solver_t data type, for get/put operations.
!> =============================================================================
MODULE bifrost_rt_mod
  ! Dispatch modules
  USE trace_mod
  ! Bifrost modules
  USE params_mod
  USE radiation_mod
  implicit none
  private
  type, public:: bifrost_rt_t
    real, dimension(:,:,:,:,:), allocatable:: SIboundxupin, SIboundxdnin
    real, dimension(:,:,:,:,:), allocatable:: SIboundyupin, SIboundydnin
    real, dimension(:,:,:,:,:), allocatable:: SIboundzupin, SIboundzdnin
    real, dimension(:,:,:,:  ), allocatable:: S
    real, dimension(:,:,:    ), allocatable:: qtot
    real(kind=8)                           :: trad
    integer                                :: nrefine
    logical :: do_rt = .false.
  contains
    procedure:: init
    procedure:: put
    procedure:: get
    procedure:: refine
  end type
  character(len=80), save:: id = &
  '$Id$'
CONTAINS

!> =============================================================================
!> Allocate patch space for Bifrost RT data
!> =============================================================================
SUBROUTINE init (self, mb)
  class(bifrost_rt_t):: self
  integer, optional:: mb(3)
  integer:: m(5), n(4), o(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('bifrost_rt_t%init')
  if (allocated(SIboundxupin)) then
self%do_rt = .true.
    m=shape(SIboundxupin); allocate (self%SIboundxupin(m(1),m(2),m(3),m(4),m(5)))
  end if
  if (allocated(SIboundyupin)) then
    m=shape(SIboundyupin); allocate (self%SIboundyupin(m(1),m(2),m(3),m(4),m(5)))
  end if
  if (allocated(SIboundzupin)) then
    m=shape(SIboundzupin); allocate (self%SIboundzupin(m(1),m(2),m(3),m(4),m(5)))
  end if
  if (allocated(SIboundxdnin)) then
    m=shape(SIboundxdnin); allocate (self%SIboundxdnin(m(1),m(2),m(3),m(4),m(5)))
  end if
  if (allocated(SIboundydnin)) then
    m=shape(SIboundydnin); allocate (self%SIboundydnin(m(1),m(2),m(3),m(4),m(5)))
  end if
  if (allocated(SIboundzdnin)) then
    m=shape(SIboundzdnin); allocate (self%SIboundzdnin(m(1),m(2),m(3),m(4),m(5)))
  end if
  if (allocated(qtot)) then
    o=shape(qtot)        ; allocate (self%qtot(o(1),o(2),o(3)))
  end if
  if (allocated(S)) then
    n=shape(S)           ; allocate (self%S(n(1),n(2),n(3),n(4)))
  end if
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Put the patch values back into the Bifrost variables
!> =============================================================================
SUBROUTINE put (self)
  class(bifrost_rt_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('bifrost_rt_t%put')
  if (self%do_rt) then
    SIboundxupin = self%SIboundxupin
    SIboundyupin = self%SIboundyupin
    SIboundzupin = self%SIboundzupin
    SIboundxdnin = self%SIboundxdnin
    SIboundydnin = self%SIboundydnin
    SIboundzdnin = self%SIboundzdnin
    S            = self%S
    qtot         = self%qtot
    trad         = self%trad
    nrefine      = self%nrefine
  end if
  call trace%end
END SUBROUTINE put

!> =============================================================================
!> Copy the Bifrost RT data to patch memory
!> =============================================================================
SUBROUTINE get (self)
  class(bifrost_rt_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('bifrost_rt_t%get')
  if (self%do_rt) then
    self%SIboundxupin = SIboundxupin
    self%SIboundyupin = SIboundyupin
    self%SIboundzupin = SIboundzupin
    self%SIboundxdnin = SIboundxdnin
    self%SIboundydnin = SIboundydnin
    self%SIboundzdnin = SIboundzdnin
    self%S            = S
    self%qtot         = qtot
    self%trad         = trad
    self%nrefine      = nrefine
  end if
  call trace%end
END SUBROUTINE get

SUBROUTINE refine (self)
  class(bifrost_rt_t):: self
  !-----------------------------------------------------------------------------
END SUBROUTINE refine

END MODULE bifrost_rt_mod
