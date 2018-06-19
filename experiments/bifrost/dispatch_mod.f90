!> =============================================================================
!> Module with startup and finalize procedures for DISPATCH
!> =============================================================================
MODULE dispatch_mod
  USE io_mod
  USE mpi_mod
  implicit none
  private
  type, public:: dispatch_t
  contains
    procedure:: init
    procedure:: end
  end type
  type(dispatch_t), public:: dispatch
CONTAINS

!> =============================================================================
!> Initialize DISPATCH program
!> =============================================================================
SUBROUTINE init (self)
  class(dispatch_t):: self
  !-----------------------------------------------------------------------------
  call io%init
  !call mpi%init
END SUBROUTINE init

!> =============================================================================
!> Finalize DISPATCH program
!> =============================================================================
SUBROUTINE end (self)
  class(dispatch_t):: self
  !-----------------------------------------------------------------------------
  !call mpi%end
END SUBROUTINE end

END MODULE dispatch_mod
