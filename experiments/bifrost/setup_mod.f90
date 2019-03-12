!===============================================================================
!> Module with startup and finalize procedures for DISPATCH
!===============================================================================
MODULE setup_mod
  USE io_mod
  USE mpi_mod
  USE extern_mod
  implicit none
  private
  type, public:: setup_t
  contains
    procedure:: init
    procedure:: end
  end type
  type(setup_t), public:: setup
CONTAINS

!===============================================================================
!> Initialize DISPATCH program
!===============================================================================
SUBROUTINE init (self)
  class(setup_t):: self
  !-----------------------------------------------------------------------------
  call io%init
  call extern%init
  !call mpi%init

END SUBROUTINE init

!===============================================================================
!> Finalize DISPATCH program
!===============================================================================
SUBROUTINE end (self)
  class(setup_t):: self
  !-----------------------------------------------------------------------------
  !call mpi%end
  call extern%end

END SUBROUTINE end

END MODULE setup_mod
