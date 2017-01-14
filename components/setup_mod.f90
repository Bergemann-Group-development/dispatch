!> =============================================================================
!> =============================================================================
MODULE setup_mod
  USE io_mod
  USE mpi_mod
  USE eos_mod
  USE opacity_mod
  USE scaling_mod
  implicit none
  private
  type, public:: setup_t
    contains
    procedure:: init
    procedure:: end
  end type
  type(setup_t), public:: setup
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self)
  class(setup_t):: self
  !.............................................................................
  call mpi%init                         ! MPI startup
  call io%init                          ! I/O system
  call scaling%init                     ! Scaling to code units
  call eos%init                         ! Equation-of-state
  call opacity%init                     ! Opacity
!> =============================================================================
END SUBROUTINE init
 
!> =============================================================================
!> =============================================================================
SUBROUTINE end (self)
  class(setup_t):: self
  !.............................................................................
  call mpi%end
!> =============================================================================
END SUBROUTINE end
 
END MODULE setup_mod
