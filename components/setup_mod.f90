!===============================================================================
!===============================================================================
MODULE setup_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE eos_mod
  USE opacity_mod
  USE scaling_mod
  USE timer_mod
  USE download_mod
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
!> Initialise basic components.
!===============================================================================
SUBROUTINE init (self)
  class(setup_t):: self
  !.............................................................................
  call mpi%init                         ! MPI startup
  call io%init                          ! I/O system
  call trace%begin ('setup_t%init')
  call scaling%init                     ! Scaling to code units
  call eos%init                         ! Equation-of-state
  call opacity%init                     ! Opacity
  call download%init                    ! Downloading and interpolation parameters
  call trace%end()
END SUBROUTINE init
 
!===============================================================================
!> The run is finished; tidy up.
!> For the moment, only MPI needs to ended; this would be a good place to put
!> hooks for other finalisation tasks.
!===============================================================================
SUBROUTINE end (self)
  class(setup_t):: self
  !.............................................................................
  call trace%begin ('setup_t%end')
  call mpi%end
  call trace%end()
END SUBROUTINE end
 
END MODULE setup_mod
