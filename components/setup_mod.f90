!===============================================================================
!> Dependence hierarchy of main service modules, initialize here:
!>
!> io
!>   mpi
!>     trace
!>       timer
!>         omp_lock
!>           omp_timer
!>           omp
!===============================================================================
MODULE setup_mod
  USE io_mod
  USE mpi_mod
  USE omp_lock_mod
  USE timer_mod
  USE trace_mod
  USE mpi_mesg_mod
  USE scaling_mod
  USE eos_mod
  USE opacity_mod
  USE refine_mod
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
!> Initialise basic components.  Note that since omp_lock_mod uses unit numbers
!> and directory names defined by io_mod and stored in io_unit_mod, its init
!> must called after io%init.
!===============================================================================
SUBROUTINE init (self)
  class(setup_t):: self
  call io%init                          ! I/O system
  call omp_lock%init ('setup')          ! needs omp, io%init
  call timer%init                       ! needs locks
  call mpi_mesg%init                    ! needs locks
  call trace%begin ('setup_t%init')     ! needs timer
  call io%header ('scaling_t%init:')
  call scaling%init                     ! Scaling to code units
  call eos%init                         ! Equation-of-state
  call opacity%init                     ! Opacity
  call refine%init                      ! adaptive mesh refinement
  call trace%end ()
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
