!===============================================================================
!> Start-up and shutdown of the run's service modules.
!>
!> The service modules form a chain of USE statements: io_mod uses mpi_mod, which
!> uses trace_mod, then timer_mod, omp_lock_mod, and omp_timer_mod and omp_mod.
!===============================================================================
MODULE setup_mod
  USE io_mod
  USE mpi_mod
  USE mpi_mesg_mod
  USE scaling_mod
  USE trace_mod
  USE job_record_mod
  implicit none
  private
  !===============================================================================
  !> Start-up and shutdown sequence of a run.
  !===============================================================================
  type, public:: setup_t
    contains
    procedure:: init
    procedure:: end
  end type
  type(setup_t), public:: setup
CONTAINS

!===============================================================================
!> Initialises the I/O system, the job record, MPI messaging, tracing and the code
!> units, in that order. omp_lock_mod uses unit numbers and directory names that
!> io_mod defines and stores in io_unit_mod, so its init must be called after
!> io%init.
!===============================================================================
SUBROUTINE init (self)
  class(setup_t):: self
  call io%init                          ! I/O system
  call job_record%start                 ! needs io%run_dir
  call mpi_mesg%init                    ! needs locks
  call trace%begin ('setup_t%init')     ! needs timer
  call scaling%init                     ! Scaling to code units
  call trace%end ()
END SUBROUTINE init
 
!===============================================================================
!> Marks the job record finished and ends MPI. A crash never gets here, so its
!> job.nml keeps state='running'.
!===============================================================================
SUBROUTINE end (self)
  class(setup_t):: self
  !.............................................................................
  call trace%begin ('setup_t%end')
  call job_record%finish                ! a normal exit; a crash never gets here
  call mpi%end
  call trace%end()
END SUBROUTINE end
 
END MODULE setup_mod
