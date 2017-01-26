!> =============================================================================
!> A generic interface to an external code, with the data to communicate between
!> DISPATCH and the code.
!>
!> With this interface, the experiments_mod does not need to access any of the
!> codes parameters or data, and the code/dispatch_mod does not need to access
!> any DISPATCH data directly.  This makes for cleaner code, and avoids potential
!> collisions between global names.  Module name collions, !> if any, can then
!> be avoided by modifying the make setup, so external files are compiled "in
!> place", rather than in the DISPATCH experiment directory.
!> =============================================================================
MODULE interface_mod
  implicit none
  public
  type:: interface_t
    real(4), pointer:: mem(:,:,:,:,:)
    real(8), pointer:: position(:), size(:)
    real(8), pointer:: time, dtime
    integer, pointer:: istep
    integer, pointer:: id
    integer, pointer:: mpi_rank
    integer, pointer:: log_unit
    integer, pointer:: input_unit
    real,    pointer:: u_max
  end type
  type(interface_t), save:: interface
  !$omp threadprivate (interface)
END MODULE interface_mod
