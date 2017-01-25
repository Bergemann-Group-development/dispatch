!> =============================================================================
!> A generic interface to an external code, which the data to communicate
!> =============================================================================
MODULE interface_mod
  implicit none
  type, public:: interface_t
    real(4), pointer:: mem(:,:,:,:)
    real(8), pointer:: position(:), size(:)
    real(8), pointer:: time, dtime
    integer, pointer:: istep
    integer, pointer:: id
    integer, pointer:: mpi_rank
    integer, pointer:: log_unit
    real,    pointer:: u_max
  end type
  type(interface_t):: interface
  !$omp threadprivate (interface)
END MODULE interface_mod
