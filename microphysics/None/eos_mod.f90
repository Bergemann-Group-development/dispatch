!===============================================================================
!> empty EOS module
!===============================================================================
MODULE eos_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE mpi_mod
  USE units_mod
  USE scaling_mod
  USE kinds_mod
  implicit none
  private
  integer, parameter:: n_lambda=4
  type, public:: eos_t
    integer :: n_lambda = n_lambda
    real    :: w_lambda(n_lambda) = 1.0
    real    :: lnd_min, lnd_max, lnd_stp
    real    ::  ee_min,  ee_max,  ee_stp
    real(8) :: gamma=1.4d0
    real    :: kappa
    logical :: on = .false.
    logical :: initialized = .false.
  end type
  type (eos_t), public:: eos
END MODULE eos_mod
