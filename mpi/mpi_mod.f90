!> =============================================================================
!> $Id$
!> =============================================================================
MODULE mpi_mod
  USE omp_mod
  USE mpi_basic
  USE mpi_coords
  USE mpi_all
  USE mpi_reduce
  USE mpi_rma
  USE mpi_reduce
  USE mpi_send
  USE mpi_timer
  implicit none
  type box_t
    real:: size(3)
    real:: num(3)
    integer:: dim(3)
  end type
  type(box_t):: box, rank, patch
  type mpi_info_t
    integer:: size, rank, coord(3)
  contains
    procedure:: init
  end type
  type(mpi_info_t):: info
  !-----------------------------------------------------------------------------
  ! Template task, with procedure to place it it MPI-space
  !-----------------------------------------------------------------------------
  type template_t
    type(box_t):: box, mpi, patch
    real(8):: position(3)
    integer:: itask
    integer:: rank
  contains
    procedure:: init => init_template
  end type
CONTAINS

SUBROUTINE init_template (self)
  class(template_t):: self
END SUBROUTINE init_template

SUBROUTINE init (mpi)
  class(mpi_info_t):: mpi
  integer:: rank_num(3)=1, box_dim=64, patch_dim=16
  real:: box_size(3)=1.0
  namelist /mpi_params/ box_size, box_dim, rank_num, patch_dim
  !.............................................................................
  box%size   = box_size
  box%dim    = box_dim
  !
  rank%num   = rank_num
  rank%dim   = box%dim/rank%num
  rank%size  = box%size/rank%num
  !
  patch%dim  = patch_dim
  patch%num  = rank%dim/patch%dim
  patch%size = rank%size/patch%dim
  !
  if (product(rank%num) /= mpi_size) then
    print*,'rank%num, mpi_size', rank%num, mpi_size
    call abort_mpi ('wrong size')
  end if
  if (all(patch%num*patch%dim /= rank%dim)) then
    print'(3(3i8,2x))', patch%num, patch%dim, rank%dim
    call abort_mpi ('wrong size')
  end if
END SUBROUTINE init

FUNCTION get_mpi_info ()
  type(mpi_info_t):: get_mpi_info
  get_mpi_info%size  = mpi_size
  get_mpi_info%rank  = mpi_rank
  get_mpi_info%coord = mpi_coord
END FUNCTION

END MODULE mpi_mod
