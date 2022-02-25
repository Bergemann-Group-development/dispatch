!===============================================================================
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE solver_mod
  USE math_mod
  USE mesh_mod
  USE kinds_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: update
    procedure:: bndries
  end type
  logical, save:: first_time = .true.
CONTAINS

!===============================================================================
!> pre-initialize the hydro solver. 
!===============================================================================
SUBROUTINE pre_init (self)
  class(experiment_t), target :: self
  !.............................................................................
  call trace%begin ('solver_t%pre_init')
  call self%solver_t%pre_init              ! calls self%mhd_t%pre_init
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Setup initial states in the patches
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t)                 :: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init')
  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  self%mhd = .true.
  call self%solver_t%init
  call trace%end
END SUBROUTINE init

!===============================================================================
!> Boundary conditions
!===============================================================================
SUBROUTINE bndries (self)
  class(experiment_t):: self
  real(var_k), dimension(:,:,:,:), pointer, contiguous:: mem
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%bndries')
  mem => self%mem(:,:,:,:,self%it,1)
  call self%boundaries%primary (mem)
  call trace%end()
END SUBROUTINE bndries

!===============================================================================
!> Solver update; the default should be built into the solver
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call self%solver_t%update
  call trace_end
END SUBROUTINE update

END MODULE
