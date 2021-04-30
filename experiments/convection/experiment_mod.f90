!===============================================================================
!> Interface on top of solvers and boundary conditions.  The init and update
!> entries are redundant in this case, but the bndries() entry is not.
!>
!> The task%bndries() procedure is called at the same time as guard zones are
!> loaded, and should do whatever is required to set the boundaries -- it is
!> not reqiured that it calls a framework procedure as part of this.  However,
!> a number of standard boundary conditions are available in the subsystem
!> that has standard_boundaries available for inheritence.
!>
!> A standard way to implement BCs is thus to create a boundaries_mod.f90, and
!> let it inherit standard_boundaries_t% and its methods.  Here, the entry
!> point %primary() is called from task%bndries(), and the actual BCs are
!> implented there.
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE solver_mod
  USE kinds_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
    procedure:: bndries
    procedure:: update
  end type
  logical, save:: first_time=.true.
CONTAINS

!===============================================================================
!> Setup initial states in the patches
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init')
  call self%solver_t%init
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Setup initial states in the patches
!===============================================================================
SUBROUTINE bndries (self)
  class(experiment_t):: self
  real(var_k), dimension(:,:,:,:), pointer, contiguous:: mem
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%bndries')
  mem => self%mem(:,:,:,:,self%new,1)
  call self%boundaries%primary (mem)
  call trace%end()
END SUBROUTINE bndries

!===============================================================================
!> Setup initial states in the patches
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call self%solver_t%update
  call trace%end ()
END SUBROUTINE update

END MODULE
