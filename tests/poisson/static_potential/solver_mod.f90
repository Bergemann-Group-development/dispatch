!===============================================================================
MODULE solver_mod
  USE io_mod
  USE kinds_mod
  USE trace_mod
  USE task_mod
  USE extras_mod
  private
  type, public, extends(extras_t):: solver_t
  contains
    procedure:: init
    procedure:: update
    procedure:: gas_pressure
    procedure:: gas_velocity_scalar
    procedure:: gas_velocity_vector
    generic:: gas_velocity => gas_velocity_scalar, gas_velocity_vector
    procedure, nopass:: cast2solver
    procedure:: compression_magnitude
    procedure:: vorticity_magnitude
  end type
  type(solver_t), public:: solver

CONTAINS

!===============================================================================
!> Organize calls to the extras and the hydro solver. 
!===============================================================================
SUBROUTINE init (self)
  class(solver_t) :: self
  !.............................................................................
  call self%extras_t%pre_init
  call self%patch_t%init
  call self%gpatch_t%init
  call self%extras_t%init
END SUBROUTINE init

!===============================================================================
!> Organize calls to the extras and the hydro solver
!===============================================================================
SUBROUTINE update (self)
  class(solver_t) :: self
  !.............................................................................
  call trace%begin ('solver_t%update')
  call self%mem_lock(self%new)%set ('update')
  self%mem(:,:,:,:,self%new,1) = self%mem(:,:,:,:,self%it,1)
  call self%mem_lock(self%new)%unset ('update')
  call self%extras_t%pre_update
  call self%extras_t%post_update
  call trace%end()
END SUBROUTINE update

!===============================================================================
FUNCTION gas_pressure (self) RESULT (pg)
  class (solver_t):: self
  real(kind=KindScalarVar), dimension(self%gn(1),self%gn(2),self%gn(3)):: pg
  !-----------------------------------------------------------------------------
END FUNCTION gas_pressure

!===============================================================================
FUNCTION gas_velocity_vector (self) RESULT (v)
  class (solver_t):: self
  real(kind=KindScalarVar), dimension(self%gn(1),self%gn(2),self%gn(3),3):: v
  !-----------------------------------------------------------------------------
END FUNCTION gas_velocity_vector

!===============================================================================
FUNCTION gas_velocity_scalar (self, idir) RESULT (v)
  class (solver_t):: self
  integer:: idir
  real(kind=KindScalarVar), dimension(self%gn(1),self%gn(2),self%gn(3)):: v
  !-----------------------------------------------------------------------------
END FUNCTION gas_velocity_scalar

!===============================================================================
!> Cast a generic task_t to patch_t
!===============================================================================
FUNCTION cast2solver (task) RESULT(solver)
  class(task_t), pointer:: task
  class(solver_t), pointer:: solver
  !.............................................................................
  select type (task)
  class is (solver_t)
  solver => task
  class default
  nullify(solver)
  call io%abort ('solver_t%cast2solver: failed to cast a task to solver_t')
  end select
END FUNCTION cast2solver

!===============================================================================
SUBROUTINE compression_magnitude (self, w)
  class (solver_t):: self
  real(kind=KindScalarVar), dimension(:,:,:):: w
  !-----------------------------------------------------------------------------
END SUBROUTINE compression_magnitude

!===============================================================================
SUBROUTINE vorticity_magnitude (self, w)
  class (solver_t):: self
  real(kind=KindScalarVar), dimension(:,:,:):: w
  !-----------------------------------------------------------------------------
END SUBROUTINE vorticity_magnitude

END MODULE solver_mod
