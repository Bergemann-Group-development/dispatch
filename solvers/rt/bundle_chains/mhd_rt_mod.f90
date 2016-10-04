!> A wrapper module for different mhd solvers. It seamlesly incorporates
!> radiative transfer into the solver/
MODULE mhd_rt_mod
  USE mhd_mod
  USE trace_mod
  USE microphysics_mod
  USE rt_data_mod
  USE rt_internal_mod
  implicit none
  private
  ! --------------------------------------------------------------------
  type, public, extends(mhd_t):: mhd_rt_t
    type(rt_data_t), pointer:: rt_data => null()
    type(rt_internal_t), pointer:: first_internal
    integer:: BC_direction=0
  contains
    procedure :: init_rt
    procedure :: update
    procedure :: bundle_chain
  end type
! ----------------------------------------------------------------------
CONTAINS

!> =====================================================================
!> This is currently a place holder for setting up the RT for one patch
!> =====================================================================
SUBROUTINE init_rt (self)
  class (mhd_rt_t) :: self
  integer :: nt = 5
  ! --------------------------------------------------------------------
  !call trace_begin ('mhd_rt_t_%init_mhd_rt')
  !call trace_end
END SUBROUTINE init_rt

!> =====================================================================
!> This updates an RMHD patch, by first computing RT data, then starting
!> the bundle chains that do the actual solutions, then storing the
!> heating inside the MHD task, which is assumed to adjuste the MHD
!> update in a method specific way.
!> =====================================================================
SUBROUTINE update(self)
  class (mhd_rt_t) :: self
  real, dimension(:,:,:), pointer:: q
  ! --------------------------------------------------------------------
  call trace_begin ('mhd_rt_t_%update')
  call self%rt_data%update
  call self%bundle_chain
  q => self%rt_data%q(:,:,:,self%it)  ! FIXME: should be interpolation in time
  call self%mhd_t%apply_heating (q)
  call self%mhd_t%update
  call trace_end
END SUBROUTINE update

!> =====================================================================
!> This gets called for all patches, but for most of them, which are
!> not boundary patches, the first_internal pointer is null, so the call
!> returns immediately.  The remaining cases have a linked list of tasks
!> that need BCs, and then execute the first task in a chain, which then
!> executes the next, recursively.
!> =====================================================================
SUBROUTINE bundle_chain (self)
  class (mhd_rt_t) :: self
  class (rt_internal_t), pointer :: task
  ! --------------------------------------------------------------------
  call trace_begin ('mhd_rt_t_%bundle_chain')
  task => self%first_internal
  do while (associated(task))
    !$omp task firstprivate (task)
    !call task%get_BC()
    call task%recurse
    !$omp end task
    task => task%next
  end do
  !$omp taskwait
  call trace_end
END SUBROUTINE bundle_chain

END MODULE mhd_rt_mod
