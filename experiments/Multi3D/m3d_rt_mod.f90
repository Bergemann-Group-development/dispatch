!===============================================================================
!> Extension of short charadcteristis RT module with Multi3D specifics
!===============================================================================
MODULE m3d_rt_mod
  USE trace_mod
  USE rt_mod
  USE multi3d_mod
  implicit none
  private
  type, public, extends(rt_t):: m3d_rt_t
    !---------------------------------------------------------------------------
    ! Before adding things here, check of the variable already exists in rt_t,
    ! and consider if it should be shared, or Multi3D specific (new name)
    !---------------------------------------------------------------------------
    integer:: n_nu, n_pop  ! ...
    real:: task_variable
    !---------------------------------------------------------------------------
    ! Add arrays here to store population numbers are whatever else
    !---------------------------------------------------------------------------
    real, pointer, contiguous:: array(:,:,:,:) => null()
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: update
    procedure:: output
  end type
CONTAINS

!===============================================================================
!> Pre-initialize: read parameters, if any.  Probably we can rely on the
!> existing parameters inside Multi3D and DISPATCH.
!===============================================================================
SUBROUTINE pre_init (self)
  class(m3d_rt_t), target:: self
  !-----------------------------------------------------------------------------
  call trace%begin ('m3d_rt_t%pre_init')
  call self%rt_t%pre_init
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize: allocate arrays etc
!===============================================================================
SUBROUTINE init (self)
  class(m3d_rt_t):: self
  integer:: m(4)
  !-----------------------------------------------------------------------------
  call trace%begin ('m3d_rt_t%update')
  call self%rt_t%init
  m = shape(self%q)
  m(4) = self%n_pop
  if (.not. associated(self%array)) then
    allocate (self%array(m(1),m(2),m(3),m(4)))
    !allocate (...)
  end if
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Update the population numbers and other dynamic task data.  Consider spliting
!> this into three steps; put data to Multi3D in pre_update, update, and get
!> data back from Mulit3D in post_update.
!===============================================================================
SUBROUTINE update (self)
  class(m3d_rt_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('m3d_rt_t%update')
  call self%rt_t%update
  call multi3d%pre_update ! (... argument list ...)
  call multi3d%update
  call multi3d%post_update ! (... argument list ...)
  call trace%end ()
END SUBROUTINE update

!===============================================================================
!> Do Multi3D specific output
!===============================================================================
SUBROUTINE output (self)
  class(m3d_rt_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('m3d_rt_t%output')
  !call multi3d%output
  call trace%end ()
END SUBROUTINE output

END MODULE m3d_rt_mod
