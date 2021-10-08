!===============================================================================
!> Extension of short charadcteristis RT module with Multi3D specifics
!===============================================================================
MODULE m3d_rt_mod
  USE io_mod
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
    real, pointer, contiguous:: Jint(:,:,:,:) => null()
    real, pointer, contiguous:: Npop(:,:,:,:) => null()
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
  integer:: iostat
  namelist /multi3d_params/ dir
  character(128), save:: dir='input_multi3d'
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('m3d_rt_t%pre_init')
  if (first_time) then
    call io%lock%set ('multi3d_t')
    rewind (io%input)
    read (io%input, multi3d_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, multi3d_params)
    end if
    call io%namelist_warning ('multi3d_params', iostat)
    multi3d%dir = dir
    first_time = .false.
    call io%lock%unset ('multi3d_t')
  end if
  call self%rt_t%pre_init
  call multi3d%pre_init_multi3d (self%n_nu, self%n_pop)
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize: allocate arrays etc
!===============================================================================
SUBROUTINE init (self)
  class(m3d_rt_t):: self
  integer:: m(4), ng(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('m3d_rt_t%update')
  call self%rt_t%init
  if (self%on) then
    ng = (self%gn-self%n)/2
    call multi3d%init (self%n, ng, self%ds)
    m = shape(self%q)
    ! if (.not. associated(self%Jint)) then
    !   m(4) = self%n_nu
    !   allocate (self%Jint(m(1),m(2),m(3),m(4)))
    !   m(4) = self%n_pop
    !   allocate (self%Npop(m(1),m(2),m(3),m(4)))
    ! end if
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
