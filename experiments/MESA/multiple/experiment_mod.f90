!===============================================================================
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE solver_mod
  USE scaling_mod
  USE units_mod
  USE link_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
    procedure:: update
    procedure:: dnload
  end type
  real, save:: l_pc=4.0, n_H2=1e-20, cs_kms=0.18
CONTAINS

!===============================================================================
!> Setup an initial state.
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  real:: d, e
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /experiment_params/ l_pc, n_H2, cs_kms
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%init')
  call self%solver_t%init
  call io%lock%set
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    call io%namelist_warning ('experiment_params', iostat)
    write (io%namelist, experiment_params)
    first_time = .false.
  end if
  call io%lock%unset
  solver%csound = cs_kms*cgs%kms/scaling%u
  d = n_H2*scaling%mu*cgs%m_u/scaling%d
  if (self%gamma == 1d0) then
    e = (cs_kms*1e5/scaling%u)**2*d
  else
    e = (cs_kms*1e5/scaling%u)**2*d/(self%gamma*(self%gamma-1.0))
  end if
  self%mem(:,:,:,self%idx%d,1,1) = d
  self%mem(:,:,:,self%idx%e,1,1) = e
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  class(link_t), pointer:: link
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  call trace%begin ('experiment_t%update')
  call self%solver_t%update
  call trace%end ()
END SUBROUTINE update

!===============================================================================
!===============================================================================
SUBROUTINE dnload (self, only)
  class(experiment_t):: self
  integer, optional:: only
  class(link_t), pointer:: link
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  call trace%begin ('experiment_t%dnload')
  call self%gpatch_t%dnload (only)
  call trace%end ()
END SUBROUTINE dnload

END MODULE
