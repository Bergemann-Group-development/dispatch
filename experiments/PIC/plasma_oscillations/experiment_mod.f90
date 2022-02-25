!===============================================================================
!> Simple plasma oscillation experiment
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE solver_mod
  USE timer_mod
  USE link_mod
  USE math_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
    procedure:: ICs
    procedure:: update
  end type
CONTAINS

!===============================================================================
!> Setup an initial state.
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%init')
  call self%solver_t%init
  call self%ICs
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Setup an initial state, with sinusoidal Ey(x)
!===============================================================================
SUBROUTINE ICs (self)
  class(experiment_t):: self
  real, pointer, contiguous, dimension(:,:,:):: Ez, By
  integer:: ix
  real:: x
  real, save:: amplitude = 0.01
  logical, save:: first_time = .true.
  namelist /experiment_params/ amplitude
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%init')
  self%mem = 1.0
  By => self%mem(:,:,:,self%idx%by,self%it,1)
  Ez => self%mem(:,:,:,self%idx%ez,self%it,1)
  do ix=1,self%gn(1)
    x = self%position(1) + self%mesh(1)%r(ix)
    By(ix,:,:) = amplitude*sin(2.*math%pi*x/self%box(1))
    Ez(ix,:,:) = 0.0
  end do
  call trace%end()
END SUBROUTINE ICs

!===============================================================================
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  class(link_t), pointer:: link
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update')
  call self%extras_t%pre_update
  call self%extras_t%post_update
  call trace%end()
END SUBROUTINE update

END MODULE
