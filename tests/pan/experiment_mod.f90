!===============================================================================
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE gpatch_mod
  USE force_mod
  implicit none
  private
  type, public, extends(gpatch_t):: experiment_t
    type(force_t):: force
  contains
    procedure:: init
    procedure:: dnload
    procedure:: update
  end type
  integer, save:: seed=3
CONTAINS

!===============================================================================
!> Initialize patch and force
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  ! Allocate patch memory and set time step
  !----------------------------------------------------------------------------
  self%nv = 3
  self%nt = 2
  call self%patch_t%init
  self%dtime = io%out_time
  call self%force%init (self%kind, self%id, self%mesh)
END SUBROUTINE init

!===============================================================================
!> Dummy download
!===============================================================================
SUBROUTINE dnload (self, only)
  class(experiment_t):: self
  integer, optional:: only
END SUBROUTINE dnload

!===============================================================================
!> Output and update
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  real, dimension(:,:,:), pointer:: d, Ux, Uy, Uz
  real, dimension(:,:,:,:), pointer:: ff, p
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  print *,''
  !----------------------------------------------------------------------------
  ! Point the force to self%mem and compuate
  !----------------------------------------------------------------------------
  ff => self%mem(:,:,:,1:3,self%it,1)
  ff = self%force%selected (self%time, d, p, Ux, Uy, Uz, self%mesh)
  !----------------------------------------------------------------------------
  ! Output the results
  !----------------------------------------------------------------------------
  print *, 'update: time, sample =', self%time, ff(3,4,5,:)
  call self%gpatch_t%output
  call trace_end
END SUBROUTINE update

END MODULE
