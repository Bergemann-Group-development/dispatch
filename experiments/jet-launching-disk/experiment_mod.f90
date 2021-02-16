!===============================================================================
!> Initializiation and update of each ZEUS patch
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE solver_mod
  USE pboundaries_mod
  USE pebble_patch_mod
  USE potential_mod
  USE kinds_mod
  USE eos_mod
  implicit none
  private
  type, public, extends(pebble_patch_t):: experiment_t
  contains
    procedure:: init => init_experiment
    procedure:: update => update_experiment
    procedure:: density_floor => enforce_d_floor
  end type
CONTAINS

!===============================================================================
!> Initialise this experiment.
!> Things that need to be initialised per-patch happens here.
!===============================================================================
SUBROUTINE init_experiment (self)
  class(experiment_t) :: self
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%init')

  self%mhd = .true.
  self%initial%mhd = .true.
  call self%initial%init(self%kind, real(eos%gamma))
  call self%solver_t%init()
  call self%pebble_patch_t%init()
  ! intialize time-independent gravitational potential of the central star
  call potential%init(self%mesh, self%idx, self%initial%star_mass, self%mem(:,:,:,:,1,1))
  call pboundaries%init_bdries (self)
  call self%density_floor(self%it)
  call self%check_density()

  call trace%end()
END SUBROUTINE init_experiment

!===============================================================================
!> Update procedure for this experiment.
!===============================================================================
SUBROUTINE update_experiment (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update')

  call pboundaries%apply(self, self%mhd, self%it)
  call self%pebble_patch_t%output()
  if (self%kind == 'zeus_mhd_patch' .or. self%kind == 'zeus_tw_mhd_patch') then
    if (self%mesh(1)%lower_boundary) self%btypes%iib = 9
    if (self%mesh(1)%upper_boundary) self%btypes%oib = 9
    if (self%mesh(2)%lower_boundary) self%btypes%ijb = 1
    if (self%mesh(2)%upper_boundary) self%btypes%ojb = 1
  end if
  call self%solver_t%update()
  call pboundaries%apply(self, self%mhd, self%new, previt=self%it)
  call self%density_floor(self%new)
  if (self%do_pebbles) call self%pebble_patch_t%update_pebbles()

  call trace%end()
END SUBROUTINE update_experiment

!===============================================================================
!> Apply a density floor to the current patch.
!> The density floor value is stored in `initial_t`.
!===============================================================================
SUBROUTINE enforce_d_floor (self, cit)
  class(experiment_t), target :: self
  integer, intent(in):: cit

  integer:: i1, i2, i3, l(3), u(3)
  real(8):: dfloor
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_mod::enforce_d_floor')

  l = self%mesh%lb
  u = self%mesh%ub

  dfloor = self%initial%density_floor
  d => self%mem(:,:,:,self%idx%d,cit,1)
  do i3=l(3),u(3)
    do i2=l(2),u(2)
      do i1=l(1),u(1)
        d(i1,i2,i3) = max(d(i1,i2,i3), dfloor)
      end do
    end do
  end do

  nullify (d)
  call trace%end()
END SUBROUTINE enforce_d_floor

END MODULE experiment_mod
