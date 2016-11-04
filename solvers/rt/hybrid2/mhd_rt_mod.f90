!> =============================================================================
!> Wrapper module, which seamlessly incorporates RT into the stagger2/ solver.
!> =============================================================================
MODULE mhd_rt_mod
  USE mhd_mod
  USE trace_mod
  USE microphysics_mod
  USE rt_data_mod
  implicit none
  private
  ! --------------------------------------------------------------------
  type, public, extends(mhd_t):: mhd_rt_t
    type(rt_data_t), pointer:: rt_data => null()
  contains
    procedure :: init => init_rt
    procedure :: update
  end type
  ! ----------------------------------------------------------------------------
CONTAINS

!> =============================================================================
!> Initialize the underlying MHD patch & allocate the extra memory needed for RT
!> =============================================================================
SUBROUTINE init_rt (self)
  class (mhd_rt_t) :: self
  integer :: nt = 5
  ! ----------------------------------------------------------------------------
  call trace%begin ('mhd_rt_t_%init')
  self%type = 'mhd_rt_t'
  ! --- initialise heating tables ---
  call self%mhd_t%init
  print *,"initializing mhd_rt", self%time
  call trace%end
END SUBROUTINE init_rt

!> =============================================================================
!> Update the MHD state and add the RT head exchange
!> =============================================================================
SUBROUTINE update (self)
  class (mhd_rt_t) :: self
  real, dimension(:,:,:), pointer:: q
  ! ----------------------------------------------------------------------------
  call trace%begin ('mhd_rt_t_%update')
  ! --- Update MHD, whatever it means ---
  call self%mhd_t%update
  !print *,'updating mhd',self%time
  ! --- Update dsdt by adding Q ---
  q => self%rt_data%q(:,:,:,self%it)  ! FIXME: should be interpolation in time
  !call self%mhd_t%apply_heating (q)
  call trace%end
END SUBROUTINE update

END MODULE mhd_rt_mod
