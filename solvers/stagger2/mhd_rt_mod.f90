!> =============================================================================
!> Wrapper module, which seamlessly incorporates RT into the stagger2/ solver.
!> =============================================================================
MODULE mhd_rt_mod
  USE mhd_mod
  USE trace_mod
  USE microphysics_mod
  implicit none
  private
  ! --------------------------------------------------------------------
  type, public:: rt_q_t
    real, dimension(:,:,:,:),   pointer :: q   => null()
    real(kind=8), dimension(:), pointer :: t   => null()
    integer,      dimension(:), pointer :: iit => null()
  end type
  ! --------------------------------------------------------------------
  type, public, extends(mhd_t):: mhd_rt_t
    type(rt_q_t), pointer:: rt_q => null()
  contains
    procedure :: init_rt
    procedure :: update
    procedure :: apply_heating
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
  allocate (self%rt_q)
  allocate (self%rt_q%q(self%mesh(1)%gn,self%mesh(2)%gn,self%mesh(3)%gn,nt))
  allocate (self%rt_q%t(nt))
  allocate (self%rt_q%iit(nt))
  self%rt_q%q = 0.
  call trace%end
END SUBROUTINE init_rt

!> =============================================================================
!> Update the MHD state and add the RT head exchange
!> =============================================================================
SUBROUTINE update (self)
  class (mhd_rt_t) :: self
  ! ----------------------------------------------------------------------------
  call trace%begin ('mhd_rt_t_%update')
  ! --- Update MHD, whatever it means ---
  call self%mhd_t%update
  ! --- Update dsdt by adding Q ---
  call self%apply_heating
  call trace%end
END SUBROUTINE update

!> =============================================================================
!> Add the heating to the MHD equations; this is the only part that is expected
!> to differ from solver to solver
!> =============================================================================
SUBROUTINE apply_heating (self)
  class(mhd_rt_t) :: self
  integer:: l(3), u(3)
  ! ----------------------------------------------------------------------------
  call trace%begin ('mhd_rt_t_%apply_heating')
  l = self%mesh%li
  u = self%mesh%ui

  self%mem(l(1):u(1),l(2):u(2),l(3):u(3),2,self%it,2) = &
  self%mem(l(1):u(1),l(2):u(2),l(3):u(3),2,self%it,2) + &
  self%rt_q%q(l(1):u(1),l(2):u(2),l(3):u(3),self%it)

  call trace%end
END SUBROUTINE apply_heating

END MODULE mhd_rt_mod
