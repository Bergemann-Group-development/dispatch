!> A wrapper module for different mhd solvers. It seamlesly incorporates
!> radiative transfer into the solver/
MODULE mhd_rt_mod
  USE mhd_mod
  USE trace_mod
  USE microphysics_mod
  implicit none
  private
  ! --------------------------------------------------------------------  
  type, public:: rt_q_t
    real, dimension(:,:,:,:), pointer :: q   => null()
    real(kind=8), dimension(:),       pointer :: t   => null()
    integer,      dimension(:),       pointer :: iit => null()
  contains
    procedure :: rotate => rotate_q
  end type
  ! --------------------------------------------------------------------
  type, public, extends(mhd_t):: mhd_rt_t
    type(rt_q_t), pointer:: rt_q => null()
  contains
    procedure :: init_rt_heating
    procedure :: update => update_mhd_rt
    procedure :: update_q
    procedure :: apply_heating
    procedure :: add_cr_heating
  end type
! ----------------------------------------------------------------------
CONTAINS


SUBROUTINE init_rt_heating (self)
  class (mhd_rt_t) :: self
  integer :: nt = 5
  ! --------------------------------------------------------------------
  call trace_begin ('mhd_rt_t_%init_mhd_rt')
  
  ! --- initialise heating tables --- 
  allocate(self%rt_q)           
  allocate (self%rt_q%q(self%mesh(1)%gn,self%mesh(2)%gn,self%mesh(3)%gn,nt))
  allocate (self%rt_q%t(nt))
  allocate (self%rt_q%iit(nt))
  self%rt_q%q = 0.

  call trace_end
END SUBROUTINE init_rt_heating

SUBROUTINE update_mhd_rt(self)
  class (mhd_rt_t) :: self
  ! --------------------------------------------------------------------
  call trace_begin ('mhd_rt_t_%update_mhd_rt')
  ! --- Update MHD, whatever it means ---
  call self%mhd_t%update
  ! --- Update dsdt by addind Q ---
  call self%apply_heating
  
  call trace_end
END SUBROUTINE update_mhd_rt

SUBROUTINE update_q(self, tt, q, p, nwb, op, temp, dirw)
  class(mhd_rt_t) :: self
  integer :: p(3), nwb, wl
  real(kind=8) :: tt
  real, dimension(nwb) :: op, q
  real :: temp, itemp, dirw
  !---------------------------------------------------------------------
  call trace_begin ('mhd_rt_t_%update_q',4)
  if (tt > self%rt_q%t(self%it)) call self%rt_q%rotate(tt)
  
  itemp = 1./temp
  do wl =1, nwb
    self%rt_q%q(p(1),p(2),p(3),self%it) = &
    self%rt_q%q(p(1),p(2),p(3),self%it) + dirw*op(wl)*q(wl)*itemp
  end do
  call trace_end
END SUBROUTINE update_q

SUBROUTINE add_cr_heating (self,tt, temp)
  class(mhd_rt_t) :: self
  real, dimension(:,:,:) :: temp
  real(kind=8) :: tt
  integer :: l(3), u(3)
  !---------------------------------------------------------------------
  call trace_begin ('mhd_rt_t_%add_cr_heating')
  l = self%mesh%li
  u = self%mesh%ui
  associate (rho => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),1,self%it,1), &
               q => self%rt_q%q(l(1):u(1),l(2):u(2),l(3):u(3),self%it), &
               t => temp)
           
  q = q + rho * microphysics%cr_heating / t
  
  end associate
  call trace_end

END SUBROUTINE add_cr_heating

SUBROUTINE apply_heating(self)
  class(mhd_rt_t) :: self
  integer:: l(3), u(3)
  !---------------------------------------------------------------------
  call trace_begin ('mhd_rt_t_%apply_heating')
  l = self%mesh%li
  u = self%mesh%ui
  !print *, 'DSDT now is... '
  !print *,self%mem(l(1):u(1),l(2):u(2),l(3):u(3),1,self%it,2)
  self%mem(l(1):u(1),l(2):u(2),l(3):u(3),1,self%it,2) = &
  self%mem(l(1):u(1),l(2):u(2),l(3):u(3),1,self%it,2) + &
  self%rt_q%q(l(1):u(1),l(2):u(2),l(3):u(3),self%it)
  !print *, 'heating.. '
  !print *,self%rt_q%q(l(1):u(1),l(2):u(2),l(3):u(3),self%it)
  !print *, 'after heating was applied:'
  !print *,self%mem(l(1):u(1),l(2):u(2),l(3):u(3),1,self%it,2)
  !print *, 'AND S IS ORIGINALY:'
  !print *,self%mem(l(1):u(1),l(2):u(2),l(3):u(3),1,self%it,1)
  !error stop "stopping"
  call trace_end
END SUBROUTINE apply_heating

SUBROUTINE rotate_q(self, tt)
  class(rt_q_t) :: self
  real (kind=8) :: tt
! ----------------------------------------------------------------------

END SUBROUTINE rotate_q

END MODULE mhd_rt_mod
