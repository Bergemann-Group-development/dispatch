!> =====================================================================
!> This module combines the MHD and RT solver solutions, to produce
!> time derivatives that include the net radiative cooling
!> =====================================================================
MODULE rmhd_mod
  USE io_mod
  USE trace_mod
  USE eos_mod
  USE mhd_mod
  USE rt_mod
  implicit none
  private
  ! --------------------------------------------------------------------
  type, public, extends(mhd_t):: rmhd_t
    type(rt_t)                      :: rt
    real, dimension(:,:,:), pointer :: lnd, lne
  contains
    procedure :: init
    procedure :: update
    procedure :: solve
    procedure :: interpolate
  end type
! ----------------------------------------------------------------------
CONTAINS

!> =====================================================================
!> Initialization of radiative transfer, preparing for two main tasks:
!> 1) most of the time, just inter-/extrapolate in the exising RT solutions
!> 2) on the rt%dtime intervals, solve the RT and add the new net 
!>    heating / cooling in a time slot
!> =====================================================================
SUBROUTINE init (self)
  class (rmhd_t) :: self
  integer :: nt = 5
  ! --------------------------------------------------------------------
  call trace%begin ('rmhd_t_%init')
  allocate (self%lnd (self%mesh(1)%ng, self%mesh(2)%ng, self%mesh(3)%ng), &
            self%lne (self%mesh(1)%ng, self%mesh(2)%ng, self%mesh(3)%ng))
  call self%rt%init                                                     ! setup the nbor connections
  call self%mhd_t%init
  call trace%end
END SUBROUTINE init

!> =====================================================================
!> =====================================================================
SUBROUTINE update (self)
  class (rmhd_t) :: self
  ! --------------------------------------------------------------------
  call trace%begin ('rmhd_t%update')
  call self%mhd_t%update                                                ! get df/dt from MHD
  if (.not.self%rt%is_ahead_of (self)) then                             ! if RT is not up-to-date
    call self%solve                                                     ! solve the RT for a new slot
  end if
  call self%rt%interpolate (self%time)                                  ! interpolate in time and space
  associate (dedt => self%mem(:,:,:,self%i_e,self%it,2))
  dedt = dedt + self%rt%dedt                                            ! add the RT contribution to de/dt
  end associate
  call trace%end
END SUBROUTINE update

!> =====================================================================
!> =====================================================================
SUBROUTINE solve (self)
  class (rmhd_t) :: self
  ! --------------------------------------------------------------------
  call trace%begin ('rmhd_t%solve')
  call self%interpolate (self%rt%time, self%i_d, self%lnd)              ! interpolate in time and space
  call self%interpolate (self%rt%time, self%i_e, self%lne)            ! interpolate in time and space
  self%lne = self%lne - self%lnd                                        ! internal energy per unit mass
  call eos%lookup (shape(self%lne), lne=self%lne, lnd=self%lnd, &       ! call the table lookup, ...
                   rk=self%rt%rk, src=self%rt%src)                      ! ... results go into rt_t
  call self%rt%solve                                                    ! all the RT solver to get Q
  call trace%end
END SUBROUTINE solve

!> =====================================================================
!> =====================================================================
SUBROUTINE interpolate (self, time, i, a)
  class (rmhd_t)                  :: self
  real(8)                         :: time
  integer                         :: i
  real, dimension(:,:,:), pointer :: a
  ! --------------------------------------------------------------------
  call trace%begin ('rmhd_t%interpolate')
  call trace%end
END SUBROUTINE interpolate

END MODULE rmhd_mod
