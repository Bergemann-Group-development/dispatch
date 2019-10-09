!=======================================================================
!> This module combines the MHD and RT solver solutions, to produce
!> time derivatives that include the net radiative cooling
!>
!> Radiative transfer module structure:
!>
!> experiment_mod                       ! generic time advance modul
!>   rmhd_mod                           ! wrapper for MHD and RT
!>     mhd_mod                          ! MHD solver
!>       eos_mod                        ! computing gas pressure
!>     rt_mod                           ! RT organizer
!>       eos_mod                        ! computing opacity and source fn
!>       rt_bundle_mod                  ! handles bundles of rays
!>         rt_solver_mod                ! solver the RT equation
!=======================================================================
MODULE rmhd_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE eos_mod
  USE solver_mod
  USE rt_mod
  USE timestep_mod
  implicit none
  private
  ! --------------------------------------------------------------------
  type, public, extends(solver_t):: rmhd_t
    class(rt_t), pointer      :: rt
  contains
    procedure :: init
    procedure :: update
  end type
! ----------------------------------------------------------------------
CONTAINS

!=======================================================================
!> Initialization of radiative transfer, preparing for two main tasks:
!> 1) most of the time, just inter-/extrapolate in the exising RT solutions
!> 2) on the rt%dtime intervals, solve the RT and add the new net
!>    heating / cooling in a time slot
!=======================================================================
SUBROUTINE init (self)
  class (rmhd_t) :: self
  !---------------------------------------------------------------------
  call trace%begin ('rmhd_t_%init')
  !---------------------------------------------------------------------
  ! Initialize the MHD data type, and allocate arrays for RT
  !---------------------------------------------------------------------
  call self%mhd_t%init
  !---------------------------------------------------------------------
  ! Initialize the rt_mod object (one could create more than one).
  ! Copy parameters that may be different for various reasons
  !---------------------------------------------------------------------
  self%mesh%p      = self%position
  call trace%end()
END SUBROUTINE init

!=======================================================================
!> If RT is on, add the radiative cooling / heating to the entropy time
!> derivative, before taking the timestep.
!=======================================================================
SUBROUTINE update (self)
  class (rmhd_t) :: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: dedt, e, qrad
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: dddt, d
  real:: u_max
  integer:: iy
  integer, save :: itimer=0
  !---------------------------------------------------------------------
  call trace%begin ('rmhd_t%update', itimer=itimer)
  call self%output
  call self%pde
  associate(l => self%li, u => self%ui)
  if (eos%do_rt.and.associated(self%rt)) then
    if (.not.self%rt%do_chains) call self%rt%update
    allocate(qrad(self%n(1), self%n(2), self%n(3)))
    qrad = 0.0
    call self%rt%get_q (self%time, qrad)                                  ! interpolate in time and space
    !-----------------------------------------------------------------------------
    ! Time step limitation for rate of change dS/dt = (ds/dt-(s/d)*dddt)/d
    !-----------------------------------------------------------------------------
    d    => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%d,self%it,1)
    e    => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e,self%it,1)
    dddt => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%d,self%it,2)
    dedt => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e,self%it,2)
    dedt = dedt + qrad
    if (io%verbose>2) then
      print *,'mk:    d', minval(d), maxval(d)
      print *,'mk: dedt', minval(dedt), maxval(dedt)
      print *,'mk: dddt', minval(dddt), maxval(dddt)
    else if (io%verbose>1) then
      print *,'rmhd_t%update: minval(qrad), maxval(qrad) =', minval(qrad), maxval(qrad)
      do iy=1,self%n(2)
        print *,'rmhd_t%update: min,aver,max(qrad) =', iy, &
          minval(qrad(:,iy,:)/e(:,iy,:)), sum(qrad(:,iy,:))/sum(e(:,iy,:)), &
          maxval(qrad(:,iy,:)/e(:,iy,:))
      end do
    end if
    u_max = maxval(self%ds)*maxval(abs((dedt-e*dddt/d)/d))*self%courant/self%cdtd
    self%u_max = max(self%u_max,u_max)
    if (io%verbose>1) print 1, 'u_max(RT)  ', u_max, self%u_max
    1 format(1x,a,1p,2g12.3)
    !---------------------------------------------------------------------------
    ! Store qrad in temperatures slot for now
    !---------------------------------------------------------------------------
    self%mem(:,:,:,self%idx%tt,self%it,1) = 0.0
    self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%tt,self%it,1) = qrad
    deallocate (qrad)
  end if                              
  call self%courant_condition
  if (io%verbose>2) &
    call self%check_density('before timestep')
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call self%check_density(' after timestep')
  end associate
  call trace%end (itimer)
END SUBROUTINE update

END MODULE rmhd_mod
