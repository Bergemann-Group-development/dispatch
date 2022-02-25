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
  call self%solver_t%init
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
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: dedt, pgas, s
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: dddt, dsdt, d
  integer, save :: itimer=0
  real:: u_max
  !---------------------------------------------------------------------
  call trace%begin ('rmhd_t%update', itimer=itimer)
  call self%mhd_t%pre_update
  call self%output

  if (eos%do_rt.and.associated(self%rt)) then
    call self%rt%update
    call self%pde
      associate(l => self%li, u => self%ui)
      allocate(dedt(self%n(1), self%n(2), self%n(3)))
      allocate(pgas(self%gn(1), self%gn(2), self%gn(3)))
      dedt=0.0
      call self%rt%get_q (self%time, dedt)                                  ! interpolate in time and space
      d    => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%d,self%it,1)
      s    => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%s,self%it,1)
      call eos%lookup (d=d, s=s, pg=pgas)
      dsdt => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%s,self%it,2)
      dsdt = dsdt + dedt*d/pgas(l(1):u(1),l(2):u(2),l(3):u(3))
      !-----------------------------------------------------------------------------
      ! Time step limitation for rate of change dS/dt = (ds/dt-(s/d)*dddt)/d
      !-----------------------------------------------------------------------------
      dddt => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%d,self%it,2)
      if (io%debug(1)) then
        print *,'mk: cdtd', self%cdtd
        print *,'mk:    d', minval(d), maxval(d)
        print *,'mk: dedt', minval(dedt), maxval(dedt)
        print *,'mk: pgas', minval(pgas(l(1):u(1),l(2):u(2),l(3):u(3))), &
                            maxval(pgas(l(1):u(1),l(2):u(2),l(3):u(3)))
        print *,'mk: dsdt', minval(dsdt), maxval(dsdt)
        print *,'mk: dddt', minval(dddt), maxval(dddt)
        print *,'mk:    s', minval(s), maxval(s)
      end if
      if (io%verbose>0) then
        if (self%id == 100) print *, self%id, 'minmax dsdt', minval(dsdt), maxval(dsdt)
        if (self%id == 100) print *, 'minmax S', &
         minval(self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%s,self%it,1)), &
         maxval(self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%s,self%it,1))
      end if
      deallocate(dedt, pgas)
      d    => self%mem(:,:,:,self%idx%d,self%it,1)
      dsdt => self%mem(:,:,:,self%idx%s,self%it,2)
      s    => self%mem(:,:,:,self%idx%s,self%it,1)
      dddt => self%mem(:,:,:,self%idx%d,self%it,2)
      u_max = maxval(self%ds)*self%fmaxval(abs((dsdt-s*dddt/d)/d))*self%courant/self%cdtd
      self%u_max = max(self%u_max,u_max)
    if (io%debug(1)) print 1, self%id, 'u_max(RT)', u_max, self%u_max
    1 format(i6,3x,a,1p,2g12.3)
    end associate
  else
    call self%pde
  end if                         
  call self%courant_condition
  if (io%verbose>2) &
    call self%check_density('before timestep')
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh, self%lock)
  call self%check_density(' after timestep')
  !end if
  call self%counter_update
  call self%mhd_t%post_update
  call trace%end (itimer)
END SUBROUTINE update

END MODULE rmhd_mod
