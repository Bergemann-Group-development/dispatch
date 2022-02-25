!===============================================================================
!> Interface on top of solvers and boundary conditions.  The init and update
!> entries are redundant in this case, but the bndries() entry is not.
!>
!> The task%bndries() procedure is called at the same time as guard zones are
!> loaded, and should do whatever is required to set the boundaries -- it is
!> not reqiured that it calls a framework procedure as part of this.  However,
!> a number of standard boundary conditions are available in the subsystem
!> that has standard_boundaries available for inheritence.
!>
!> A standard way to implement BCs is thus to create a boundaries_mod.f90, and
!> let it inherit standard_boundaries_t% and its methods.  Here, the entry
!> point %primary() is called from task%bndries(), and the actual BCs are
!> implented there.
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE solver_mod
  USE kinds_mod
  USE eos_mod
  USE benchmark_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
    procedure:: bndries
    procedure:: update
  end type
CONTAINS

!===============================================================================
!> Setup initial states in the patches
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer:: iostat
  real, save:: z0    = 0.2
  real, save:: z1    = 1.0
  real, save:: mu    = 0.0
  real, save:: grav  = 5.
  real, save:: kappa = 0.07
  real, save:: mpoly = 4.
  real, save:: q_rate = 0.0
  real, save:: sigma = 1.0
  real, save:: theta = 1.
  real, save:: stencil = 2.0
  real, save:: courant_v = 1.0
  real, save:: courant_c = 1.0
  real, save:: amplitude(2) = 0.0
  integer, save:: seed = 11
  logical, save:: kinematic = .true.
  logical, save:: compressive = .false.
  logical, save:: first_time=.true.
  namelist /experiment_params/ &
    z0, z1, theta, kappa, sigma, mpoly, theta, stencil, courant_v, courant_c, &
    kinematic, amplitude, grav, compressive, q_rate, seed
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init')
  call io%lock%set ('experiment_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, experiment_params)
    end if
    call io%namelist_warning ('experiment_params', iostat)
    write (io%namelist, experiment_params)
    benchmark%n_poly  = mpoly
    benchmark%q_rate  = q_rate
    benchmark%z0      = z0
    benchmark%z1      = z1
    benchmark%perturb = amplitude(1)
    first_time = .false.
  end if
  !---------------------------------------------------------------------------
  ! Compute and set viscosity and conduction parameters
  !---------------------------------------------------------------------------
  mu = (eos%gamma-1d0) * kappa / eos%gamma * sigma
  eos%gamma                     = 1.0 + 1.0/mpoly
  grav                          = (mpoly+1)/(z1+1)
  self%viscosity%mu             = mu
  self%conduction%kappa         = kappa
  self%conduction%stencil       = stencil
  self%conduction%kinematic     = kinematic
  self%conduction%courant       = courant_c
  self%viscosity%kinematic      = kinematic
  self%viscosity%compressive    = compressive
  self%viscosity%courant        = courant_v
  self%gravity%const_grav       = grav
  self%initial_conditions%seed  = seed
  self%initial_conditions%theta = theta
  self%initial_conditions%mpoly = mpoly
  self%initial_conditions%amplitude = amplitude
  call self%solver_t%init
  self%boundaries%theta         = theta
  self%boundaries%mpoly         = mpoly
  self%boundaries%g             = grav
  self%mhd                      = .false.
  call io%lock%unset ('experiment_t%init')
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Setup initial states in the patches
!===============================================================================
SUBROUTINE bndries (self)
  class(experiment_t):: self
  real(var_k), dimension(:,:,:,:), pointer, contiguous:: mem
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%bndries')
  call dbg (self, self%it, 'before bndries')
  mem => self%mem(:,:,:,:,self%it,1)
  call self%boundaries%primary (mem)
  call dbg (self, self%it, ' after bndries')
  call trace%end()
END SUBROUTINE bndries

!===============================================================================
!> Setup initial states in the patches
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  real(var_k), dimension(:,:,:,:), pointer, contiguous:: mem
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call dbg (self, self%it, 'before update')
  call self%solver_t%update
  call dbg (self, self%new, ' after update')
  !----------------------------------------------------------------------------
  ! This call is intentionally redundant, and is there only to check that
  !----------------------------------------------------------------------------
  mem => self%mem(:,:,:,:,self%new,1)
  call self%boundaries%primary (mem)
  call dbg (self, self%new, ' after boundaries')
  call trace%end ()
END SUBROUTINE update

SUBROUTINE dbg (self, it, label)
  class(experiment_t):: self
  integer:: it
  character(*):: label
  if (self%id==0) print 1,self%time, log(self%mem(4,4,2:6,self%idx%d ,it,1)), &
                                         self%mem(4,4,2:6,self%idx%pz,it,1)  &
                                        /self%mem(4,4,2:6,self%idx%d ,it,1), &
                         (eos%gamma-1d0)*self%mem(4,4,2:6,self%idx%e ,it,1)  &
                                        /self%mem(4,4,2:6,self%idx%d ,it,1), &
                                        label
  1 format(f10.6,1x,5f10.4,1x,5f8.4,1x,5f8.4,1x,a)
END SUBROUTINE dbg

END MODULE
