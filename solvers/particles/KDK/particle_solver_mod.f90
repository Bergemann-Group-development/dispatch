!===============================================================================
!> Kick-Drift-Kick (KDK) N-body solver
!===============================================================================
MODULE particle_solver_mod
  USE io_unit_mod
  USE trace_mod
  USE particle_mod
  USE particle_list_mod
  USE dll_mod
  implicit none
  private
  type, public, extends(particle_list_t):: particle_solver_t
  contains
    procedure:: init
    procedure:: update
    procedure:: courant_time
  end type
  real(8), save:: courant=0.5
  type(particle_solver_t), target, public:: particle_solver
CONTAINS

!===============================================================================
!> Initialize a particle list
!===============================================================================
SUBROUTINE init (self, name)
  class(particle_solver_t):: self
  character(len=*), optional:: name
  logical, save:: first_time=.true.
  namelist /kdk_params/ courant
  !-----------------------------------------------------------------------------
  call self%particle_list_t%init
  if (first_time) then
    !$omp critical (input_cr)
    if (first_time) then
      first_time = .false.
      rewind (io_unit%input)
      read (io_unit%input, kdk_params)
      write (io_unit%output, kdk_params)
    end if
    !$omp end critical (input_cr)
  end if
END SUBROUTINE init

!===============================================================================
!> KDK update
!===============================================================================
SUBROUTINE update (self, p)
  class(particle_solver_t):: self
  class(particle_t), pointer:: p
  real(8):: dt, a(3)
  integer:: new, nt
  !.............................................................................
  new = 1 + mod(p%it,self%nt)
  associate (v=>p%v(:,new), r=>p%r(:,new))
  dt = self%courant_time (p)                    ! estimate dt
  a = self%force(p)/p%mass                      ! acceleration
  v = v + 0.5_8*a*dt                            ! kick 1
  r = r + v*dt                                  ! drift
  a = self%force(p)/p%mass                      ! acceleration
  v = v + 0.5_8*a*dt                            ! kick 2
  end associate
  p%time = p%time+dt
  p%t(new) = p%time
  p%it = new
  nt = self%nt
  p%iit(1:nt-1) = p%iit(2:nt)
  p%iit(nt) = new
END SUBROUTINE update

!===============================================================================
!> Courant conditions:  Estimate, as inexpensively as possible, an approximately
!> reflexive time step. If an estimate of the acceleration results in it not
!> being important for the timestep size, then choose a fraction of the cell
!> travel time (whatever that means in this context)
!===============================================================================
FUNCTION courant_time (self, p) RESULT (dt)
  class(particle_solver_t):: self
  class(particle_t), pointer:: p
  real(8):: dt
  !.............................................................................
  real(8):: a(3), ac, vc, dt1, dt2
  integer:: new, i0, i1
  !-----------------------------------------------------------------------------
  new = 1 + mod(p%it,self%nt)
  i1 = p%iit(self%nt-1)
  i0 = p%iit(self%nt-2)
  associate (v=>p%v(:,i1))
  a = (p%v(:,i1)-p%v(:,i0))/(p%t(i1)-p%t(i0))
  ac = 1.0/sqrt(sum(a**2)+tiny(1.0_8))
  vc = 1.0/sqrt(sum(v**2)+tiny(1.0_8))
  dt1 = p%ds*vc*courant
  dt2 = p%ds*ac/vc
  if (dt2 < dt1) then
    a = self%force(p)/p%mass
    ac = 1.0/sqrt(sum(a**2)+tiny(1.0_8))
    dt2 = p%ds*ac/vc*courant
  end if
  dt = min(dt1,dt2)
  end associate
END FUNCTION courant_time

END MODULE particle_solver_mod
