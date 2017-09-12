
MODULE particle_mod
  USE trace_mod
  implicit none
  private
  type, public:: particle_t
    integer:: id=0
    real, dimension(3):: r=0.0, v=0.0, a=0.0
    real:: t=0.0, dt=1.0, m=1.0
  contains
    procedure:: kick
    procedure:: drift
  end type
CONTAINS

SUBROUTINE kick (self)
  class(particle_t):: self
  self%v = self%v + 0.5*self%dt*self%a
END SUBROUTINE kick

SUBROUTINE drift (self)
  class(particle_t):: self
  self%r = self%r + self%dt*self%v
  self%t = self%t + self%dt
END SUBROUTINE drift

END MODULE particle_mod
