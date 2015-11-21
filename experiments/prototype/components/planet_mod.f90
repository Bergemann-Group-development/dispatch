MODULE planet_mod
  USE particle_mod
  implicit none
  private
  type, public, extends(particle_t):: planet_t
  end type
END MODULE planet_mod
