MODULE particle_mod
  implicit none
  type, public:: particle_t
    type(particle_t), pointer:: next
    real:: r(3)                         ! position
    real:: v(3)                         ! velocity
    real:: m                            ! mass
    real(8):: t=0d0                     ! time increment
  end type
END MODULE particle_mod
