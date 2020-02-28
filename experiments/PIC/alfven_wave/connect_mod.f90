!===============================================================================
!> Module holding anonymous pointers back to extras features
!===============================================================================
MODULE connect_mod
  implicit none
  private
  type, public:: connect_t
    !class(*), pointer:: trace_particles
    !class(*), pointer:: radiative_transfer
    class(*), pointer:: pic_patch
  end type
END MODULE connect_mod
