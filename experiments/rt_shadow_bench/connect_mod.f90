!===============================================================================
!> Module holding anonymous pointers back to extras features
!===============================================================================
MODULE connect_mod
  implicit none
  private
  type, public:: connect_t
    !class(*), pointer:: trace_particles
    class(*), pointer:: rt
  end type
END MODULE connect_mod
