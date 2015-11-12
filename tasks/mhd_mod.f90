!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for MHD task.  Copy this to your experiment directory, or
!> to a specialized solver directory, and modify with relevant parameters and
!> procedures..
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE mhd_mod
  USE patch_mod
  USE trace_mod
  implicit none
  type, extends(patch_t):: mhd_t
    real:: courant                  ! Courant number
    contains
    procedure:: init       => init_mhd
    procedure:: update     => update_mhd
    procedure:: allocate   => allocate_mhd
    procedure:: deallocate => deallocate_mhd
    procedure:: print      => print_mhd
    procedure:: check      => check_mhd
    procedure:: send_to    => send_to_mhd
  end type
PRIVATE
PUBLIC mhd_t
CONTAINS

SUBROUTINE init_mhd (self)
  class(mhd_t):: self
  integer:: i
  !.............................................................................
  call trace_begin ('init_mhd')
  if (verbose >0) print*,'init_mhd:', self%id, self%time
  do i=1,3
    allocate (self%mesh(i))
    call self%mesh(i)%init
  end do
  call self%allocate
  call trace_end
END SUBROUTINE

SUBROUTINE update_mhd (self)
  class(mhd_t):: self
  !.............................................................................
  if (verbose >0) print*,'update_mhd:', self%id, self%time
  self%time = self%time + self%dtime
END SUBROUTINE

SUBROUTINE check_mhd (self)
  class(mhd_t):: self
  !.............................................................................
  if (verbose >0) print*,'check_mhd:', self%id, self%time
END SUBROUTINE

SUBROUTINE send_to_mhd (self, rank)
  class(mhd_t):: self
  integer:: rank
  !.............................................................................
  if (verbose >0) print*,'check_mhd:', self%id, self%time
END SUBROUTINE

SUBROUTINE allocate_mhd  (self)
  class(mhd_t):: self
  integer:: m(3)
  !.............................................................................
  if (verbose >0) print*,'allocate_mhd:', self%id, self%time
  if (associated(self%mem)) deallocate (self%mem)
  m = self%mesh%gn
  allocate (self%mem(m(1),m(2),m(3),self%nvar,self%nt,2))
END SUBROUTINE

SUBROUTINE deallocate_mhd  (self)
  class(mhd_t):: self
  !.............................................................................
  if (verbose >0) print*,'deallocate_mhd:', self%id, self%time
  if (associated(self%mem)) deallocate (self%mem)
END SUBROUTINE

SUBROUTINE print_mhd  (self)
  class(mhd_t):: self
  if (verbose >0) print*,'print_mhd:', self%id, self%time
END SUBROUTINE

END MODULE
