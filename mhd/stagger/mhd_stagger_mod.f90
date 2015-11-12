!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for MHD task.  Copy this to your experiment directory, or
!> to a specialized solver directory, and modify with relevant parameters and
!> procedures..
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE mhd_solver_mod
  USE patch_mod
  USE trace_mod
  USE mhd_mod
  implicit none
  type, extends(mhd_t):: mhd_solver_t
    contains
    procedure:: init       => init_stagger
    procedure:: update     => update_stagger
    procedure:: allocate   => allocate_stagger
    procedure:: deallocate => deallocate_stagger
    procedure:: print      => print_stagger
    procedure:: check      => check_stagger
    procedure:: send_to    => send_to_stagger
  end type
PRIVATE
PUBLIC mhd_solver_t
CONTAINS

!> =============================================================================
!> Initialize -- reading parameters etc
!> =============================================================================
SUBROUTINE init_stagger (self)
  class(mhd_solver_t):: self
  integer:: i
  !.............................................................................
  call trace_begin ('init_stagger')
  if (verbose >0) print*,'init_stagger:', self%id, self%time
  do i=1,3
    allocate (self%mesh(i))
    call self%mesh(i)%init
  end do
  call self%allocate
  call trace_end
END SUBROUTINE

!> =============================================================================
!> Take a local time step
!> =============================================================================
SUBROUTINE update_stagger (self)
  class(mhd_solver_t):: self
  !.............................................................................
  if (verbose >0) print*,'update_stagger:', self%id, self%time
  self%time = self%time + self%dtime
END SUBROUTINE

!> =============================================================================
!> Check if neighbors are ready, so this task can execute
!> =============================================================================
SUBROUTINE check_stagger (self)
  class(mhd_solver_t):: self
  !.............................................................................
  if (verbose >0) print*,'check_stagger:', self%id, self%time
END SUBROUTINE

!> =============================================================================
!> Send data to another rank -- will likely be handled generically
!> =============================================================================
SUBROUTINE send_to_stagger (self, rank)
  class(mhd_solver_t):: self
  integer:: rank
  !.............................................................................
  if (verbose >0) print*,'check_stagger:', self%id, self%time
END SUBROUTINE

!> =============================================================================
!> Any extra allocations needed
!> =============================================================================
SUBROUTINE allocate_stagger  (self)
  class(mhd_solver_t):: self
  integer:: m(3)
  !.............................................................................
  if (verbose >0) print*,'allocate_stagger:', self%id, self%time
  if (associated(self%mem)) deallocate (self%mem)
  m = self%mesh%gn
  allocate (self%mem(m(1),m(2),m(3),self%nvar,self%nt,2))
END SUBROUTINE

!> =============================================================================
!> Any extra deallocations needed
!> =============================================================================
SUBROUTINE deallocate_stagger  (self)
  class(mhd_solver_t):: self
  !.............................................................................
  if (verbose >0) print*,'deallocate_stagger:', self%id, self%time
  if (associated(self%mem)) deallocate (self%mem)
END SUBROUTINE

!> =============================================================================
!> Print statistics
!> =============================================================================
SUBROUTINE print_stagger  (self)
  class(mhd_solver_t):: self
  if (verbose >0) print*,'print_stagger:', self%id, self%time
END SUBROUTINE

END MODULE
