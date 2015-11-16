!> =============================================================================
!> $Id$
!> This type exists to be able to organize task lists by time.  It corresponds
!> 0-dimensional tasks, which just evolve in time
!> =============================================================================
MODULE item_mod
  USE bits_mod
  USE trace_mod
  USE io_mod
  implicit none
  private
  type, public:: item_t
    real(8):: time = 0d0
    real(8):: dtime = 0.1d0
    integer:: status = 0
    integer:: id = 0
  contains
    procedure:: init
    procedure:: update
    procedure:: set
    procedure:: clear
    procedure:: print
    procedure:: is_set
    procedure:: is_clear
    procedure:: check => void
    procedure:: allocate => void
    procedure:: deallocate => void
    procedure:: execute => void
  end type
  integer:: id=0
CONTAINS

!> =============================================================================
SUBROUTINE void (self)
  class (item_t):: self
  print*, 'ERROR: item%void called'
END SUBROUTINE
!> =============================================================================
SUBROUTINE init (self)
  class (item_t):: self
  !.............................................................................
  call trace_begin ('item_t%init')
  id = id+1
  self%id = id
  call trace_end
END SUBROUTINE
!> =============================================================================
SUBROUTINE print (self)
  class (item_t):: self
  if (io%master) print*, 'item%print:', self%id, self%time
END SUBROUTINE
!> =============================================================================
SUBROUTINE update (self)
  class (item_t):: self
  !.............................................................................
  call trace_begin ('item_t%update')
  self%time = self%time + self%dtime
  call trace_end
END SUBROUTINE

!> =============================================================================
!> Set and clear status bits
!> =============================================================================
SUBROUTINE set (self, bit)
  class (item_t):: self
  integer:: bit
  !.............................................................................
  call trace_begin ('item_t%set')
  self%status = ior(self%status, bit)
  call trace_end
END SUBROUTINE
SUBROUTINE clear (self, bit)
  class (item_t):: self
  integer:: bit
  !.............................................................................
  call trace_begin ('item_t%clear')
  self%status = iand(self%status, not(bit))
  call trace_end
END SUBROUTINE
LOGICAL FUNCTION is_set (self, bit)
  class (item_t):: self
  integer:: bit
  !.............................................................................
  call trace_begin ('item_t%is_set')
  is_set = (ior(self%status, bit) /= 0)
  call trace_end
END FUNCTION
LOGICAL FUNCTION is_clear (self, bit)
  class (item_t):: self
  integer:: bit
  !.............................................................................
  call trace_begin ('item_t%is_clear')
  is_clear = (iand(self%status, bit) == 0)
  call trace_end
END FUNCTION

END MODULE item_mod
