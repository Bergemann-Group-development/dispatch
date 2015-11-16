!> =============================================================================
!> $Id$
!> This type exists to be able to organize task lists by time.  It corresponds
!> 0-dimensional tasks, which just evolve in time
!> =============================================================================
MODULE item_mod
  USE bits_mod
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

SUBROUTINE void (self)
  class (item_t):: self
  print*, 'ERROR: item%void called'
END SUBROUTINE

SUBROUTINE init (self)
  class (item_t):: self
  id = id+1
  self%id = id
  print*, 'item%init:', self%id, self%time
END SUBROUTINE

SUBROUTINE print (self)
  class (item_t):: self
  print*, 'item%print:', self%id, self%time
END SUBROUTINE

SUBROUTINE update (self)
  class (item_t):: self
  self%time = self%time + self%dtime
  print*,'item%update:', self%id, self%time
END SUBROUTINE

!> =============================================================================
!> Set and clear status bits
!> =============================================================================
SUBROUTINE set (self, bit)
  class (item_t):: self
  integer:: bit
  !.............................................................................
  self%status = ior(self%status, bit)
END SUBROUTINE
SUBROUTINE clear (self, bit)
  class (item_t):: self
  integer:: bit
  !.............................................................................
  self%status = iand(self%status, not(bit))
END SUBROUTINE
LOGICAL FUNCTION is_set (self, bit)
  class (item_t):: self
  integer:: bit
  !.............................................................................
  is_set = (ior(self%status, bit) /= 0)
END FUNCTION
LOGICAL FUNCTION is_clear (self, bit)
  class (item_t):: self
  integer:: bit
  !.............................................................................
  is_clear = (iand(self%status, bit) == 0)
END FUNCTION

END MODULE item_mod
