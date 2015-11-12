!> =============================================================================
!> $Id$
!> =============================================================================
MODULE item_mod
  type, abstract:: item_t
    real(8):: time = 0d0
  contains
    procedure(one), deferred:: init
    procedure(one), deferred:: execute
    procedure(one), deferred:: update
    procedure(one), deferred:: allocate
    procedure(one), deferred:: deallocate
    procedure(one), deferred:: print
  end type
  abstract interface
    subroutine one (self)
    import item_t
    class(item_t):: self
    end subroutine
  end interface
PRIVATE
PUBLIC item_t
END MODULE item_mod
