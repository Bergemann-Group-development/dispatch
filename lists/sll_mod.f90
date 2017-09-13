!> =============================================================================
!> Singly link list module, carrying any kind of 
!> =============================================================================
MODULE sll_mod
  implicit none
  private
  type, public:: sll_node_t
    type(sll_node_t), pointer:: next => null()
    class(*), pointer:: car
  end type
  type, public:: sll_t
    type(sll_node_t), pointer:: head => null()
    type(sll_node_t), pointer:: tail => null()
    integer:: n=0
  contains
    procedure:: append
    procedure:: prepend
    procedure:: insert_after
    procedure:: remove
    procedure:: find
    procedure:: test
  end type
CONTAINS

!> =============================================================================
SUBROUTINE append (self, new)
  class(sll_t):: self
  class(sll_node_t), pointer:: new
  !.............................................................................
  if (.not.associated(self%head)) self%head => new
  if (associated(self%tail)) self%tail%next => new
  self%tail => new
  !$omp atomic
  self%n = self%n+1
END SUBROUTINE append

!> =============================================================================
SUBROUTINE prepend (self, new)
  class(sll_t):: self
  class(sll_node_t), pointer:: new
  !.............................................................................
  new%next => self%head
  self%head => new
  !$omp atomic
  self%n = self%n+1
END SUBROUTINE prepend

!> =============================================================================
SUBROUTINE find (self, old, prev, next)
  class(sll_t):: self
  class(sll_node_t), pointer:: prev, next
  class(*), pointer:: old
  !.............................................................................
  next => self%head
  nullify (prev)
  do while (associated(next))
    if (associated(next%car, old)) then
      return
    end if
    prev => next
    next => next%next
  end do
  nullify (prev, next)
END SUBROUTINE find

!> =============================================================================
SUBROUTINE insert_after (self, old, new)
  class(sll_t):: self
  class(sll_node_t), pointer:: old, new
  !.............................................................................
  new%next => old%next
  old%next => new
  !$omp atomic
  self%n = self%n+1
END SUBROUTINE insert_after

!> =============================================================================
SUBROUTINE remove (self, prev)
  class(sll_t):: self
  class(sll_node_t), pointer:: prev
  !.............................................................................
  if (associated(prev)) then
    prev%next => prev%next%next
    !$omp atomic
    self%n = self%n-1
  end if
END SUBROUTINE remove

!> =============================================================================
SUBROUTINE test (self)
  class(sll_t):: self
  class(*), pointer:: car, save
  integer, pointer:: ip
  real, pointer:: rp
  class(sll_node_t), pointer:: prev, next
  integer:: i, n
  !-----------------------------------------------------------------------------
  ! Append five numbers
  !-----------------------------------------------------------------------------
  print *,'===================== Singly linked list test ======================='
  n = 5
  print *, 'fill:', n
  do i=1,n
    if (i == 2) then
      allocate (next, ip)
      ip = i
      next%car => ip
    else
      allocate (next, rp)
      rp = i
      next%car => rp
    end if
    if (i == 3) save => next%car
    call self%append (next)
  end do
  !-----------------------------------------------------------------------------
  ! List them
  !-----------------------------------------------------------------------------
  print *, 'list:', self%n
  next => self%head
  do while (associated(next))
    call print (next%car)
    next => next%next
  end do
  !-----------------------------------------------------------------------------
  ! Find number 3
  !-----------------------------------------------------------------------------
  print *, 'find:'
  call print (save)
  call self%find (save, prev, next)
  call print (next%car)
  print *, 'previous:'
  call print (prev%car)
contains
  subroutine print (car)
    class(*), pointer:: car
    select type (car)
    type is (integer)
      print *, car
    type is (real)
      print *, car
    end select
  end subroutine print
END SUBROUTINE test

END MODULE sll_mod
