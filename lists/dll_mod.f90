!> =============================================================================
!> Doubly linked list, carrying anything, as simply as possible
!> =============================================================================
MODULE dll_mod
  implicit none
  private
  type, public:: dll_node_t
    type(dll_node_t), pointer:: prev => null()
    type(dll_node_t), pointer:: next => null()
    class(*), pointer:: car
  end type
  type, public:: dll_t
    type(dll_node_t), pointer:: head => null()
    type(dll_node_t), pointer:: tail => null()
    integer:: n=0
  contains
    procedure:: init
    procedure:: append
    procedure:: prepend
    procedure:: find
    procedure:: insert_before
    procedure:: insert_after
    procedure:: remove
    procedure:: test
  end type
CONTAINS

!> =============================================================================
SUBROUTINE init (self)
  class(dll_t):: self
  !.............................................................................
  nullify (self%head, self%tail)
  !$omp atomic write
  self%n = 0
END SUBROUTINE init

!> =============================================================================
SUBROUTINE append (self, new)
  class(dll_t):: self
  class(dll_node_t), pointer:: new
  !.............................................................................
  new%prev => self%tail
  if (.not.associated(self%head)) self%head => new
  if (associated(self%tail)) self%tail%next => new
  self%tail => new
  !$omp atomic
  self%n = self%n+1
END SUBROUTINE append

!> =============================================================================
SUBROUTINE prepend (self, new)
  class(dll_t):: self
  class(dll_node_t), pointer:: new
  !.............................................................................
  if (.not.associated(self%tail)) self%tail => new
  new%next => self%head
  self%head => new
  !$omp atomic
  self%n = self%n+1
END SUBROUTINE prepend

!> =============================================================================
SUBROUTINE find (self, old, prev, next)
  class(dll_t):: self
  class(dll_node_t), pointer:: prev, next
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
SUBROUTINE insert_before (self, old, new)
  class(dll_t):: self
  class(dll_node_t), pointer:: old, new
  !.............................................................................
  new%prev => old%prev
  new%next => old
  if (associated(new%prev)) new%prev%next => new
  if (associated(new%next)) new%next%prev => new
  !$omp atomic
  self%n = self%n+1
END SUBROUTINE insert_before

!> =============================================================================
SUBROUTINE insert_after (self, old, new)
  class(dll_t):: self
  class(dll_node_t), pointer:: old, new
  !.............................................................................
  new%next => old%next
  old%next => new
  if (associated(new%prev)) new%prev%next => new
  if (associated(new%next)) new%next%prev => new
  !$omp atomic
  self%n = self%n+1
END SUBROUTINE insert_after

!> =============================================================================
SUBROUTINE remove (self, old)
  class(dll_t):: self
  class(dll_node_t), pointer:: old
  !.............................................................................
  if (associated(old%prev)) old%prev%next => old%next
  if (associated(old%next)) old%next%prev => old%prev
  !$omp atomic
  self%n = self%n-1
END SUBROUTINE remove

!> =============================================================================
SUBROUTINE test (self)
  class(dll_t):: self
  class(*), pointer:: car, save
  integer, pointer:: ip
  real, pointer:: rp
  class(dll_node_t), pointer:: prev, next
  integer:: i, n
  !-----------------------------------------------------------------------------
  ! Append five numbers
  !-----------------------------------------------------------------------------
  print *,'===================== Doubly linked list test ======================='
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

END MODULE dll_mod
