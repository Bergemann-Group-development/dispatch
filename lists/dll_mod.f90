!===============================================================================
!> Doubly linked list (DLL), carrying anything, as simply as possible
!===============================================================================
MODULE dll_mod
  USE trace_mod
  USE io_mod
  USE omp_lock_mod
  implicit none
  private
  type, public:: dll_node_t
    class(dll_node_t), pointer:: prev => null()
    class(dll_node_t), pointer:: next => null()
    class(*), pointer:: car
    integer:: id=0
  end type
  type, public:: dll_t
    class(dll_node_t), pointer:: head => null()
    class(dll_node_t), pointer:: tail => null()
    integer:: n=0
    character(len=32):: name
    type(lock_t), pointer:: lock => null()
    logical:: initialized = .false.
  contains
    procedure:: init
    procedure:: dealloc
    procedure:: append
    procedure:: prepend
    procedure:: find
    procedure:: find_id
    procedure:: insert_before
    procedure:: insert_after
    procedure:: remove
    procedure:: delete
    procedure:: pop
    procedure:: test
    procedure, private:: check_initialized
  end type
CONTAINS

!===============================================================================
!> Initialize a doubly-linked-list (DLL), with anonymous carry
!===============================================================================
SUBROUTINE init (self, name)
  class(dll_t):: self
  character(len=*), optional:: name
  !.............................................................................
  call trace%begin ('dll_t%init')
  call io%lock%set ('dll_t%init')
  nullify (self%head, self%tail)
  self%n = 0
  if (present(name)) then
    self%name = name
  else
    self%name = 'dll'
  end if
  if (.not.associated(self%lock)) then
    allocate (self%lock)
    call self%lock%init ('dll')
  end if
  self%initialized = .true.
  call io%lock%unset ('dll_t%init')
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Initialize a doubly-linked-list (DLL), with anonymous carry
!===============================================================================
SUBROUTINE dealloc (self)
  class(dll_t):: self
  !.............................................................................
  call trace%begin ('dll_t%dealloc')
  call io%lock%set ('dll_t%dealloc')
  if (associated(self%lock)) then
    call self%lock%destroy
    deallocate (self%lock)
  end if
  call io%lock%unset ('dll_t%dealloc')
  call trace%end ()
END SUBROUTINE dealloc

!===============================================================================
!> Add an item at the tail
!===============================================================================
SUBROUTINE append (self, new)
  class(dll_t):: self
  class(dll_node_t), pointer:: new
  !.............................................................................
  call self%check_initialized
  new%prev => self%tail
  nullify(new%next)
  if (.not.associated(self%head)) self%head => new
  if (associated(self%tail)) self%tail%next => new
  self%tail => new
  !$omp atomic
  self%n = self%n+1
END SUBROUTINE append

!===============================================================================
!> Add an item at the head
!===============================================================================
SUBROUTINE prepend (self, new)
  class(dll_t):: self
  class(dll_node_t), pointer:: new
  !.............................................................................
  call self%check_initialized
  if (.not.associated(self%tail)) self%tail => new
  new%next => self%head
  nullify (new%prev)
  self%head => new
  !$omp atomic
  self%n = self%n+1
END SUBROUTINE prepend

!===============================================================================
!> Find an item, by matching with a known member
!===============================================================================
SUBROUTINE find (self, old, prev, next)
  class(dll_t):: self
  class(dll_node_t), pointer:: prev, next
  class(*), pointer:: old
  !.............................................................................
  call self%check_initialized
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

!===============================================================================
!> Find an item, by matching with a known member
!===============================================================================
SUBROUTINE find_id (self, id, node)
  class(dll_t):: self
  class(dll_node_t), pointer:: node
  integer:: id
  !.............................................................................
  call self%check_initialized
  node => self%head
  do while (associated(node))
    if (node%id==id) then
      return
    end if
    node => node%next
  end do
END SUBROUTINE find_id

!===============================================================================
!> Insert a new item before an old one
!===============================================================================
SUBROUTINE insert_before (self, old, new)
  class(dll_t):: self
  class(dll_node_t), pointer:: old, new
  !.............................................................................
  call self%check_initialized
  new%prev => old%prev
  new%next => old
  if (associated(new%prev)) new%prev%next => new
  if (associated(new%next)) new%next%prev => new
  !$omp atomic
  self%n = self%n+1
END SUBROUTINE insert_before

!===============================================================================
!> Insert a new item after an old one
!===============================================================================
SUBROUTINE insert_after (self, old, new)
  class(dll_t):: self
  class(dll_node_t), pointer:: old, new
  !.............................................................................
  call self%check_initialized
  new%next => old%next
  old%next => new
  if (associated(new%prev)) new%prev%next => new
  if (associated(new%next)) new%next%prev => new
  !$omp atomic
  self%n = self%n+1
END SUBROUTINE insert_after

!===============================================================================
!> Remove an item from the DLL; note that the item is not deallocated!
!===============================================================================
SUBROUTINE remove (self, old)
  class(dll_t):: self
  class(dll_node_t), pointer:: old
  !.............................................................................
  call self%check_initialized
  if (associated(old%prev)) old%prev%next => old%next
  if (associated(old%next)) old%next%prev => old%prev
  if (associated(old, self%head)) then
    self%head => old%next
  end if
  if (associated(old, self%tail)) then
    self%tail => old%prev
  end if
  !$omp atomic
  self%n = self%n-1
END SUBROUTINE remove

!===============================================================================
!> Pop the last appended item off the list.  The list node and its carry should
!> be deallocated after use.
!===============================================================================
FUNCTION pop (self)
  class(dll_t):: self
  class(dll_node_t), pointer:: pop
  !.............................................................................
  call self%check_initialized
  pop => self%tail
  if (associated(pop%prev)) then
    self%tail => pop%prev
    nullify(self%tail%next)
  end if
  !$omp atomic
  self%n = self%n-1
END FUNCTION pop

!===============================================================================
!> Run a self test
!===============================================================================
SUBROUTINE test (self)
  class(dll_t):: self
  !.............................................................................
  integer:: n
  class(*), pointer:: car, save
  integer, pointer:: ip
  real, pointer:: rp
  class(dll_node_t), pointer:: prev, next
  integer:: i
  !-----------------------------------------------------------------------------
  ! Append five numbers
  !-----------------------------------------------------------------------------
  print *,'===================== Doubly linked list test ======================='
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

!===============================================================================
!> Delete a double-linked list, including carried content
!===============================================================================
SUBROUTINE delete (self)
  class(dll_t):: self
  !.............................................................................
  class (dll_node_t), pointer:: item, next
  !-----------------------------------------------------------------------------
  call self%check_initialized
  item => self%head
  do while (associated(item))
    next => item%next
    if (associated(item%car)) deallocate (item%car)
    deallocate (item)
    item => next
  end do
  nullify (self%head)
  nullify (self%tail)
  self%n = 0
  call self%init
END SUBROUTINE delete

!===============================================================================
!> Check that the dll has been initialized
!===============================================================================
SUBROUTINE check_initialized (self)
  class(dll_t):: self
  !-----------------------------------------------------------------------------
  if (.not. self%initialized) then
    call io%abort ('dll_t instance not initialized')
  end if
END SUBROUTINE check_initialized

END MODULE dll_mod
