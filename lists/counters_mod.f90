!> =============================================================================
!> Help keep track of when all patches have passed some counter, by decrementing
!> a counter, from start to 0.  Typical use:
!>
!> if (self%time > some_time) then
!>   count = counters%count (id, io%ntask)
!>   if (count==0) then
!>     print *, 'all tasks have passed time =', some_time
!>     call counters%remove (id)
!>   end if
!> end if
!> =============================================================================
MODULE counters_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! counter data type, to keep track of output
  !-----------------------------------------------------------------------------
  type:: counter_t
    type(counter_t), pointer:: prev=>null(), next=>null()
    integer:: id=0, count=0, start=0
  end type
  type, public:: counters_t
    type(counter_t), pointer:: head=>null(), tail=>null()
    integer:: n=0
    integer:: id=0
  contains
    procedure:: decrement ! (id, start)
    procedure:: reset     ! (id)
  end type
  type(counters_t), public:: counters
CONTAINS

!> =============================================================================
!> Find the counter instance with a given id and decrement counter.  If a
!> counter with that id is not in the list, create one, and decrement it
!> =============================================================================
FUNCTION decrement (self, id, start) RESULT (count)
  integer:: count
  class(counters_t):: self
  type(counter_t), pointer:: counter
  integer:: id, start
  call trace%begin ('counters_t%decrement')
  !-----------------------------------------------------------------------------
  ! Search for the counter instance, return a new one if not found
  !-----------------------------------------------------------------------------
  !$omp critical (counter_cr)
  counter => find (self, id, start)
  !-----------------------------------------------------------------------------
  ! Decrement to zero
  !-----------------------------------------------------------------------------
  if (counter%count>0) then
    counter%count = counter%count-1    
  end if
  count = counter%count
  !$omp end critical (counter_cr)
  call trace %end
END FUNCTION decrement

!> =============================================================================
!> Find a counter with the given id, or create one
!> =============================================================================
FUNCTION find (self, id, start) RESULT (counter)
  class(counters_t):: self
  integer:: id, start
  optional:: start
  type(counter_t), pointer:: counter
  !-----------------------------------------------------------------------------
  call trace%begin ('counters_t%find')
  counter => self%head
  do while (associated(counter))
    if (counter%id==id) then
      exit
    end if
    counter => counter%next
  end do
  !-----------------------------------------------------------------------------
  ! If a counter instance does not exist, append a new one to the list
  !-----------------------------------------------------------------------------
  if (.not.associated(counter)) then
    allocate (counter)
    !---------------------------------------------------------------------------
    ! If id is zero on entry, give it a unique value, otherwise use the given id
    !---------------------------------------------------------------------------
    if (id<=0) then
      self%id = self%id + 1
      id = self%id
    else
      self%id = id
    end if
    call io%assert (present(start), &
      'counters_t%find: start parameter must be present')
    counter%id    = self%id
    counter%start = start
    counter%count = start
    call append (self, counter)
  end if
  call trace%end
END FUNCTION find

!> =============================================================================
!> Append a new counter record counter to the list 
!> =============================================================================
SUBROUTINE append (self, counter)
  class(counters_t):: self
  type(counter_t), pointer:: counter
  !-----------------------------------------------------------------------------
  call trace%begin ('counters_t%append')
  counter%prev => self%tail
  if (associated(self%tail)) then
    self%tail%next => counter
  end if
  self%tail => counter
  self%n = self%n+1
  if (.not.associated(self%head)) self%head => counter
  call trace%end
END SUBROUTINE append

!> =============================================================================
!> Remove the counter instance from the list
!> =============================================================================
SUBROUTINE reset (self, id)
  class(counters_t):: self
  integer:: id
  type(counter_t), pointer:: counter
  !-----------------------------------------------------------------------------
  call trace%begin ('counters_t%remove')
  !$omp critical (counter_cr)
  counter => find (self, id)
  if (associated(counter%prev)) then
    counter%prev%next => counter%next
  else
    counters%head => counter%next
  end if
  if (associated(counter%next)) then
    counter%next%prev => counter%prev
  else
    counters%tail => counter%prev
  end if
  deallocate(counter)
  nullify (counter)
  counters%n = counters%n-1
  !$omp end critical (counter_cr)
  call trace%end
END SUBROUTINE reset

END MODULE counters_mod
