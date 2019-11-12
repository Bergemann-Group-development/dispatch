!===============================================================================
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
!===============================================================================
MODULE counters_mod
  USE io_mod
  USE trace_mod
  USE omp_lock_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! counter data type, to keep track of output
  !-----------------------------------------------------------------------------
  type:: counter_t
    type(counter_t), pointer:: prev=>null(), next=>null()
    integer:: id=0, count=0, start=0
    type(lock_t):: lock
  end type
  type, public:: counters_t
    type(counter_t), pointer:: head=>null(), tail=>null()
    integer:: n=0
    integer:: id=0
    type(lock_t):: lock
  contains
    procedure:: init
    procedure:: decrement ! (id, start)
    procedure:: increment ! (id, start)
    procedure:: update    ! (id, incr)
    procedure:: remove    ! (id)
  end type
  integer:: verbose=0
  type(counters_t), public:: counters
CONTAINS

!===============================================================================
!> Allow changing verbosity
!===============================================================================
SUBROUTINE init(self)
  class(counters_t):: self
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /counters_params/ verbose
  !-----------------------------------------------------------------------------
  if (first_time) then
    call trace%begin ('counters_t%init')
    !$omp critical (input_cr)
    if (first_time) then
      self%id = 0
      rewind (io_unit%input)
      read (io_unit%input, counters_params, iostat=iostat)
      write (io_unit%output, counters_params)
      call self%lock%init ('coun')
      first_time = .false.
    end if
    !$omp end critical (input_cr)
    call trace%end()
  end if
END SUBROUTINE init

!===============================================================================
!> Find the counter instance with a given id and decrement counter.  If a
!> counter with that id is not in the list, create one, and decrement it
!===============================================================================
FUNCTION decrement (self, id, start) RESULT (count)
  integer:: count
  class(counters_t):: self
  type(counter_t), pointer:: counter
  integer:: id, start
  !-----------------------------------------------------------------------------
  call trace%begin ('counters_t%decrement')
  if (self%id == -1) call self%init
  !-----------------------------------------------------------------------------
  ! Search for the counter instance, return a new one if not found
  !-----------------------------------------------------------------------------
  call self%lock%set ('counters_t%decrement')
  counter => find (self, id, start)
  call self%lock%unset ('counters_t%decrement')
  !-----------------------------------------------------------------------------
  ! Decrement to zero
  !-----------------------------------------------------------------------------
  if (counter%count>0) then
    !$omp atomic capture
    counter%count = counter%count-1    
    count = counter%count
    !$omp end atomic
  else
    count = counter%count
  end if
  if (verbose > 0) &
    write (io_unit%output,*) 'counter_t%decrement: id, count =', counter%id, count
  call trace%end()
END FUNCTION decrement

!===============================================================================
!> Find the counter instance with a given id and decrement counter.  If a
!> counter with that id is not in the list, create one, with count=start.
!===============================================================================
FUNCTION increment (self, id, start) RESULT (count)
  integer:: count
  class(counters_t):: self
  type(counter_t), pointer:: counter
  integer:: id, start
  !-----------------------------------------------------------------------------
  call trace%begin ('counters_t%increment')
  if (self%id == -1) call self%init
  !-----------------------------------------------------------------------------
  ! Search for the counter instance, return a new one if not found
  !-----------------------------------------------------------------------------
  call self%lock%set ('counters_t%increment')
  counter => find (self, id, start)
  call self%lock%unset ('counters_t%increment')
  !-----------------------------------------------------------------------------
  ! Increment away from zero
  !-----------------------------------------------------------------------------
  !$omp atomic capture
  counter%count = counter%count+1    
  count = counter%count
  !$omp end atomic
  if (verbose > 0) &
    write (io_unit%output,*) 'counter_t%increment: id, count =', counter%id, count
  call trace%end()
END FUNCTION increment

!===============================================================================
!> Find the counter instance with a given id and update counter.  If a counter
!> with the given ID does not exist, create one with count=start initially.
!===============================================================================
SUBROUTINE update (self, id, start, incr, count)
  class(counters_t):: self
  integer:: id, start, incr, count
  !.............................................................................
  type(counter_t), pointer:: counter
  !-----------------------------------------------------------------------------
  call trace%begin ('counters_t%update')
  call self%lock%set ('counters_t%update')
  counter => find (self, id, start)
  call self%lock%unset ('counters_t%update')
  !$omp atomic capture
  counter%count = counter%count + incr
  count = counter%count
  !$omp end atomic
  if (verbose > 0) &
    write (stdout,*) 'counter_t%update: id, start, count =', id, start, count
  call trace%end()
END SUBROUTINE update

!===============================================================================
!> Find a counter with the given id, or create one
!===============================================================================
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
    call counter%lock%init ('counter')
    call counter%lock%set ('init')
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
      'counters_t%find: start parameter must be present on new counters')
    counter%id    = self%id
    counter%start = start
    counter%count = start
    call counter%lock%set ('init')
    call append (self, counter)
    if (verbose > 0) then
      write (stdout,*) 'counters_t%find: new id, start =', counter%id, start
    end if
  end if
  call trace%end()
END FUNCTION find

!===============================================================================
!> Append a new counter record counter to the list 
!===============================================================================
SUBROUTINE append (self, counter)
  class(counters_t):: self
  type(counter_t), pointer:: counter
  !-----------------------------------------------------------------------------
  call trace%begin ('counters_t%append')
  call self%lock%set ('append')
  counter%prev => self%tail
  if (associated(self%tail)) then
    self%tail%next => counter
  end if
  self%tail => counter
  self%n = self%n+1
  if (.not.associated(self%head)) self%head => counter
  call self%lock%unset ('append')
  call trace%end()
END SUBROUTINE append

!===============================================================================
!> Remove the counter instance from the list
!===============================================================================
SUBROUTINE remove (self, id)
  class(counters_t):: self
  integer:: id
  type(counter_t), pointer:: counter
  !-----------------------------------------------------------------------------
  call trace%begin ('counters_t%remove')
  call self%lock%set ('counters_t%remove')
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
  call self%lock%unset ('counters_t%remove')
  call trace%end()
END SUBROUTINE remove

END MODULE counters_mod
