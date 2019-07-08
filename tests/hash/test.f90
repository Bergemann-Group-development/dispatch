!===============================================================================
!> Test of storing anonymous pointers in a hash tale
!===============================================================================
PROGRAM test
  USE hash_mod
  USE omp_timer_mod
  implicit none
  integer, parameter:: n=1000
  type(hash_table_t):: ht
  integer(8):: key(1:2)
  class(*), pointer:: ptr
  integer, target:: value
  integer:: i
  integer, target:: a(n)
  real(8):: start
  logical:: ok
  !-----------------------------------------------------------------------------
  ! Make a table, and set the hash table ptr values
  !-----------------------------------------------------------------------------
  call ht%init (n)
  do i=1,n
    key = [i,i]
    a(i) = i
    ptr => a(i)
    call ht%set (key, ptr)
  end do
  !-----------------------------------------------------------------------------
  ! Get the ptr values from the hash table
  !-----------------------------------------------------------------------------
  ok = .true.
  do i=1,n
    key = [i,i]
    call ht%get (key, ptr)
    select type (ptr)
    type is (integer)
    if (ptr /= i) then
      ok = .false.
      print *, i, ptr
    end if
    end select
  end do
  if (ok) then
    print *, '1000 values retrieved OK'
  else
    print *, '1000 values NOT retrieved OK'
  end if
  !-----------------------------------------------------------------------------
  ! Time the retrieval (no print)
  !-----------------------------------------------------------------------------
  start = wallclock()
  do i=1,n
    key = [i,i]
    call ht%get (key, ptr)
    select type (ptr)
    type is (integer)
    end select
  end do
  print '(" time to get 1000 results:",6p,f6.2," mus")', wallclock()-start
END PROGRAM
