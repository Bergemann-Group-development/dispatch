!===============================================================================
!> Test demonstrating how to loop in parallel over the list
!===============================================================================
PROGRAM test
  USE hash_mod
  USE omp_timer_mod
  implicit none
  type(hash_table_t):: ht
  integer(8):: key(1:2)
  integer:: value
  !-----------------------------------------------------------------------------
  call ht%init (1000)
  key = [111,2]
  value = 118
  call ht%set (key, value)
  value = 0
  value = ht%get (key)
  print *, value
END PROGRAM
