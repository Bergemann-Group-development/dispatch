!> =============================================================================
!> Program to test functionality and speed of anonymous linked lists for tasks
!> =============================================================================

MODULE test_mod
  integer, save:: id=0
CONTAINS
SUBROUTINE make_task_at_time (task, time)
  USE task_mod
  implicit none
  class(task_t), pointer:: task
  real:: time
  !.............................................................................
  allocate (task)
  id = id+1
  task%id = id
  task%time = time
END SUBROUTINE
END MODULE

PROGRAM test
  USE list_mod
  USE task_mod
  USE timer_mod
  USE test_mod
  implicit none
  class(task_t), pointer:: task
  type(list_t), pointer:: task_list
  type(link_t), pointer:: link
  class(*), pointer:: item
  integer:: i, n, m
  !> ---------------------------------------------------------------------------
  !> Test task list routines
  !> ---------------------------------------------------------------------------
  print*,'-------------------- testing task list routines -----------------'
  allocate (task_list)
  call task_list%init ('task_list')

  print*,'first item'
  call make_task_at_time (task,4.9)
  call task_list%append (task)
  call task_list%print

  print*,'new before'
  call make_task_at_time (task,0.9)
  call task_list%add_by_time (task)
  call task_list%print

  print*,'new inside'
  call make_task_at_time (task,1.9)
  call task_list%add_by_time (task)
  call task_list%print

  print*,'new after'
  call make_task_at_time (task,7.9)
  call task_list%add_by_time (task)
  call task_list%print

  call tic
  n = 1000
  do i=1,n
    allocate (task_t::task)
    task%time = i
    call task_list%append (task)
  end do
  call toc ('     append @  1000', n)

  call tic
  m = 200
  do i=1,m
    allocate (task_t::task)
    task%time = n-m+(i-0.5)
    call task_list%add_by_time (task)
  end do
  call toc ('add_by_time @  1000', m)

  call tic
  n = 10000
  do i=1,n
    allocate (task_t::task)
    task%time = i
    call task_list%append (task)
  end do
  call toc ('     append @ 10000', n)

  call tic
  m = 1000
  do i=1,m
    allocate (task_t::task)
    task%time = n-m+(i-0.5)
    call task_list%add_by_time (task)
  end do
  call toc ('add_by_time @ 10000', m)

  call task_list%reset
  deallocate (task_list)
END PROGRAM
