!> =============================================================================
!> Program to test functionality and speed of anonymous linked lists for tasks
!> =============================================================================

MODULE test_mod
CONTAINS
SUBROUTINE add_task_link_at (link, time)
  USE list_mod
  USE task_mod
  USE timer_mod
  implicit none
  class(task_t), pointer:: task
  type(link_t),  pointer:: link
  real:: time
  !.............................................................................
  allocate (link)
  allocate (task)
  link%item => task
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
  call add_task_link_at (link,4.9)
  call task_list%append (link)
  call task_list%print

  print*,'new before'
  call add_task_link_at (link,0.9)
  call task_list%add_by_time (link)
  call task_list%print

  print*,'new inside'
  call add_task_link_at (link,1.9)
  call task_list%add_by_time (link)
  call task_list%print

  print*,'new after'
  call add_task_link_at (link,7.9)
  call task_list%add_by_time (link)
  call task_list%print

  call tic
  n = 1000
  do i=1,n
    allocate (link)
    allocate (task_t::task)
    link%item => task
    task%time = i
    call task_list%append (link)
  end do
  call toc ('     append @  1000', n)

  call tic
  m = 200
  do i=1,m
    allocate (link)
    allocate (task_t::task)
    link%item => task
    task%time = n-m+(i-0.5)
    call task_list%add_by_time (link)
  end do
  call toc ('add_by_time @  1000', m)

  call tic
  n = 10000
  do i=1,n
    allocate (link)
    allocate (task_t::task)
    link%item => task
    task%time = i
    call task_list%append (link)
  end do
  call toc ('     append @ 10000', n)

  call tic
  m = 1000
  do i=1,m
    allocate (link)
    allocate (task_t::task)
    link%item => task
    task%time = n-m+(i-0.5)
    call task_list%add_by_time (link)
  end do
  call toc ('add_by_time @ 10000', m)

END PROGRAM
