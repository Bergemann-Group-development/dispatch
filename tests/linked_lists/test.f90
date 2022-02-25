!===============================================================================
!> Program to test functionality and speed of anonymous linked lists for tasks
!===============================================================================

MODULE test_mod
  USE list_mod
  integer, save:: id=0
CONTAINS
SUBROUTINE make_task_at_time (link, time)
  USE mhd_mod
  implicit none
  class(link_t), pointer:: link
  class(mhd_t),   pointer:: task
  real:: time
  !.............................................................................
  allocate (link)
  allocate (task)
  link%item => task
  id = id+1
  task%id = id
  task%time = time
END SUBROUTINE
END MODULE

PROGRAM test
  USE list_mod
  USE mhd_mod
  USE mpi_timer
  USE test_mod
  implicit none
  type(mhd_t), pointer:: task
  class(list_t), pointer:: task_list
  class(link_t), pointer:: link
  integer:: i, n, m
  !> ---------------------------------------------------------------------------
  !> Test task list routines
  !> ---------------------------------------------------------------------------
  print*,'-------------------- testing task list routines -----------------'
  allocate (task_list)
  call task_list%init ('task_list')

  print*,'first item'
  call make_task_at_time (link,4.9)
  call task_list%append (link)
  call task_list%print

  print*,'new before'
  call make_task_at_time (link,0.9)
  call task_list%add_by_time (link)
  call task_list%print

  print*,'new inside'
  call make_task_at_time (link,1.9)
  call task_list%add_by_time (link)
  call task_list%print

  print*,'new after'
  call make_task_at_time (link,7.9)
  call task_list%add_by_time (link)
  call task_list%print

  call tic
  n = 1000
  do i=1,n
    allocate (task)
    call make_task_at_time (link,real(i))
    call task_list%append (link)
  end do
  call timer%toc ('     append @  1000', n)

  call tic
  m = 200
  do i=1,m
    call make_task_at_time (link,n-m+(i-0.5))
    call task_list%add_by_time (link)
  end do
  call timer%toc ('add_by_time @  1000', m)

  call tic
  n = 10000
  do i=1,n
    call make_task_at_time (link,real(i))
    call task_list%append (link)
  end do
  call timer%toc ('     append @ 10000', n)

  call tic
  m = 1000
  do i=1,m
    call make_task_at_time (link, n-m+(i-0.5))
    call task_list%add_by_time (link)
  end do
  call timer%toc ('add_by_time @ 10000', m)

  call task_list%reset
  deallocate (task_list)
END PROGRAM
